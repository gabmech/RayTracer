// Basic includes to get this file to work.  
// Basic includes to get this file to work.  
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <deque>
#include <stack>

#include "Transform.h"

using namespace std ; 
#include "FreeImage.h"
#include "variables.h"
#include "readfile.h"
#include "display.h"
#include "Geometry.h"


Intersection intersect(Ray);
void rayTrace(vec3);
vec3 findColor(Intersection, int, vec3);
int isPointNotInShadow(int lightIndex, int indexOfMinT);
int isDirectionalNotInShadow(int lightIndex, int indexOfMinT);
float getMagnitude(vec3 vec);
void printVec(vec3 vec, string str);

vec3 _u, _v, _w;
RGBQUAD color;


// New helper transformation function to transform vector by modelview 
// May be better done using newer glm functionality.
// Provided for your convenience.  Use is optional.  
// Some of you may want to use the more modern routines in readfile.cpp 
// that can also be used.  
void transformvec (const float input[4], float output[4]) 
{
  glm::vec4 inputvec(input[0], input[1], input[2], input[3]);
  glm::vec4 outputvec = modelview * inputvec;
  output[0] = outputvec[0];
  output[1] = outputvec[1];
  output[2] = outputvec[2];
  output[3] = outputvec[3];
}


/*
	Transform each objects vertices in order to place it in spacial coordinates
		- 	since we use rays to intersect object, all we need is the object's
			vertices to be transformed. Then, when we calculate the intersection
			using the ray's equation and shape's equation, the transformed
			vertices will be sufficient.
*/
void display() {

	//modelview = Transform::lookAt(eye,center,up); 

	// Transformations for objects, involving translation and scaling 
	mat4 sc(1.0) , tr(1.0), transf(1.0); 
	sc = Transform::scale(sx,sy,1.0); 
	tr = Transform::translate(tx,ty,0.0); 

	// transformation applied
	transf = tr * sc;

	//apply transformation to each vertex in object
  	for (int i = 0 ; i < numobjects ; i++) {

    	object* obj = &(objects[i]); // Grabs an object struct.
    	//apply object transform to each vertex
    	for (int i = 0; i < obj->shapeVertices.size(); ++i)
    	{
    		obj->shapeVertices[i] = obj->transform * obj->shapeVertices[i];
    	}
    }

    //construct the camera
	_w = glm::normalize(eye-center);        	// eye
	_u = glm::normalize( glm::cross(up, _w) ); 	// direction from eye to center
	_v = glm::cross(_w, _u);                  	// up direction

	vec3 camera = glm::normalize( vec3(eye) );

	//cerr << "center: " << center.x <<center.y <<center.z << " eye: " << eye.x << eye.y<< eye.z << " up: " << up.x << up.y<< up.z<< endl;

	//cerr << "u: " << _u.x <<_u.y <<_u.z << " v: " << _v.x << _v.y<< _v.z << " w: " << _w.x << _w.y<< _w.z<< endl;

    rayTrace(camera);
}


void rayTrace(vec3 camera) {

	FreeImage_Initialise();
	FIBITMAP* bitmap = FreeImage_Allocate(w, h, 24);

	//memory allocation check
	if(!bitmap){
		cout << "Can't allocate image???" << endl; 
		exit(1);
	}

	// shoot a ray through every pixel on the image
	for (float y = 0.5; y < h; ++y)
	{
		for (float x = 0.5; x < w; ++x)
		{
			float fovx, beta, alpha;
			vec3 direction, foundColor;
			Intersection hit;

			// float x = 105;
			// float y = 250;

			//generate weights
			fovx = 2.0 * (atan(tan(fovy/2.0) * (float)w/h));
			alpha = tan(fovx/2.0) * (((float)x-w/2.0)/((float)w/2.0));
			beta = tan(fovy/2.0) * (((float)y-h/2.0)/((float)h/2.0));

			//calculate ray equation in world coordinates NEW from Office Hours
			direction = vec3( alpha*_u + beta*_v - _w );
			direction = glm::normalize(direction);

			// cerr << "fovx: " << fovx << endl;
			// cerr << "alpha: " << alpha << endl;
			// cerr << "beta: " << beta << endl;
			// cerr << "fovy: " << fovy << endl;

			Ray ray(camera, direction);	

			// find out if ray intersects object geometry
			hit = intersect(ray);
			foundColor = findColor(hit, 1, camera);

			color.rgbRed = foundColor[0];
			color.rgbGreen = foundColor[1];
			color.rgbBlue = foundColor[2];

			//cerr << foundColor[0] << " " << foundColor[1] << " " << foundColor[2] << endl;

			FreeImage_SetPixelColor(bitmap, x, y, &color);			

		}
	}



		if(FreeImage_Save(FIF_PNG, bitmap, fileName, 0)){
			cout << "Image successfully saved in filename" << endl;
		}

	//housekeeping
	FreeImage_DeInitialise();
}

//determine whether the ray through pixel x, y intersects geometry
Intersection intersect(Ray ray) {

	bool hit = false;
	float minT;
	int indexOfMinT=-1;
	Intersection intersection;
	vec3 AP, BP, CP, P, cross, n, A, B, C;
	float Aw, Bw, Cw, alpha, beta, gamma, t;

	//check against each object geometry
	for (int ind = 0; ind < numobjects; ++ind)
	{
		object obj = objects[ind];

		//triangles
		if (obj.type == tri)
		{
			//define variables for easier understanding
			A = vec3(obj.shapeVertices[0].x, obj.shapeVertices[0].y, obj.shapeVertices[0].z);
			B = vec3(obj.shapeVertices[1].x, obj.shapeVertices[1].y, obj.shapeVertices[1].z );
			C = vec3(obj.shapeVertices[2].x, obj.shapeVertices[2].y, obj.shapeVertices[2].z); 

			//calculate normal and normalize it
			cross = glm::cross( (C-A), (B-A));
			n = glm::normalize( cross );

			//calculating ray plane intersection
			t = ((float)(glm::dot(A, n) - glm::dot(ray.p0, n))) / (float)(glm::dot(ray.p1, n));
			P = ray.p0 + ray.p1 * t;

			objects[ind].point = P;
			objects[ind].normal = glm::normalize(glm::cross(C-B, A-C));
			
			//calculate barycentric coordinates
			AP = (glm::cross(n, C-B)) / (glm::dot(glm::cross(n, C-B), A-C));
			Aw = glm::dot(AP, C) * -1;
			BP = (glm::cross(n, A-C)) / (glm::dot(glm::cross(n, A-C), B-A));
			Bw = glm::dot(BP, A) * -1;
			CP = (glm::cross(n, B-A)) / (glm::dot(glm::cross(n, B-A), C-B));
			Cw = glm::dot(CP, B) * -1;

			//calculate weights
			alpha = glm::dot(AP, P) + Aw;
			beta = glm::dot(BP, P) + Bw;
			gamma = glm::dot(CP, P) + Cw;

			//if intersection, weights all between 0 and 1
			if (alpha >= 0 && beta >= 0 && gamma >= 0 && alpha <= 1 && beta <= 1 && gamma <= 1 ){

				if(t>0 && (indexOfMinT==-1 || t < minT)){
					hit = true;
					minT = t;
					indexOfMinT = ind;
				}
			}
		} 
		else if (obj.type == sphere)
		{
			float a, b, c, radius, pos, neg;
			vec3 center, P, n;

			//extract center
			center = vec3(obj.shapeVertices[0].x, obj.shapeVertices[0].y, obj.shapeVertices[0].z);
			radius = obj.shapeVertices[0].w;

			//cast to appropriate types
			mat4 inverseTransf = inverse(obj.transform);
			vec3 p0Transf = vec3(inverseTransf * vec4(ray.p0,1.0));
			vec3 p1Transf = vec3(inverseTransf * vec4(ray.p1,0.0));
			vec3 centerTransf = vec3(inverseTransf * vec4(center,1.0));

			// calculating a,b,c to be used on quadratic equation
			a = glm::dot(p1Transf, p1Transf);
			b = 2.0 * glm::dot(p1Transf, (p0Transf - centerTransf) );
			c = glm::dot( (p0Transf-centerTransf), (p0Transf-centerTransf) ) - (radius*radius);

			//solving for positive and negative signed versions of quadratic equation
			pos = (-b + sqrt( b*b - 4.0*a*c )) / (2.0*a);
			neg = (-b - sqrt( b*b - 4.0*a*c )) / (2.0*a);

			//determining what to do based on roots
			//complex roots
			if ( ( b*b - 4.0*a*c ) < 0 ) {
				continue;
			}
			//2 positive roots and not equal
			else if (pos > 0 && neg > 0 && pos != neg){
				//pick smaller root
				if (pos>neg){ t = neg;	} 
				else 		{ t = pos; 	}
				hit = true;
			} 
			//roots are equal
			else if (pos == neg){
				t = pos;
				hit = true;
			}
			// one positive, one negative root
			else if ((pos > 0 && neg < 0) || (pos < 0 && neg > 0)) {
				// pick positive root
				if (pos > 0){	t = pos;	} 
				else		{	t = neg;	}
				hit = true;
			} 

			//calculate intersection point & normal on sphere (not ellipse)
			P = p0Transf + t * p1Transf;
			n = glm::normalize(P - centerTransf);

			//transform point and normal back to world coords
			objects[ind].point = vec3(obj.transform * vec4(P, 1));
			objects[ind].normal = vec3(inverse(transpose(obj.transform)) * vec4(n, 1));

			//find if t was minimum
			if(t>0 && (indexOfMinT==-1 || t < minT)){
				hit = true;
				minT = t;
				indexOfMinT = ind;
			}	
		}
	}

	if(hit){
		intersection = Intersection(objects[indexOfMinT].point, objects[indexOfMinT].normal, indexOfMinT);
		return intersection;
	}

	intersection = Intersection();
	return intersection;
}


vec3 findColor(Intersection hit, int depth, vec3 eye) {

        float red, green, blue;
        vec3 result, lambert, intersectionPoint, reflective, reflectiveDirection, reflectedColor; 
        int indexOfMinT;
        object obj;
        Intersection reflectiveHit;

        if (depth > maxDepth)
        {
        	return vec3(0,0,0);
        }

        indexOfMinT = hit.ind;

        if (indexOfMinT >= 0)
        {
        	vec3 sumOfLighting;
        	int s;

        	// check all point lights
        	for (int lightIndex = 0; lightIndex < numPointLights; ++lightIndex)
        	{
        		float distToLight, atten, notInShadow, nDotL, nDotH;
        		vec3 res, half, dirToLight, intersectionPoint, myDiffuse, phong, lambert, lightColor, mySpecular;

        		//find if intersection point in shadow
        		notInShadow = (float)isPointNotInShadow(lightIndex, indexOfMinT);

        		// calculate vector from intersection pt to light
        		dirToLight = pointLights[lightIndex] - hit.point;
        		half = glm::normalize((vec3(eye) - hit.point) + dirToLight);
        		distToLight = getMagnitude(dirToLight);
        		dirToLight = glm::normalize(dirToLight);
        		atten = (float)1.0 / (objects[hit.ind].attenuation[0] + objects[hit.ind].attenuation[1] * 
        				distToLight + objects[hit.ind].attenuation[2] * pow(distToLight, 2));
        		lightColor = pointColors[lightIndex];

        		//lambert
        		myDiffuse = (float)255.0 * vec3(objects[hit.ind].diffuse[0],objects[hit.ind].diffuse[1],objects[hit.ind].diffuse[2]);
        		nDotL = glm::dot(hit.normal, dirToLight);
        		lambert = myDiffuse * max(nDotL, (float)0.0);

        		//phong
        		mySpecular = (float)255.0 * vec3(objects[hit.ind].specular[0],objects[hit.ind].specular[1],objects[hit.ind].specular[2]);;
        		nDotH = glm::dot(hit.normal, half);
        		shininess = objects[hit.ind].shininess;
        		phong = mySpecular * pow( max(nDotH, (float)0.0), shininess );

        		//summation part 
        		res = vec3( (lambert + phong) * atten );
        		
        		//increment lighting sum so far
        		sumOfLighting = sumOfLighting + notInShadow*res;
        	}

        	// check all directional lights
        	for (int lightIndex = 0; lightIndex < numDirectionalLights; ++lightIndex)
        	{
        		float distToLight, atten, notInShadow, nDotL, nDotH;
        		vec3 res, half, dirToLight, intersectionPoint, myDiffuse, phong, lambert, lightColor, mySpecular;

        		//find if intersection point in shadow
        		notInShadow = (float)isDirectionalNotInShadow(lightIndex, indexOfMinT);

        		// calculate vector from intersection pt to light
        		dirToLight = directionalLights[lightIndex];
        		half = glm::normalize((vec3(eye) - hit.point) + dirToLight);
        		distToLight = getMagnitude(dirToLight);
        		dirToLight = glm::normalize(dirToLight);
        		atten = (float)1.0 / (objects[hit.ind].attenuation[0] + objects[hit.ind].attenuation[1] * 
        				distToLight + objects[hit.ind].attenuation[2] * pow(distToLight, 2));
        		lightColor = directionalColors[lightIndex];

        		//lambert
        		myDiffuse = (float)255.0 * vec3(objects[hit.ind].diffuse[0],objects[hit.ind].diffuse[1],objects[hit.ind].diffuse[2]);
        		nDotL = glm::dot(hit.normal, dirToLight);
        		lambert = myDiffuse * max(nDotL, (float)0.0);

        		//phong
        		mySpecular = (float)255.0 * vec3(objects[hit.ind].specular[0],objects[hit.ind].specular[1],objects[hit.ind].specular[2]);;
        		nDotH = glm::dot(hit.normal, half);
        		shininess = objects[hit.ind].shininess;
        		phong = mySpecular * pow( max(nDotH, (float)0.0), shininess );

        		// cerr << "atten: " << atten << endl;
        		// printVec(myDiffuse, "myDiffuse");
        		// printVec(mySpecular, "mySpecular");


        		//summation part 
        		res = vec3( (lambert + phong) * atten );
        		
        		//increment lighting sum so far
        		sumOfLighting = sumOfLighting + notInShadow*res;
        	}

        	// calculate color from intersection
        	red = objects[indexOfMinT].ambient[0]*255.0 + objects[indexOfMinT].emission[0]*255.0 + sumOfLighting.x;
        	red = red > 255 ? 255: red;
        	green = objects[indexOfMinT].ambient[1]*255.0 + objects[indexOfMinT].emission[1]*255.0 + sumOfLighting.y;
        	green = green > 255 ? 255: green;
        	blue = objects[indexOfMinT].ambient[2]*255.0 + objects[indexOfMinT].emission[2]*255.0 + sumOfLighting.z;
        	blue = blue > 255 ? 255: blue;

        	//save clipped colors
			result[0] = red;
			result[1] = green;
			result[2] = blue;

			// //calculate new recursive hit point from pt of intersection's mirror direction
			// intersectionPoint = objects[indexOfMinT].point + objects[indexOfMinT].normal*(float)0.00001;
			// //calculate ray direction
			// reflectiveDirection = ;

			// //shoot new ray to calculate reflective lighting
			// Ray reflectiveRay(intersectionPoint, reflectiveDirection);
			// reflectiveHit = intersect(reflectiveRay);
			// reflectedColor = findColor(reflectiveHit, depth+1, intersectionPoint);

			// // calculate sum of reflected colors on this object
			// reflective = vec3( 	objects[indexOfMinT].specular[0] * reflectedColor[0] * 255.0, 
			// 					objects[indexOfMinT].specular[1] * reflectedColor[1] * 255.0, 
			// 					objects[indexOfMinT].specular[2] * reflectedColor[2] * 255.0  );

        }
        else {
			result[0] = 0.0;
			result[1] = 0.0;
			result[2] = 0.0;

			reflective = vec3(0,0,0);
		}


		return result /*+ reflective*/;

}

int isPointNotInShadow(int lightIndex, int indexOfMinT) {

	vec3 intersectionPoint, direction, normalizedDirection, objectVec;
	Intersection objectHit = Intersection();

	//create ray
	intersectionPoint = objects[indexOfMinT].point + objects[indexOfMinT].normal*(float)0.00001;
	direction = pointLights[lightIndex] - intersectionPoint;
	normalizedDirection = glm::normalize(direction);
	Ray ray(intersectionPoint, normalizedDirection);

	//intersect ray with objects 
	objectHit = intersect(ray);

	if (objectHit.ind == -1)
	{
		//not in shadow
		// cerr << "not in shadow\n";
		return 1;
	}

	//compare manitude of object intersection with point light
	objectVec = objectHit.point - intersectionPoint;

	//check if the length to the object is smaller than length to light
	if (getMagnitude(direction) > getMagnitude(objectVec))
	{
		// in shadow
		// cerr << "in shadow\n";
		return 0;
	}

	//not in shadow
	// cerr << "not in shadow\n";
	return 1;

}


int isDirectionalNotInShadow(int lightIndex, int indexOfMinT) {

	vec3 intersectionPoint, direction, normalizedDirection, objectVec;
	Intersection objectHit = Intersection();

	//create ray
	intersectionPoint = objects[indexOfMinT].point;
	direction = directionalLights[lightIndex];
	normalizedDirection = glm::normalize(direction);
	Ray ray(intersectionPoint, normalizedDirection);

	//intersect ray with objects 
	objectHit = intersect(ray);

	//didnt hit any objects
	if (objectHit.ind == -1)
	{
		//not in shadow
		return 1;
	}

	//in shadow
	return 0;

}

float getMagnitude(vec3 vec) {
	return sqrt(vec.x*vec.x + vec.y*vec.y + vec.z*vec.z);
}

void printVec(vec3 vec, string str){
	cerr << str << ": " << vec.x << " " << vec.y << " " << vec.z << endl;
}

















