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


bool intersect(Ray);
void rayTrace(vec3);
bool drawSquare(int, int);

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
	_u = glm::normalize(glm::cross(up, _w)); 	// direction from eye to center
	_v = glm::cross(_w, _u);                  	// up direction

	vec3 camera = vec3 (eye);

	//cerr << "center: " << center.x <<center.y <<center.z << " eye: " << eye.x << eye.y<< eye.z << " up: " << up.x << up.y<< up.z<< endl;

	//cerr << "u: " << _u.x <<_u.y <<_u.z << " v: " << _v.x << _v.y<< _v.z << " w: " << _w.x << _w.y<< _w.z<< endl;

    rayTrace(camera);
}

/*
	questions:
		- how do we draw the image?
			freeImage
		- what is the camera (eye?)
			up, eye, and center vectors to construct cam... still not sure
*/

void rayTrace(vec3 camera) {

	FreeImage_Initialise();
	FIBITMAP* bitmap = FreeImage_Allocate(w, h, 24);

	//memory allocation check
	if(!bitmap){
		cout << "Can't allocate image???" << endl; 
		exit(1);
	}

	// shoot a ray through every pixel on the image
	for (int y = 0; y < h; ++y)
	{
		for (int x = 0; x < w; ++x)
		{
			bool hit;

			//generate weights
			float fovx = 2 * (atan(tan(fovy/2) * (float)w/h));
			float alpha = tan(fovx/2) * (((float)x-w/2)/((float)w/2));
			float beta = tan(fovy/2) * (((float)y-h/2)/((float)h/2));

			//calculate ray equation in world coordinates NEW from Office Hours
			vec3 direction = vec3( alpha*_u + beta*_v - _w );
			direction = glm::normalize(direction);
			Ray ray(camera, direction);	

			// find out if ray intersects object geometry
			if(intersect(ray))
			{
				FreeImage_SetPixelColor(bitmap, x, y, &color);			
			} 
			else {
				FreeImage_SetPixelColor(bitmap, x, y, &color);			
			}

		}
	}



	//draw image
	if(FreeImage_Save(FIF_PNG, bitmap, "test11.png", 0)){
		cout << "Image successfully saved" << endl;
	}
	//housekeeping
	FreeImage_DeInitialise();
}

//determine whether the ray through pixel x, y intersects geometry
bool intersect(Ray ray) {

	bool hit = false;
	float minT;
	int indexOfMinT=-1;

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

			//cerr << alpha << " " << beta << " " << gamma << endl;

			//if intersection, weights all between 0 and 1
			if (alpha >= 0 && beta >= 0 && gamma >= 0 && alpha <= 1 && beta <= 1 && gamma <= 1 ){

				if(t>0 && (indexOfMinT==-1 || t < minT)){
					hit = true;
					minT = t;
					indexOfMinT = ind;
				}
				//cerr << "Intersected triangle " << ind+1 << ". HIT!!!" << endl;

			}
			//cerr << "didn't intersect triangle " << ind+1 << endl;
			
		} 
		else if (obj.type == sphere)
		{

			float a, b, c, radius, pos, neg;
			vec3 center;

			center = vec3(obj.shapeVertices[0].x, obj.shapeVertices[0].y, obj.shapeVertices[0].z);

			mat4 inverseTransf = inverse(obj.transform);
			vec3 p0Transf = vec3(inverseTransf * vec4(ray.p0,1));
			vec3 p1Transf = vec3(inverseTransf * vec4(ray.p1,0));
			vec3 centerTransf = vec3(inverseTransf * vec4(center,1));

			// extract center
			radius = obj.shapeVertices[0].w;

			// calculating a,b,c to be used on quadratic equation
			a = glm::dot(p1Transf, p1Transf);
			b = 2.0 * glm::dot(p1Transf, (p0Transf - centerTransf) );
			c = glm::dot( (p0Transf-centerTransf), (p0Transf-centerTransf) ) - (radius*radius);

			//solving for positive and negative signed versions of quadratic equation
			pos = (-b + sqrt( b*b - 4.0*a*c )) / (2*a);
			neg = (-b - sqrt( b*b - 4.0*a*c )) / (2*a);

			//determining what to do based on roots
			//complex roots
			if ( ( b*b - 4.0*a*c ) < 0 )
			{
				continue;
			}
			//2 positive roots and not equal
			else if (pos > 0 && neg > 0 && pos != neg)
			{
				//pick smaller root
				if (pos>neg)
				{
					t = neg;	
				} else {
					t = pos;
				}
				hit = true;
			} 
			//roots are equal
			else if (pos == neg)
			{
				t = pos;
				hit = true;
			}
			// one positive, one negative root
			else if ((pos > 0 && neg < 0) || (pos < 0 && neg > 0))
			{
				// pick positive root
				if (pos > 0)
				{
					t = pos;
				} else	{
					t = neg;
				}
				hit = true;
			} 

			//find if t was minimum
			if(t>0 && (indexOfMinT==-1 || t < minT)){
				hit = true;
				minT = t;
				indexOfMinT = ind;
			}

			
		} else {
			cerr << "Incorrect way to tell which type of object it is while intersecting in display.cpp\n";
		}
	}

	if(hit){
		color.rgbRed = objects[indexOfMinT].ambient[0]*255.0;
		color.rgbGreen = objects[indexOfMinT].ambient[1]*255.0;
		color.rgbBlue = objects[indexOfMinT].ambient[2]*255.0;
		return 1;
	}

	color.rgbRed = 0.0;
	color.rgbGreen = 0.0;
	color.rgbBlue = 0.0;
	return 0;
}



















