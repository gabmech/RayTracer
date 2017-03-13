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
	for (int x = 0; x < w; ++x)
	{
		for (int y = 0; y < h; ++y)
		{
			bool hit;
			RGBQUAD color;

			//generate weights
			float fovx = 2 * (atan(tan(fovy/2) * (float)w/h));
			float alpha = tan(fovx/2) * (((float)x-w/2)/((float)w/2));
			float beta = tan(fovy/2) * ((h/2-y)/(h/2));
			
			 	cerr << _w.x << _w.y << _w.z << endl;
			 	cerr << center.x << center.y << center.z << endl;
			 	cerr << eye.x << eye.y << eye.z << endl;

			//calculate ray equation in world coordinates NEW from Office Hours
			vec3 direction = vec3( alpha*_u + beta*_v - _w );
			direction = glm::normalize(direction);
			Ray ray(camera, direction);	

			// if (x == w/2 && y == h/2)
			// {
			// 	cerr << _w.x << _w.y << _w.z << endl;
			// 	cerr << "Direction: " << direction.x << " " << direction.y << " " << direction.z << endl;
			// 	int a1 = intersect(ray);
			// 	cerr << "intersection: " << a1 << endl;

			// } else {
			// 	continue;
			// }

			//find out if ray intersects object geometry
			if(intersect(ray))
			{
				color.rgbRed = 0.0;
				color.rgbGreen = (double) 255.0 * x/w;
				color.rgbBlue = (double)  255.0 * x/w;
				FreeImage_SetPixelColor(bitmap, x, y, &color);			
			} 
			else {
				color.rgbRed = 0.0;
				color.rgbGreen = 0.0;
				color.rgbBlue = 0.0;
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

	vec3 AP, BP, CP, P, cross, n, A, B, C;
	float Aw, Bw, Cw, alpha, beta, gamma, t;

	cerr << "NUMO: " << numobjects << endl;

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

			cerr << "triangle: " << ind+1 << endl;

			//calculate normal and normalize it
			cross = glm::cross( (C-A), (B-A));
			n = glm::normalize( cross );

			//calculating ray plane intersection
			t = (glm::dot(A, n) - glm::dot(ray.p0, n)) / (glm::dot(ray.p1, n));
			P = ray.p0 + ray.p1 * t;

			//calculate barycentric coordinates
			AP = glm::normalize((glm::cross(n, C-B)) / (glm::dot(glm::cross(n, C-B), A-C)));
			Aw = glm::dot(AP, C) * -1;
			BP = glm::normalize((glm::cross(n, A-C)) / (glm::dot(glm::cross(n, A-C), B-A)));
			Bw = glm::dot(BP, A) * -1;
			CP = glm::normalize((glm::cross(n, B-A)) / (glm::dot(glm::cross(n, B-A), C-B)));
			Cw = glm::dot(CP, B) * -1;

			//calculate weights
			alpha = glm::dot(AP, P) + Aw;
			beta = glm::dot(BP, P) + Bw;
			gamma = glm::dot(CP, P) + Cw;

			cerr << alpha << " " << beta << " " << gamma << endl;

			//if intersection, weights all between 0 and 1
			if (alpha >= 0 && beta >= 0 && gamma >= 0 && alpha <= 1 && beta <= 1 && gamma <= 1 ){
				cerr << "Intersected triangle " << ind+1 << ". HIT!!!" << endl;
				return 1;
			}
			cerr << "didn't intersect triangle " << ind+1 << endl;
			
		} 
		else if (obj.type == sphere)
		{
			
		} else {
			cerr << "Incorrect way to tell which type of object it is while intersecting in display.cpp\n";
		}
	}
	return 0;
}



















