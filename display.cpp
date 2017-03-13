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



bool intersect(Ray, int, int);
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

    // construct camera
	_w = glm::normalize(eye-center);        // eye
	_u = glm::normalize(glm::cross(up, _w)); // direction from eye to center
	_v = glm::cross(_w, _u);                  // up direction

	vec3 camera(0, 0, -1);

    rayTrace(camera	);
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

	// camera, scene, width, height
	for (int x = 0; x < w; ++x)
	{
		for (int y = 0; y < h; ++y)
		{
			bool hit;
			RGBQUAD color;

			float fovx = 2 * (atan(tan(fovy/2) * (float)w/h));
			float alpha = tan(fovx/2) * (((float)x-w/2)/((float)w/2));
			float beta = tan(fovy/2) * ((h/2-y)/(h/2));
			float gamma = 1-beta-alpha;

			vec3 direction = vec3( alpha*_u + beta*_v -_w );
			direction = glm::normalize(direction);

			Ray ray(camera, direction);	

			// if (drawSquare(i, j))
			if(intersect(ray, x, y))
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

	if(FreeImage_Save(FIF_PNG, bitmap, "test11.png", 0)){
		cout << "Image successfully saved" << endl;
	}
	FreeImage_DeInitialise();
}


bool intersect(Ray ray, int x, int y) {
	for (int ind = 0; ind < numobjects; ++ind)
	{
		object obj = objects[ind];
		if (obj.type == tri)
		{
			//define variables for easier understanding
			vec3 A = vec3(obj.shapeVertices[1].x, obj.shapeVertices[1].y, obj.shapeVertices[1].z),
				 C = vec3(obj.shapeVertices[0].x, obj.shapeVertices[0].y, obj.shapeVertices[0].z), 
				 B = vec3(obj.shapeVertices[2].x, obj.shapeVertices[2].y, obj.shapeVertices[2].z );

			//calculate normal and normalize it
			vec3 cross = glm::cross( (C-A), (B-A));
			vec3 n = glm::normalize( cross );

			//use normal in plane equation
			float t = (glm::dot(A, n) - glm::dot(ray.p0, n)) / (glm::dot(ray.p1, n));
			vec3 P = ray.p0 + ray.p1 * t;


			//find weights			
			float fovx = 2 * (atan(tan(fovy/2) * (float)w/h));
			float alpha = tan(fovx/2) * (((float)x-w/2)/((float)w/2));
			float beta = tan(fovy/2) * ((h/2-y)/(h/2));
			float gamma = 1-beta-alpha;
			//plugging into formula given
			vec3 calcedP = alpha * A + beta * B + gamma * C;

			cerr << "alpha: " << alpha << "\tA: " << A.x << ", " << A.y << ", " << A.z << endl;
			cerr << "beta: " << beta << "\tB: " << B.x << ", " << B.y << ", " << B.z << endl;
			cerr << "gamma: " << gamma << "\tC: " << C.x << ", " << C.y << ", " << C.z << endl;
			//cerr << "calcedP: " << calcedP.x << ", " << calcedP.y << ", " << calcedP.z << endl;
			//cerr << "P: " << P.x << ", " << P.y << ", " << P.z << endl;

			cerr << "Sum weights: " << alpha + beta + gamma << endl; 

			//if intersection, weights all greater than 0
			if (alpha <= 0 || beta <= 0 || gamma <= 0 ){
				cerr << "MISS! ray intersected plane outside of triangle\n";
				return 0;
			}

			cerr << "HIT! ray intersects plane inside triange" << x << endl;
			return 1;
			

		} 
		else if (obj.type == sphere)
		{
			
		} else {
			cerr << "Incorrect way to tell which type of object it is while intersecting in display.cpp\n";
		}
	}
	return 1;
}



















