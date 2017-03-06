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
#include "variables.h"
#include "readfile.h"
#include "display.h"


class Ray {
    vec3 p0, p1;
  public:
    Ray (vec3,vec3);
};

Ray::Ray(vec3 start_point, vec3 direction) {
	p0 = start_point;
	p1 = direction;
}

bool intersect(Ray);



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
    		obj->shapeVertices[i] = obj->shapeVertices[i] * obj->transform;
    	}
    }
}

/*
	questions:
		- how do we draw the image?
			freeImage
		- what is the camera (eye?)
			up, eye, and center vectors to construct cam... still not sure
*/

void rayTrace() {
	//Image img(w, h);
	// camera, scene, width, height
	for (int i = 0; i < h; ++i)
	{
		for (int j = 0; j < w; ++j)
		{
			bool hit;
			vec3 pixel(i, j, 0);
			vec3 camera(0, 0, -1);
			Ray ray(camera, pixel);
			if (intersect(ray))
			{
				//img[i][j] = GREEN;
			} 
			else {
				//img[i][j] = BLACK;
			}

		}
	}
}

bool intersect(Ray ray) {
	for (int i = 0; i < numobjects; ++i)
	{
		object obj = objects[i];
		if (obj.type == tri)
		{
			
		} 
		else if (obj.type == sphere)
		{
			
		} else {
			cerr << "Incorrect way to tell which type of object it is while intersecting in display.cpp\n";
		}
	}
}



















