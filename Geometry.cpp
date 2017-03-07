#include "Geometry.h"



Ray::Ray(vec3 start_point, vec3 direction) {
	p0 = start_point;
	p1 = direction;
}




Sphere::Sphere(vec3 c, float r) {
	center = c;
	radius = r;
}




Triangle::Triangle(vec3 p_0, vec3 p_1, vec3 p_2)  {
	p0 = p_0;
	p1 = p_1;
	p2 = p_2;
}