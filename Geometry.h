//#include "Transform.h"
#include <glm/glm.hpp>

// glm provides vector, matrix classes like glsl
// Typedefs to make code more readable 
typedef glm::vec3 vec3 ; 


class Ray {
  public:
    vec3 p0, p1;
    Ray (vec3,vec3);
};




class Sphere {
	vec3 center;
	float radius;
  public:
  	Sphere(vec3 c, float r);
};



class Triangle {
	vec3 p0, p1, p2;
  public:
  	Triangle(vec3, vec3, vec3);
};


class Intersection {
public:
  vec3 point, normal;
  int ind;
  Intersection(vec3, vec3, int);
  Intersection();
};