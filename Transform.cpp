// Transform.cpp: implementation of the Transform class.

// Note: when you construct a matrix using mat4() or mat3(), it will be COLUMN-MAJOR
// Keep this in mind in readfile.cpp and display.cpp
// See FAQ for more details or if you're having problems.

#include "Transform.h"
#include <stdio.h>
#include <iostream>

using namespace std;



// Helper rotation function.  Please implement this.  
mat3 Transform::rotate(const float degrees, const vec3& axis) 
{
  //copied from HW1
    float pi = 3.1415926535897932384626433832795;
    float radians = degrees * pi / 180.0; //degrees in radians    


    mat3 identity = mat3(); // 3x3 identity matrix
    mat3 a_a_transpose = mat3(
      axis[0]*axis[0],  axis[0]*axis[1],  axis[0]*axis[2],
      axis[0]*axis[1],  axis[1]*axis[1],  axis[1]*axis[2],
      axis[0]*axis[2],  axis[1]*axis[2],  axis[2]*axis[2]
      ); //as per in class formula

    mat3 a_star = mat3(
      0,  -axis[2], axis[1],
      axis[2],  0,  -axis[0],
      -axis[1], axis[0],  0
      );  // as per in class formula
    
    // using the formula seen in class: 
    mat3 result = glm::cos(radians) * identity + (1-glm::cos(radians)) * a_a_transpose + glm::sin(radians) * a_star;

  return result;
}

void Transform::left(float degrees, vec3& eye, vec3& up) 
{
  //copied from HW1
  mat3 rMat = rotate(degrees, up); //rotation matrix
  eye = eye * rMat;   // update eye
  up = up * rMat; 
  
}

void Transform::up(float degrees, vec3& eye, vec3& up) 
{
  //copied from hw1
  vec3 normalized = glm::normalize(glm::cross(eye,up));
  mat3 rMat = rotate(-degrees, normalized); // rotation matrix

  eye = rMat * eye; // update eye
  
  up = glm::normalize(rMat*up); //update up
}

mat4 Transform::lookAt(const vec3 &eye, const vec3 &center, const vec3 &up) 
{
  //copied from hw1
  vec3 w = glm::normalize(eye);               // eye
  vec3 u = glm::normalize(glm::cross(up, w)); // direction from eye to center
  vec3 v = glm::cross(w, u);                  // up direction

  //transform into 4x4 with extra column for transformation and extra row for w.
  mat4 rot = mat4(u.x, v.x, w.x, 0, u.y, v.y, w.y, 0, u.z, v.z, w.z, 0, 0, 0, 0, 1);

  // fill transformation matrix as per in slides.
  mat4 trans = mat4(1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, -eye.x, -eye.y, -eye.z, 1);

  return rot * trans;
}

mat4 Transform::perspective(float fovy, float aspect, float zNear, float zFar)
{

  mat4 ret;
  float pi = 3.1415926535897932384626433832795;

  double convDegToRad = (pi / 180) * (fovy / 2); 
  double tanOfRad = glm::tan(convDegToRad);
  double d = 1/tanOfRad;
  double din = d / aspect;
  //from lecture slides
  double a = ((zFar + zNear) / (-1 * zNear+ zFar))*-1;
  double b = ((zFar * zNear * 2) / (-1 * zNear+ zFar))*-1;

  ret =  mat4(din, 0, 0, 0,
                      0, d, 0, 0,
                      0, 0, a, -1,
                      0, 0, b, 0);

  return ret;
}

mat4 Transform::scale(const float &sx, const float &sy, const float &sz) 
{
  //from lecture slides
  mat4 ret = mat4( sx, 0, 0, 0,
                    0, sy, 0, 0,
                    0, 0, sz, 0,
                    0, 0, 0, 1);
  return ret;
}

mat4 Transform::translate(const float &tx, const float &ty, const float &tz) 
{
  //from lecture slides
  mat4 ret;
  ret = mat4( 1, 0, 0, tx,
              0, 1, 0, ty,
              0, 0, 1, tz,
              0, 0, 0, 1);
  ret = glm::transpose(ret);
  return ret;
}

// To normalize the up direction and construct a coordinate frame.  
// As discussed in the lecture.  May be relevant to create a properly 
// orthogonal and normalized up. 
// This function is provided as a helper, in case you want to use it. 
// Using this function (in readfile.cpp or display.cpp) is optional.  

vec3 Transform::upvector(const vec3 &up, const vec3 & zvec) 
{
  cerr << "up: " << up.x <<up.y <<up.z << endl;
  cerr << "zvec: " << zvec.x <<zvec.y <<zvec.z << endl;
  vec3 x = glm::cross(up,zvec); 
  cerr << "x: " << x.x <<x.y <<x.z << endl;

  vec3 y = glm::cross(zvec,x); 
    cerr << "y: " << y.x <<y.y <<y.z << endl;

  vec3 ret = glm::normalize(y); 
  cerr << "ret: " << ret.x <<ret.y <<ret.z << endl;
  return ret; 
}


Transform::Transform()
{

}

Transform::~Transform()
{

}
