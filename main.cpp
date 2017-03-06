/******************************************************************************/
/* This is the program skeleton for homework 2 in CSE167 by Ravi Ramamoorthi  */
/* Extends HW 1 to deal with shading, more transforms and multiple objects    */
/******************************************************************************/

// You shouldn't have to edit this file at all!

#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <deque>
#include <stack>
#ifdef __APPLE__
#include <OpenGL/gl3.h>
#include <OpenGL/glext.h>
#include <GLUT/glut.h>
#else
#include <GL/glew.h>
#include <GL/glut.h>
#endif
#include "shaders.h"
#include "Transform.h"
#include <FreeImage.h>
#include "UCSD/grader.h"
#include "Geometry.h"

using namespace std; 

// Main variables in the program.  
#define MAINPROGRAM 
#include "variables.h" 
#include "readfile.h" // prototypes for readfile.cpp  
void display(void);  // prototype for display function.  

Grader grader;
bool allowGrader = false;

// Reshapes the window
void reshape(int width, int height){
  w = width;
  h = height;

  glViewport(0, 0, w, h);

  float aspect = (float) w / (float) h, zNear = 0.1, zFar = 99.0 ;
  // I am changing the projection matrix to fit with the new window aspect ratio
  if (useGlu) projection = glm::perspective(glm::radians(fovy),aspect,zNear,zFar) ; 
  else {
	  projection = Transform::perspective(fovy,aspect,zNear,zFar) ;
  }
  // Now send the updated projection matrix to the shader
  glUniformMatrix4fv(projectionPos, 1, GL_FALSE, &projection[0][0]);
}

void saveScreenshot(string fname) {
  int pix = w * h;
  BYTE *pixels = new BYTE[3*pix];	
  glReadBuffer(GL_FRONT);
  glReadPixels(0,0,w,h,GL_BGR,GL_UNSIGNED_BYTE, pixels);

  FIBITMAP *img = FreeImage_ConvertFromRawBits(pixels, w, h, w * 3, 24, 0xFF0000, 0x00FF00, 0x0000FF, false);

  std::cout << "Saving screenshot: " << fname << "\n";

  FreeImage_Save(FIF_PNG, img, fname.c_str(), 0);
  delete[] pixels;
}


void printHelp() {
  std::cout << "\npress 'h' to print this message again.\n" 
    << "press '+' or '-' to change the amount of rotation that\noccurs with each arrow press.\n" 
    << "press 'i' to run image grader test cases\n"
    << "press 'g' to switch between using glm::lookAt and glm::Perspective or your own LookAt.\n"       
    << "press 'r' to reset the transformations.\n"
    << "press 'v' 't' 's' to do view [default], translate, scale.\n"
    << "press ESC to quit.\n" ;      
}


void init() {
  // Initialize shaders
  vertexshader = initshaders(GL_VERTEX_SHADER, "shaders/light.vert.glsl") ;
  fragmentshader = initshaders(GL_FRAGMENT_SHADER, "shaders/light.frag.glsl") ;
  shaderprogram = initprogram(vertexshader, fragmentshader) ; 
  // Get locations of all uniform variables.
  enablelighting = glGetUniformLocation(shaderprogram,"enablelighting") ;
  lightpos = glGetUniformLocation(shaderprogram,"lightposn") ;       
  lightcol = glGetUniformLocation(shaderprogram,"lightcolor") ;       
  numusedcol = glGetUniformLocation(shaderprogram,"numused") ;       
  ambientcol = glGetUniformLocation(shaderprogram,"ambient") ;       
  diffusecol = glGetUniformLocation(shaderprogram,"diffuse") ;       
  specularcol = glGetUniformLocation(shaderprogram,"specular") ;       
  emissioncol = glGetUniformLocation(shaderprogram,"emission") ;       
  shininesscol = glGetUniformLocation(shaderprogram,"shininess") ;    
  projectionPos = glGetUniformLocation(shaderprogram, "projection");
  modelviewPos = glGetUniformLocation(shaderprogram, "modelview");
  // Initialize geometric shapes
  initBufferObjects();
  initTeapot(); initCube(); initSphere();
}

int main(int argc, char* argv[]) {

  if (argc < 2) {
    cerr << "Usage: transforms scenefile [grader input (optional)]\n"; 
    exit(-1); 
  }

  FreeImage_Initialise();
  glutInit(&argc, argv);
// OSX systems require an extra window init flag
#ifdef __APPLE__
  glutInitDisplayMode(GLUT_3_2_CORE_PROFILE | GLUT_DOUBLE | GLUT_RGBA | GLUT_DEPTH);
#else
  glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA | GLUT_DEPTH);
#endif
  glutCreateWindow("HW2: Scene Viewer");

#ifndef __APPLE__ // GLew not needed on OSX systems
  GLenum err = glewInit() ; 
  if (GLEW_OK != err) { 
    std::cerr << "Error: " << glewGetString(err) << std::endl; 
  } 
#endif

  init();
  readfile(argv[1]) ; 
  glutDisplayFunc(display);
  glutSpecialFunc(specialKey);
  glutKeyboardFunc(keyboard);
  glutReshapeFunc(reshape);
  glutReshapeWindow(w, h);

  if (argc > 2) {
    allowGrader = true;
    stringstream tcid;
    tcid << argv[1] << "." << argv[2];
    grader.init(tcid.str());
    grader.loadCommands(argv[2]);
    grader.bindDisplayFunc(display);
    grader.bindSpecialFunc(specialKey);
    grader.bindKeyboardFunc(keyboard);
    grader.bindScreenshotFunc(saveScreenshot);
  }

  printHelp();
  glutMainLoop();
  FreeImage_DeInitialise();
  destroyBufferObjects();
  return 0;
}
