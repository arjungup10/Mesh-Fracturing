/* Base code for texture mapping lab */
/* includes three images and three meshes - Z. Wood 2016 */
#include <iostream>
#define GLEW_STATIC
#define GLM_ENABLE_EXPERIMENTAL
#include <GL/glew.h>
#include <GLFW/glfw3.h>

#include "GLSL.h"
#include "Program.h"
#include "MatrixStack.h"
#include "Shape.h"
#include "Texture.h"


// value_ptr for glm
#include <glm/gtc/type_ptr.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtx/quaternion.hpp>

///* to use glee */
//#define GLEE_OVERWRITE_GL_FUNCTIONS
//#include "glee.hpp"

using namespace std;
using namespace glm;

GLFWwindow *window; // Main application window
string RESOURCE_DIR = ""; // Where the resources are loaded from
shared_ptr<Program> prog0;
shared_ptr<Program> prog1;
shared_ptr<Program> prog2;
shared_ptr<Shape> world;
shared_ptr<Shape> shape;

shared_ptr<Texture> texture0;
shared_ptr<Texture> texture1;
shared_ptr<Texture> texture2;

vec3 eye = vec3(0, 0, 0);
vec3 lookAtPoint = vec3(0, 0, -5);
static const vec3 upV = vec3(0, 1, 0);

float alpha = 0;
float beta = -90;
bool firstMove = true;
float lastX = 512, lastY = 512;

int g_GiboLen;
int g_width, g_height;
float cTheta = 0;
float cHeight = 0;

float myTime = 0;
bool timeForward = true;
bool timePause = 0;
float timeSpeed = .001;

vector<vec3> positions(100);
vector<float> explosions(100);
vector<quat> orientations(100);
string message;


//global data for ground plane
GLuint GrndBuffObj, GrndNorBuffObj, GrndTexBuffObj, GIndxBuffObj;

void ScreenPosToWorldRay(
	int mouseX, int mouseY,             // Mouse position, in pixels, from bottom-left corner of the window
	int screenWidth, int screenHeight,  // Window size, in pixels
	glm::mat4 ViewMatrix,               // Camera position and orientation
	glm::mat4 ProjectionMatrix,         // Camera parameters (ratio, field of view, near and far planes)
	glm::vec3& out_origin,              // Ouput : Origin of the ray. /!\ Starts at the near plane, so if you want the ray to start at the camera's position instead, ignore this.
	glm::vec3& out_direction            // Ouput : Direction, in world space, of the ray that goes "through" the mouse.
) {

	// The ray Start and End positions, in Normalized Device Coordinates (Have you read Tutorial 4 ?)
	glm::vec4 lRayStart_NDC(
		((float)mouseX / (float)screenWidth - 0.5f) * 2.0f, // [0,1024] -> [-1,1]
		((float)mouseY / (float)screenHeight - 0.5f) * 2.0f, // [0, 768] -> [-1,1]
		-1.0, // The near plane maps to Z=-1 in Normalized Device Coordinates
		1.0f
	);
	glm::vec4 lRayEnd_NDC(
		((float)mouseX / (float)screenWidth - 0.5f) * 2.0f,
		((float)mouseY / (float)screenHeight - 0.5f) * 2.0f,
		0.0,
		1.0f
	);


	// The Projection matrix goes from Camera Space to NDC.
	// So inverse(ProjectionMatrix) goes from NDC to Camera Space.
	glm::mat4 InverseProjectionMatrix = glm::inverse(ProjectionMatrix);

	// The View Matrix goes from World Space to Camera Space.
	// So inverse(ViewMatrix) goes from Camera Space to World Space.
	glm::mat4 InverseViewMatrix = glm::inverse(ViewMatrix);

	glm::vec4 lRayStart_camera = InverseProjectionMatrix * lRayStart_NDC;    lRayStart_camera /= lRayStart_camera.w;
	glm::vec4 lRayStart_world = InverseViewMatrix       * lRayStart_camera; lRayStart_world /= lRayStart_world.w;
	glm::vec4 lRayEnd_camera = InverseProjectionMatrix * lRayEnd_NDC;      lRayEnd_camera /= lRayEnd_camera.w;
	glm::vec4 lRayEnd_world = InverseViewMatrix       * lRayEnd_camera;   lRayEnd_world /= lRayEnd_world.w;


	// Faster way (just one inverse)
	//glm::mat4 M = glm::inverse(ProjectionMatrix * ViewMatrix);
	//glm::vec4 lRayStart_world = M * lRayStart_NDC; lRayStart_world/=lRayStart_world.w;
	//glm::vec4 lRayEnd_world   = M * lRayEnd_NDC  ; lRayEnd_world  /=lRayEnd_world.w;


	glm::vec3 lRayDir_world(lRayEnd_world - lRayStart_world);
	lRayDir_world = glm::normalize(lRayDir_world);


	out_origin = glm::vec3(lRayStart_world);
	out_direction = glm::normalize(lRayDir_world);
}


bool TestRayOBBIntersection(
	glm::vec3 ray_origin,        // Ray origin, in world space
	glm::vec3 ray_direction,     // Ray direction (NOT target position!), in world space. Must be normalize()'d.
	glm::vec3 aabb_min,          // Minimum X,Y,Z coords of the mesh when not transformed at all.
	glm::vec3 aabb_max,          // Maximum X,Y,Z coords. Often aabb_min*-1 if your mesh is centered, but it's not always the case.
	glm::mat4 ModelMatrix,       // Transformation applied to the mesh (which will thus be also applied to its bounding box)
	float& intersection_distance // Output : distance between ray_origin and the intersection with the OBB
) {

	// Intersection method from Real-Time Rendering and Essential Mathematics for Games

	float tMin = 0.0f;
	float tMax = 100000.0f;

	glm::vec3 OBBposition_worldspace(ModelMatrix[3].x, ModelMatrix[3].y, ModelMatrix[3].z);

	glm::vec3 delta = OBBposition_worldspace - ray_origin;

	// Test intersection with the 2 planes perpendicular to the OBB's X axis
	{
		glm::vec3 xaxis(ModelMatrix[0].x, ModelMatrix[0].y, ModelMatrix[0].z);
		float e = glm::dot(xaxis, delta);
		float f = glm::dot(ray_direction, xaxis);

		if (fabs(f) > 0.001f) { // Standard case

			float t1 = (e + aabb_min.x) / f; // Intersection with the "left" plane
			float t2 = (e + aabb_max.x) / f; // Intersection with the "right" plane
											 // t1 and t2 now contain distances betwen ray origin and ray-plane intersections

											 // We want t1 to represent the nearest intersection, 
											 // so if it's not the case, invert t1 and t2
			if (t1 > t2) {
				float w = t1; t1 = t2; t2 = w; // swap t1 and t2
			}

			// tMax is the nearest "far" intersection (amongst the X,Y and Z planes pairs)
			if (t2 < tMax)
				tMax = t2;
			// tMin is the farthest "near" intersection (amongst the X,Y and Z planes pairs)
			if (t1 > tMin)
				tMin = t1;

			// And here's the trick :
			// If "far" is closer than "near", then there is NO intersection.
			// See the images in the tutorials for the visual explanation.
			if (tMax < tMin)
				return false;

		}
		else { // Rare case : the ray is almost parallel to the planes, so they don't have any "intersection"
			if (-e + aabb_min.x > 0.0f || -e + aabb_max.x < 0.0f)
				return false;
		}
	}


	// Test intersection with the 2 planes perpendicular to the OBB's Y axis
	// Exactly the same thing than above.
	{
		glm::vec3 yaxis(ModelMatrix[1].x, ModelMatrix[1].y, ModelMatrix[1].z);
		float e = glm::dot(yaxis, delta);
		float f = glm::dot(ray_direction, yaxis);

		if (fabs(f) > 0.001f) {

			float t1 = (e + aabb_min.y) / f;
			float t2 = (e + aabb_max.y) / f;

			if (t1 > t2) { float w = t1; t1 = t2; t2 = w; }

			if (t2 < tMax)
				tMax = t2;
			if (t1 > tMin)
				tMin = t1;
			if (tMin > tMax)
				return false;

		}
		else {
			if (-e + aabb_min.y > 0.0f || -e + aabb_max.y < 0.0f)
				return false;
		}
	}


	// Test intersection with the 2 planes perpendicular to the OBB's Z axis
	// Exactly the same thing than above.
	{
		glm::vec3 zaxis(ModelMatrix[2].x, ModelMatrix[2].y, ModelMatrix[2].z);
		float e = glm::dot(zaxis, delta);
		float f = glm::dot(ray_direction, zaxis);

		if (fabs(f) > 0.001f) {

			float t1 = (e + aabb_min.z) / f;
			float t2 = (e + aabb_max.z) / f;

			if (t1 > t2) { float w = t1; t1 = t2; t2 = w; }

			if (t2 < tMax)
				tMax = t2;
			if (t1 > tMin)
				tMin = t1;
			if (tMin > tMax)
				return false;

		}
		else {
			if (-e + aabb_min.z > 0.0f || -e + aabb_max.z < 0.0f)
				return false;
		}
	}

	intersection_distance = tMin;
	return true;
}


static void error_callback(int error, const char *description)
{
	cerr << description << endl;
}

static void key_callback(GLFWwindow *window, int key, int scancode, int action, int mods)
{

	float speed = .5;
	vec3 w = normalize(lookAtPoint);
	w *= -1;
	vec3 u = normalize(cross(upV, w));

	if (key == GLFW_KEY_ESCAPE && action == GLFW_PRESS) {
		glfwSetWindowShouldClose(window, GL_TRUE);
	}
	else if (key == GLFW_KEY_W) {
		eye -= speed * w;
	} else if (key == GLFW_KEY_A) {
		eye -= speed * u;
	} else if (key == GLFW_KEY_S) {
		eye += speed * w;
	} else if (key == GLFW_KEY_D) {
		eye += speed * u;
	} else if (key == GLFW_KEY_SPACE && action == GLFW_PRESS) {
		if (timePause) {
			timePause = 0;
		} else {
			timePause++;
		}
	} else if (key == GLFW_KEY_RIGHT && action == GLFW_PRESS) {
		timeForward = true;
	} else if (key == GLFW_KEY_LEFT && action == GLFW_PRESS) {
		timeForward = false;
	}
	else if (key == GLFW_KEY_UP && action == GLFW_PRESS) {
		timeSpeed *= 10;
	}
	else if (key == GLFW_KEY_DOWN && action == GLFW_PRESS) {
		timeSpeed /= 10;
	}
	else if (key == GLFW_KEY_J && action == GLFW_PRESS) {
      cTheta += 5;
   } else if (key == GLFW_KEY_L && action == GLFW_PRESS) {
      cTheta -= 5;
   } else if (key == GLFW_KEY_I && action == GLFW_PRESS) {
      cHeight += .5;
   } else if (key == GLFW_KEY_K && action == GLFW_PRESS) {
      cHeight -= 0.5;
   }
}

void mouse_callback(GLFWwindow* window, double xpos, double ypos) {
	if (firstMove) {
		lastX = xpos;
		lastY = ypos;
		firstMove = false;
	}
	float xoffset = xpos - lastX;
	float yoffset = lastY - ypos;
	lastX = xpos;
	lastY = ypos;

	float radius = 3;
	float sens = .05;

	xoffset *= sens;
	yoffset *= sens;

	alpha += yoffset;
	beta += xoffset;


	if (alpha > 80.0f)
		alpha = 80.0f;
	if (alpha < -80.0f)
		alpha = -80.0f;

	vec3 newLook;
	newLook.x = radius * cos(radians(alpha)) * cos(radians(beta));
	newLook.y = radius * sin(radians(alpha));
	newLook.z = radius * cos(radians(alpha)) * cos(radians(90.0 - beta));
	lookAtPoint = normalize(newLook);
}




float p2wx(double in_x, float left) {
	float c = (-2*left)/(g_width-1.0);
	float d = left;
   return c*in_x+d;
}

float p2wy(double in_y, float bot) {
	//flip y
  	in_y = g_height -in_y;
	float e = (-2*bot)/(g_height-1.0);
	float f = bot;
   return e*in_y + f;
}

static void resize_callback(GLFWwindow *window, int width, int height) {
	g_width = width;
	g_height = height;
	glViewport(0, 0, width, height);
}

/* code to define the ground plane */
static void initGeom() {

   float g_groundSize = 20;
   float g_groundY = -1.5;

  // A x-z plane at y = g_groundY of dimension [-g_groundSize, g_groundSize]^2
    float GrndPos[] = {
    -g_groundSize, g_groundY, -g_groundSize,
    -g_groundSize, g_groundY,  g_groundSize,
     g_groundSize, g_groundY,  g_groundSize,
     g_groundSize, g_groundY, -g_groundSize
    };

    float GrndNorm[] = {
     0, 1, 0,
     0, 1, 0,
     0, 1, 0,
     0, 1, 0,
     0, 1, 0,
     0, 1, 0
    };

  static GLfloat GrndTex[] = {
      0, 0, // back
      0, 15,
      15, 15,
      15, 0 };

   unsigned short idx[] = {0, 1, 2, 0, 2, 3};

   GLuint VertexArrayID;
	//generate the VAO
	glGenVertexArrays(1, &VertexArrayID);
	glBindVertexArray(VertexArrayID);

    g_GiboLen = 6;
    glGenBuffers(1, &GrndBuffObj);
    glBindBuffer(GL_ARRAY_BUFFER, GrndBuffObj);
    glBufferData(GL_ARRAY_BUFFER, sizeof(GrndPos), GrndPos, GL_STATIC_DRAW);

    glGenBuffers(1, &GrndNorBuffObj);
    glBindBuffer(GL_ARRAY_BUFFER, GrndNorBuffObj);
    glBufferData(GL_ARRAY_BUFFER, sizeof(GrndNorm), GrndNorm, GL_STATIC_DRAW);
    
	 glGenBuffers(1, &GrndTexBuffObj);
    glBindBuffer(GL_ARRAY_BUFFER, GrndTexBuffObj);
    glBufferData(GL_ARRAY_BUFFER, sizeof(GrndTex), GrndTex, GL_STATIC_DRAW);

    glGenBuffers(1, &GIndxBuffObj);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, GIndxBuffObj);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(idx), idx, GL_STATIC_DRAW);


	for (int i = 0; i<100; i++) {
		positions[i] = glm::vec3(rand() % 20 - 10, rand() % 20, rand() % 20 - 10);
		//orientations[i] = glm::quat(glm::vec3(rand() % 360, rand() % 360, rand() % 360));
		explosions[i] = 0;

	}
}

static void init()
{
	GLSL::checkVersion();

	// Set background color.
	glClearColor(0.5f, 0.5f, 1.0f, 1.0f);
	// Enable z-buffer test.
	glEnable(GL_DEPTH_TEST);

  // Initialize mesh.
   shape = make_shared<Shape>();
   shape->loadMesh(RESOURCE_DIR + "dog.obj");
   shape->resize();
   shape->init();
   
	world = make_shared<Shape>();
   world->loadMesh(RESOURCE_DIR + "sphere.obj");
   world->resize();
   world->init();

	// Initialize the GLSL programs
	prog0 = make_shared<Program>();
	prog0->setVerbose(true);
	prog0->setShaderNames(RESOURCE_DIR + "tex_vert.glsl", RESOURCE_DIR + "simple_geom.glsl", RESOURCE_DIR + "tex_frag0.glsl");
	prog0->init(true);
	
	prog1 = make_shared<Program>();
	prog1->setVerbose(true);
	prog1->setShaderNames(RESOURCE_DIR + "tex_vert.glsl", RESOURCE_DIR + "simple_geom.glsl", RESOURCE_DIR + "tex_frag1.glsl");
	prog1->init(true);
  
	prog2 = make_shared<Program>();
	prog2->setVerbose(true);
	prog2->setShaderNames(RESOURCE_DIR + "tex_vert.glsl", RESOURCE_DIR + "tex_frag2.glsl");
	prog2->init(false);
	
	//////////////////////////////////////////////////////
   // Intialize textures
   //////////////////////////////////////////////////////
	texture0 = make_shared<Texture>();
   texture0->setFilename(RESOURCE_DIR + "fur.jpg");
   texture0->init();
   texture0->setUnit(0);
   texture0->setWrapModes(GL_CLAMP_TO_EDGE, GL_CLAMP_TO_EDGE);

   texture1 = make_shared<Texture>();
   texture1->setFilename(RESOURCE_DIR + "world.jpg");
   texture1->init();
   texture1->setUnit(1);
   texture1->setWrapModes(GL_CLAMP_TO_EDGE, GL_CLAMP_TO_EDGE);
   
	texture2 = make_shared<Texture>();
   texture2->setFilename(RESOURCE_DIR + "grass.jpg");
   texture2->init();
   texture2->setUnit(2);
   texture2->setWrapModes(GL_REPEAT, GL_REPEAT);

	/// Add uniform and attributes to each of the programs
	prog0->addUniform("P");
	prog0->addUniform("M");
	prog0->addUniform("V");
	prog0->addAttribute("vertPos");
   prog0->addAttribute("vertNor");
	prog0->addAttribute("vertTex");
	prog0->addUniform("myTime");
	prog0->addUniform("timePause");
	prog0->addUniform("explosion");
   prog0->addUniform("Texture0");
	
	prog1->addUniform("P");
	prog1->addUniform("M");
	prog1->addUniform("V");
	prog1->addAttribute("vertPos");
   prog1->addAttribute("vertNor");
	prog1->addAttribute("vertTex");
	prog1->addUniform("myTime");
	prog1->addUniform("timePause");
	prog1->addUniform("explosion");
   prog1->addUniform("Texture1");
	
	prog2->addUniform("P");
	prog2->addUniform("M");
	prog2->addUniform("V");
	prog2->addAttribute("vertPos");
   prog2->addAttribute("vertNor");
	prog2->addAttribute("vertTex");
   prog2->addUniform("Texture2");

}


/****DRAW
This is the most important function in your program - this is where you
will actually issue the commands to draw any geometry you have set up to
draw
********/
static void render()
{
	// Get current frame buffer size.
	int width, height;
	glfwGetFramebufferSize(window, &width, &height);
	float aspect = width/(float)height;
	glViewport(0, 0, width, height);

	// Clear framebuffer.
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);


	// Create the matrix stacks 
	auto P = make_shared<MatrixStack>();
	auto M = make_shared<MatrixStack>();
	mat4 V = glm::lookAt(eye, eye + lookAtPoint, upV);

	P->pushMatrix();
   P->perspective(45.0f, aspect, 0.01f, 100.0f);

	//draw the dog mesh 

	prog0->bind();
	texture0->bind(prog0->getUniform("Texture0"));
    glUniform1f(prog0->getUniform("myTime"), myTime);
	glUniformMatrix4fv(prog0->getUniform("P"), 1, GL_FALSE, value_ptr(P->topMatrix()));
	glUniformMatrix4fv(prog0->getUniform("V"), 1, GL_FALSE, value_ptr(V));

	M->pushMatrix();
	M->loadIdentity();
	
	for (int i = 0; i < 100; i++) {
		M->pushMatrix();
		M->translate(positions[i]);
		glUniform1f(prog0->getUniform("explosion"), explosions[i]);
		glUniformMatrix4fv(prog0->getUniform("M"), 1, GL_FALSE, value_ptr(M->topMatrix()));
		shape->draw(prog0);
		M->popMatrix();
	}


	prog0->unbind();

	//draw the world sphere	
	prog1->bind();
   texture1->bind(prog1->getUniform("Texture1"));
   glUniform1f(prog0->getUniform("myTime"), myTime);
	glUniformMatrix4fv(prog0->getUniform("P"), 1, GL_FALSE, value_ptr(P->topMatrix()));
	glUniformMatrix4fv(prog0->getUniform("V"), 1, GL_FALSE, value_ptr(V));
	
	M->pushMatrix();
	M->translate(vec3(1, 0, -5));
	glUniform1f(prog0->getUniform("myTime"), myTime);
	glUniform1i(prog0->getUniform("timePause"), timePause);
	glUniformMatrix4fv(prog0->getUniform("M"), 1, GL_FALSE, value_ptr(M->topMatrix()));
	glUniformMatrix4fv(prog0->getUniform("V"), 1, GL_FALSE, value_ptr(V));
    
	world->draw(prog1);
	
	glDisableVertexAttribArray(0);
	glDisableVertexAttribArray(1);
	glDisableVertexAttribArray(2);
	M->popMatrix();
	prog1->unbind();

	//draw the ground plane	
	prog2->bind();
	M->pushMatrix();
   texture2->bind(prog2->getUniform("Texture2"));
   glUniform1f(prog0->getUniform("myTime"), myTime);
   glUniform1i(prog0->getUniform("timePause"), timePause);
	glUniformMatrix4fv(prog0->getUniform("P"), 1, GL_FALSE, value_ptr(P->topMatrix()));
	glUniformMatrix4fv(prog0->getUniform("V"), 1, GL_FALSE, value_ptr(V));
	glUniformMatrix4fv(prog0->getUniform("M"), 1, GL_FALSE, value_ptr(M->topMatrix()));

	glEnableVertexAttribArray(0);
   glBindBuffer(GL_ARRAY_BUFFER, GrndBuffObj);
   glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, 0);

	glEnableVertexAttribArray(1);
   glBindBuffer(GL_ARRAY_BUFFER, GrndNorBuffObj);
   glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 0, 0);
	 
	glEnableVertexAttribArray(2);
   glBindBuffer(GL_ARRAY_BUFFER, GrndTexBuffObj);
   glVertexAttribPointer(2, 2, GL_FLOAT, GL_FALSE, 0, 0);

   // draw!
   glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, GIndxBuffObj);
   glDrawElements(GL_TRIANGLES, g_GiboLen, GL_UNSIGNED_SHORT, 0);

	glDisableVertexAttribArray(0);
	glDisableVertexAttribArray(1);
	glDisableVertexAttribArray(2);
	
	M->popMatrix();
	prog2->unbind();

	P->popMatrix();

	if (!timePause) {
		for (int i = 0; i < 100; i++) {
			if (explosions[i] > 0) {
				if (timeForward) {
					explosions[i] += timeSpeed;
				}
				else {
					explosions[i] -= timeSpeed;
				}
			}
		}
		
	}

	if (glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_LEFT)){
			
			glm::vec3 ray_origin;
			glm::vec3 ray_direction;
			ScreenPosToWorldRay(
				1024/2, 768/2,
				1024, 768, 
				V, 
				P->topMatrix(),
				ray_origin, 
				ray_direction
			);	
			
			//ray_direction = ray_direction*20.0f;

			message = "background";

			// Test each each Oriented Bounding Box (OBB).
			// A physics engine can be much smarter than this, 
			// because it already has some spatial partitionning structure, 
			// like Binary Space Partitionning Tree (BSP-Tree),
			// Bounding Volume Hierarchy (BVH) or other.
			for(int i=0; i<100; i++){

				float intersection_distance; // Output of TestRayOBBIntersection()
				glm::vec3 aabb_min(-1.0f, -1.0f, -1.0f);
				glm::vec3 aabb_max( 1.0f,  1.0f,  1.0f);

				// The ModelMatrix transforms :
				// - the mesh to its desired position and orientation
				// - but also the AABB (defined with aabb_min and aabb_max) into an OBB
				//mat4 RotationMatrix = toMat4(orientations[i]);
				mat4 TranslationMatrix = translate(mat4(), positions[i]);
				mat4 ModelMatrix = TranslationMatrix;


				if (TestRayOBBIntersection(
					ray_origin, 
					ray_direction, 
					aabb_min, 
					aabb_max,
					ModelMatrix,
					intersection_distance)
				){
					std::ostringstream oss;
					oss << "mesh " << i;
					message = oss.str();
					cout << message << endl;
					cout << "this sucks" << endl;
					explosions[i] += timeSpeed;
					break;
				}
			}


		}


}

int main(int argc, char **argv)
{

	g_width = 1920;
	g_height = 1080;
	/* we will always need to load external shaders to set up where */
	if(argc < 2) {
      cout << "Please specify the resource directory." << endl;
      return 0;
   }
   RESOURCE_DIR = argv[1] + string("/");

	/* your main will always include a similar set up to establish your window
      and GL context, etc. */

	// Set error callback as openGL will report errors but they need a call back
	glfwSetErrorCallback(error_callback);
	// Initialize the library.
	if(!glfwInit()) {
		return -1;
	}
	//request the highest possible version of OGL - important for mac
	
	glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
	glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 2);

	// Create a windowed mode window and its OpenGL context.
	window = glfwCreateWindow(g_width, g_height, "textures", NULL, NULL);
	if(!window) {
		glfwTerminate();
		return -1;
	}
	// Make the window's context current.
	glfwMakeContextCurrent(window);
	// Initialize GLEW.
	glewExperimental = true;
	if(glewInit() != GLEW_OK) {
		cerr << "Failed to initialize GLEW" << endl;
		return -1;
	}

	glGetError();
	cout << "OpenGL version: " << glGetString(GL_VERSION) << endl;
	cout << "GLSL version: " << glGetString(GL_SHADING_LANGUAGE_VERSION) << endl;

	// Set vsync.
	glfwSwapInterval(1);
	// Set keyboard callback.
	glfwSetKeyCallback(window, key_callback);
	//set the window resize call back
	glfwSetFramebufferSizeCallback(window, resize_callback);

	glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_DISABLED);
	glfwSetCursorPosCallback(window, mouse_callback);

	/* This is the code that will likely change program to program as you
		may need to initialize or set up different data and state */
	// Initialize scene.
	initGeom();
	cout << "done initializing geometry" << endl;
	init();
	cout << "done initializing shaders" << endl;

	// Loop until the user closes the window.
	while(!glfwWindowShouldClose(window)) {
		// Render scene.
		render();
		// Swap front and back buffers.
		glfwSwapBuffers(window);
		// Poll for and process events.
		glfwPollEvents();
	}
	// Quit program.
	glfwDestroyWindow(window);
	glfwTerminate();
	return 0;
}
