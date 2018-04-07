#include <getopt.h>
#include <iostream>
#include <fstream>
#include <nanogui/nanogui.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <unordered_set>

#include "CGL/CGL.h"
#include "collision/plane.h"
#include "collision/sphere.h"
#include "fluid.h"
#include "fluidSimulator.h"
#include "json.hpp"

typedef uint32_t gid_t;

using namespace std;
using namespace nanogui;

using json = nlohmann::json;

#define msg(s) cerr << "[FluidSim] " << s << endl;

const string SPHERE = "sphere";
const string PLANE = "plane";
const string PARTICLE = "particle";
const string FLUID = "fluid";

const unordered_set<string> VALID_KEYS = {SPHERE, PLANE, PARTICLE, FLUID};

FluidSimulator *app = nullptr;
GLFWwindow *window = nullptr;
Screen *screen = nullptr;

void error_callback(int error, const char* description) {
  puts(description);
}

void createGLContexts() {
  if (!glfwInit()) {
    return;
  }

  glfwSetTime(0);

  glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
  glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
  glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
  glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

  glfwWindowHint(GLFW_SAMPLES, 0);
  glfwWindowHint(GLFW_RED_BITS, 8);
  glfwWindowHint(GLFW_GREEN_BITS, 8);
  glfwWindowHint(GLFW_BLUE_BITS, 8);
  glfwWindowHint(GLFW_ALPHA_BITS, 8);
  glfwWindowHint(GLFW_STENCIL_BITS, 8);
  glfwWindowHint(GLFW_DEPTH_BITS, 24);
  glfwWindowHint(GLFW_RESIZABLE, GL_TRUE);

  // Create a GLFWwindow object
  window = glfwCreateWindow(800, 800, "Fluid Simulator", nullptr, nullptr);
  if (window == nullptr) {
    std::cout << "Failed to create GLFW window" << std::endl;
    glfwTerminate();
    return;
  }
  glfwMakeContextCurrent(window);

  if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress)) {
    throw std::runtime_error("Could not initialize GLAD!");
  }
  glGetError(); // pull and ignore unhandled errors like GL_INVALID_ENUM

  glClearColor(0.2f, 0.25f, 0.3f, 1.0f);
  glClear(GL_COLOR_BUFFER_BIT);

  // Create a nanogui screen and pass the glfw pointer to initialize
  screen = new Screen();
  screen->initialize(window, true);

  int width, height;
  glfwGetFramebufferSize(window, &width, &height);
  glViewport(0, 0, width, height);
  glfwSwapInterval(1);
  glfwSwapBuffers(window);
}

void setGLFWCallbacks() {
  glfwSetCursorPosCallback(window, [](GLFWwindow *, double x, double y) {
    if (!screen->cursorPosCallbackEvent(x, y)) {
      app->cursorPosCallbackEvent(x / screen->pixelRatio(),
                                  y / screen->pixelRatio());
    }
  });

  glfwSetMouseButtonCallback(
      window, [](GLFWwindow *, int button, int action, int modifiers) {
        if (!screen->mouseButtonCallbackEvent(button, action, modifiers) ||
            action == GLFW_RELEASE) {
          app->mouseButtonCallbackEvent(button, action, modifiers);
        }
      });

  glfwSetKeyCallback(
      window, [](GLFWwindow *, int key, int scancode, int action, int mods) {
        if (!screen->keyCallbackEvent(key, scancode, action, mods)) {
          app->keyCallbackEvent(key, scancode, action, mods);
        }
      });

  glfwSetCharCallback(window, [](GLFWwindow *, unsigned int codepoint) {
    screen->charCallbackEvent(codepoint);
  });

  glfwSetDropCallback(window,
                      [](GLFWwindow *, int count, const char **filenames) {
                        screen->dropCallbackEvent(count, filenames);
                        app->dropCallbackEvent(count, filenames);
                      });

  glfwSetScrollCallback(window, [](GLFWwindow *, double x, double y) {
    if (!screen->scrollCallbackEvent(x, y)) {
      app->scrollCallbackEvent(x, y);
    }
  });

  glfwSetFramebufferSizeCallback(window,
                                 [](GLFWwindow *, int width, int height) {
                                   screen->resizeCallbackEvent(width, height);
                                   app->resizeCallbackEvent(width, height);
                                 });
}

void usageError(const char *binaryName) {
  printf("Usage: %s [options]\n", binaryName);
  printf("Required program options:\n");
  printf("  -f     <STRING>    Filename of scene");
  printf("\n");
  exit(-1);
}

void incompleteObjectError(const char *object, const char *attribute) {
  cout << "Incomplete " << object << " definition, missing " << attribute << endl;
  exit(-1);
}

void loadObjectsFromFile(string filename, Fluid *fluid, FluidParameters *cp, vector<CollisionObject *>* objects) {
  // Read JSON from file
  ifstream i(filename);
  json j;
  i >> j;

  // Loop over objects in scene
  for (json::iterator it = j.begin(); it != j.end(); ++it) {
    string key = it.key();

    // Check that object is valid
    unordered_set<string>::const_iterator query = VALID_KEYS.find(key);
    if (query == VALID_KEYS.end()) {
      cout << "Invalid scene object found: " << key << endl;
      exit(-1);
    }

    // Retrieve object
    json object = it.value();

    // Parse object depending on type (fluid, sphere, or plane)

    if (key == PARTICLE) {
      Vector3D origin;
      double radius, friction;

      auto it_origin = object.find("origin");
      if (it_origin != object.end()) {
        vector<double> vec_origin = *it_origin;
        origin = Vector3D(vec_origin[0], vec_origin[1], vec_origin[2]);
      } else {
        incompleteObjectError("particle", "origin");
      }

      auto it_radius = object.find("radius");
      if (it_radius != object.end()) {
        radius = *it_radius;
      } else {
        incompleteObjectError("particle", "radius");
      }

      auto it_friction = object.find("friction");
      if (it_friction != object.end()) {
        friction = *it_friction;
      } else {
        incompleteObjectError("particle", "friction");
      }

      Particle *p = new Particle(origin, radius, friction);
      objects->push_back(p);
    } else if (key == PLANE) { // PLANE
      Vector3D point, normal;
      double friction;

      auto it_point = object.find("point");
      if (it_point != object.end()) {
        vector<double> vec_point = *it_point;
        point = Vector3D(vec_point[0], vec_point[1], vec_point[2]);
      } else {
        incompleteObjectError("plane", "point");
      }

      auto it_normal = object.find("normal");
      if (it_normal != object.end()) {
        vector<double> vec_normal = *it_normal;
        normal = Vector3D(vec_normal[0], vec_normal[1], vec_normal[2]);
      } else {
        incompleteObjectError("plane", "normal");
      }

      auto it_friction = object.find("friction");
      if (it_friction != object.end()) {
        friction = *it_friction;
      } else {
        incompleteObjectError("plane", "friction");
      }

      Plane *p = new Plane(point, normal, friction);
      objects->push_back(p);
    } else if (key == FLUID){
      int num_width_points, num_height_points, num_length_points;
      double width, height, length;
      
      auto it_num_width_points = object.find("num_width_points");
      if (it_num_width_points != object.end()) {
        num_width_points = *it_num_width_points;
      } else {
        incompleteObjectError("fluid", "num_width_points");
      }

      auto it_num_height_points = object.find("num_height_points");
      if (it_num_height_points != object.end()) {
        num_height_points = *it_num_height_points;
      } else {
        incompleteObjectError("fluid", "num_height_points");
      }
      
      auto it_num_length_points = object.find("num_length_points");
      if (it_num_length_points != object.end()) {
        num_length_points = *it_num_length_points;
      } else {
        incompleteObjectError("fluid", "num_length_points");
      }
      
      auto it_width = object.find("width");
      if (it_width != object.end()) {
        width = *it_width;
      } else {
        incompleteObjectError("fluid", "width");
      }

      auto it_height = object.find("height");
      if (it_height != object.end()) {
        height = *it_height;
      } else {
        incompleteObjectError("fluid", "height");
      }
      
      auto it_length = object.find("length");
      if (it_length != object.end()) {
        length = *it_length;
      } else {
        incompleteObjectError("fluid", "length");
      }
      
      fluid->num_width_points = num_width_points;
      fluid->num_height_points = num_height_points;
      fluid->num_length_points = num_length_points;

      fluid->width = width;
      fluid->height = height;
      fluid->length = length;
    }
  }

  i.close();
}

int main(int argc, char **argv) {
  Fluid fluid;
  FluidParameters fp;
  vector<CollisionObject *> objects;

  if (argc == 1) { // No arguments, default initialization
    // string default_file_name = "../scene/pinned2.json";
    // loadObjectsFromFile(default_file_name, &fluid, &fp, &objects);
  } else {
    int c;

    while ((c = getopt (argc, argv, "f:")) != -1) {
      switch (c) {
        case 'f':
          loadObjectsFromFile(optarg, &fluid, &fp, &objects);
          break;
        default:
          usageError(argv[0]);
      }
    }
  }

  glfwSetErrorCallback(error_callback);

  createGLContexts();

  // Initialize the Fluid object
  fluid.buildGrid();

  // Initialize the FluidSimulator object
  app = new FluidSimulator(screen);
  app->loadFluid(&fluid);
  app->loadFluidParameters(&fp);
  app->loadCollisionObjects(&objects);
  app->init();

  // Call this after all the widgets have been defined

  screen->setVisible(true);
  screen->performLayout();

  // Attach callbacks to the GLFW window

  setGLFWCallbacks();

  while (!glfwWindowShouldClose(window)) {
    glfwPollEvents();

    glClearColor(0.25f, 0.25f, 0.25f, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    app->drawContents();

    // Draw nanogui
    screen->drawContents();
    screen->drawWidgets();

    glfwSwapBuffers(window);

    if (!app->isAlive()) {
      glfwSetWindowShouldClose(window, 1);
    }
  }

  return 0;
}
