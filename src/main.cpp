#include <cfloat>
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
#include "collision/triangle.h"
#include "collision/sphere.h"
#include "fluid.h"
#include "fluidSimulator.h"
#include "json.hpp"
#include "shader.hpp"
#include "OBJ_Loader.h"


typedef uint32_t gid_t;

using namespace std;
using namespace nanogui;

using json = nlohmann::json;

#define msg(s) cerr << "[FluidSim] " << s << endl;

const string SPHERE = "sphere";
const string PLANE = "plane";
const string PARTICLE = "particle";
const string TRIANGLE = "triangle";
const string OBJECT = "object";
const string FLUID = "fluid";
const string BOUNDINGBOX = "boundingBox";

const unordered_set<string> VALID_KEYS = {SPHERE, PLANE, PARTICLE, TRIANGLE, FLUID, OBJECT, BOUNDINGBOX};

FluidSimulator *app = nullptr;
GLFWwindow *window = nullptr;
Screen *screen = nullptr;

void error_callback(int error, const char* description) {
  puts(description);
}

void createGLContexts() {
  // if (!glfwInit()) {
  //   return;
  // }
  //
  // glfwSetTime(0);
  //
  // glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
  // glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
  // glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
  // glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
  //
  // glfwWindowHint(GLFW_SAMPLES, 0);
  // glfwWindowHint(GLFW_RED_BITS, 8);
  // glfwWindowHint(GLFW_GREEN_BITS, 8);
  // glfwWindowHint(GLFW_BLUE_BITS, 8);
  // glfwWindowHint(GLFW_ALPHA_BITS, 8);
  // glfwWindowHint(GLFW_STENCIL_BITS, 8);
  // glfwWindowHint(GLFW_DEPTH_BITS, 24);
  // glfwWindowHint(GLFW_RESIZABLE, GL_TRUE);

  // Create a GLFWwindow object
  // window = glfwCreateWindow(800, 800, "Fluid Simulator", nullptr, nullptr);
  // if (window == nullptr) {
  //   std::cout << "Failed to create GLFW window" << std::endl;
  //   glfwTerminate();
  //   return;
  // }
  // glfwMakeContextCurrent(window);
  //
  // if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress)) {
  //   throw std::runtime_error("Could not initialize GLAD!");
  // }
  // glGetError(); // pull and ignore unhandled errors like GL_INVALID_ENUM
  //
  // glClearColor(0.2f, 0.25f, 0.3f, 1.0f);
  // glClear(GL_COLOR_BUFFER_BIT);

  // Create a nanogui screen and pass the glfw pointer to initialize
  // screen = new Screen();
  // screen->initialize(window, true);
  //
  // int width, height;
  // glfwGetFramebufferSize(window, &width, &height);
  // glViewport(0, 0, width, height);
  // glfwSwapInterval(1);
  // glfwSwapBuffers(window);
}

void setGLFWCallbacks() {
  // glfwSetCursorPosCallback(window, [](GLFWwindow *, double x, double y) {
  //   if (!screen->cursorPosCallbackEvent(x, y)) {
  //     app->cursorPosCallbackEvent(x / screen->pixelRatio(),
  //                                 y / screen->pixelRatio());
  //   }
  // });
  //
  // glfwSetMouseButtonCallback(
  //     window, [](GLFWwindow *, int button, int action, int modifiers) {
  //       if (!screen->mouseButtonCallbackEvent(button, action, modifiers) ||
  //           action == GLFW_RELEASE) {
  //         app->mouseButtonCallbackEvent(button, action, modifiers);
  //       }
  //     });
  //
  // glfwSetKeyCallback(
  //     window, [](GLFWwindow *, int key, int scancode, int action, int mods) {
  //       if (!screen->keyCallbackEvent(key, scancode, action, mods)) {
  //         app->keyCallbackEvent(key, scancode, action, mods);
  //       }
  //     });
  //
  // glfwSetCharCallback(window, [](GLFWwindow *, unsigned int codepoint) {
  //   screen->charCallbackEvent(codepoint);
  // });
  //
  // glfwSetDropCallback(window,
  //                     [](GLFWwindow *, int count, const char **filenames) {
  //                       screen->dropCallbackEvent(count, filenames);
  //                       app->dropCallbackEvent(count, filenames);
  //                     });
  //
  // glfwSetScrollCallback(window, [](GLFWwindow *, double x, double y) {
  //   if (!screen->scrollCallbackEvent(x, y)) {
  //     app->scrollCallbackEvent(x, y);
  //   }
  // });

  // glfwSetFramebufferSizeCallback(window,
  //                                [](GLFWwindow *, int width, int height) {
  //                                  screen->resizeCallbackEvent(width, height);
  //                                  app->resizeCallbackEvent(width, height);
  //                                });
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

  Vector3D minBoundaries = Vector3D(DBL_MAX, DBL_MAX, DBL_MAX);
  Vector3D maxBoundaries = Vector3D(-DBL_MAX, -DBL_MAX, -DBL_MAX);

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

      /*Particle *p = new Particle(origin, radius, friction);
      objects->push_back(p);*/

      fluid->radius = radius;
      fluid->friction = friction;
    } else if (key == TRIANGLE) {
        Vector3D a,b,c;
        double friction;
        auto it_a = object.find("a");
        auto it_b = object.find("b");
        auto it_c = object.find("c");
        if (it_a != object.end()) {
          vector<double> vec_a = *it_a;
          a = Vector3D(vec_a[0], vec_a[1], vec_a[2]);
        } else {
          incompleteObjectError("Triangle", "Vertex a");
        }
        if (it_b != object.end()) {
          vector<double> vec_b = *it_b;
          b = Vector3D(vec_b[0], vec_b[1], vec_b[2]);
        } else {
          incompleteObjectError("Triangle", "Vertex b");
        }
        if (it_a != object.end()) {
          vector<double> vec_c = *it_c;
          c = Vector3D(vec_c[0], vec_c[1], vec_c[2]);
        } else {
          incompleteObjectError("Triangle", "Vertex c");
        }
        auto it_friction = object.find("friction");
        if (it_friction != object.end()) {
          friction = *it_friction;
        } else {
          incompleteObjectError("particle", "friction");
        }
        Triangle* tri = new Triangle(a,b,c,friction);
        objects->push_back(tri);
    } else if (key == OBJECT) {
      // COLLADA objects
      double object_friction = 0.0;
      auto fric = object.find("friction");
      if (fric != object.end()) object_friction = *fric;

      Vector3D ob_origin;
      auto ob_it_origin = object.find("origin");
      if (ob_it_origin != object.end()) {
        vector<double> vec_origin = *ob_it_origin;
        ob_origin = Vector3D(vec_origin[0], vec_origin[1], vec_origin[2]);
      }

      auto object_file_name = object.find("file");
      if (object_file_name != object.end()) {
        std::vector<Vector3D *> vertices;
        std::vector<Triangle *> tris;
        std::string objfilename = object_file_name->get<std::string>();;

        objl::Loader loader;
        loader.LoadFile(objfilename);

        objl::Mesh mesh = loader.LoadedMeshes[0];
        for (auto vertex: mesh.Vertices) {
          Vector3D *new_v = new Vector3D(vertex.Position.X, vertex.Position.Y, vertex.Position.Z);
          *new_v = *new_v + ob_origin;
          vertices.emplace_back(new_v);
        }
        for (int j = 0; j < mesh.Indices.size(); j+=3) {
          Triangle *new_t = new Triangle(*vertices[j],*vertices[j+1],*vertices[j+2],  object_friction);
          objects->push_back(new_t);
        }
      } else {
        incompleteObjectError("Object", "File");
      }
      //load collada into a list of traingles

      // insert trinagles into objects


    } else if (key == BOUNDINGBOX) { // PLANE
      for (auto plane : object){
        Vector3D point, normal;
        double friction;

        auto it_point = plane.find("point");
        if (it_point != plane.end()) {
          vector<double> vec_point = *it_point;
          point = Vector3D(vec_point[0], vec_point[1], vec_point[2]);
        } else {
          incompleteObjectError("plane", "point");
        }

        if (minBoundaries.x > point.x){
          minBoundaries.x = point.x;
        }
        if (minBoundaries.y > point.y){
          minBoundaries.y = point.y;
        }
        if (minBoundaries.z > point.z){
          minBoundaries.z = point.z;
        }

        if (maxBoundaries.x < point.x){
          maxBoundaries.x = point.x;
        }
        if (maxBoundaries.y < point.y){
          maxBoundaries.y = point.y;
        }
        if (maxBoundaries.z < point.z){
          maxBoundaries.z = point.z;
        }

        auto it_normal = plane.find("normal");
        if (it_normal != plane.end()) {
          vector<double> vec_normal = *it_normal;
          normal = Vector3D(vec_normal[0], vec_normal[1], vec_normal[2]);
        } else {
          incompleteObjectError("plane", "normal");
        }

        auto it_friction = plane.find("friction");
        if (it_friction != plane.end()) {
          friction = *it_friction;
        } else {
          incompleteObjectError("plane", "friction");
        }

        Plane *p = new Plane(point, normal, friction);
        objects->push_back(p);
      }
    } else if (key == FLUID){
      int num_width_points, num_height_points, num_length_points;
      int neighborhood_particle;
      double width, height, length, R;
      int num_width_voxels, num_height_voxels, num_length_voxels;

      auto it_num_width_points = object.find("num_width_points");
      if (it_num_width_points != object.end()) {
        num_width_points = *it_num_width_points;
      } else {
        incompleteObjectError("fluid", "num_width_points");
      }

      auto _neighborhood_particle = object.find("neighborhood_particle");
      if (_neighborhood_particle != object.end()) {
        neighborhood_particle = *_neighborhood_particle;
      } else {
        incompleteObjectError("fluid", "neighborhood_particle");
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

      auto it_R = object.find("r");
      if (it_R != object.end()) {
        R = *it_R;
      } else {
        incompleteObjectError("fluid", "r");
      }

      auto it_num_width_voxels = object.find("num_width_voxels");
      if (it_num_width_voxels != object.end()) {
        num_width_voxels = *it_num_width_voxels;
      } else {
        incompleteObjectError("fluid", "num_width_voxels");
      }

      auto it_num_height_voxels = object.find("num_height_voxels");
      if (it_num_height_voxels != object.end()) {
        num_height_voxels = *it_num_height_voxels;
      } else {
        incompleteObjectError("fluid", "num_height_voxels");
      }

      auto it_num_length_voxels = object.find("num_length_voxels");
      if (it_num_length_voxels != object.end()) {
        num_length_voxels = *it_num_length_voxels;
      } else {
        incompleteObjectError("fluid", "num_length_voxels");
      }


      auto it_rho = object.find("rho_o");
      if (it_rho != object.end()) {
        fluid->RHO_O = *it_rho;
      } else {
        incompleteObjectError("rho_o", "r");
      }

      auto it_si = object.find("solver_iter");
      if (it_si != object.end()) {
        fluid->solver_iters = *it_si;
      } else {
        incompleteObjectError("solver_iter", "r");
      }

      auto it_fps = object.find("fps");
      if (it_fps != object.end()) {
        fluid->fps = *it_fps;
      } else {
        incompleteObjectError("fps", "r");
      }

      auto it_sf = object.find("sf");
      if (it_sf != object.end()) {
        fluid->sf = *it_sf;
      } else {
        incompleteObjectError("sf", "r");
      }

      auto it_vis = object.find("viscosity");
      if (it_vis != object.end()) {
        fluid->viscosity = *it_vis;
      } else {
        incompleteObjectError("viscosity", "r");
      }

      auto it_vor = object.find("vorticity");
      if (it_vor != object.end()) {
        fluid->vorticity = *it_vor;
      } else {
        incompleteObjectError("vorticity", "r");
      }

      auto it_numberCube = object.find("numberCube");
      if (it_numberCube != object.end()) {
        fluid->numberCube = *it_numberCube;
      } else {
        incompleteObjectError("numberCube", "r");
      }

      fluid->maxBoundaries = maxBoundaries + 1.0;
      fluid->minBoundaries = minBoundaries - 1.0;



      fluid->R = R;
      fluid->W_CONSTANT =  315.0/(64.0*PI*pow(R,9));
      fluid->W_DEL_CONSTANT = 45.0/(PI*pow(R,6)); //TODO this maybe negated

      fluid->num_width_points = num_width_points;
      fluid->num_height_points = num_height_points;
      fluid->num_length_points = num_length_points;
      fluid->neighborhood_particle = neighborhood_particle;

      fluid->width = width;
      fluid->height = height;
      fluid->length = length;

      fluid->num_cells = Vector3D(num_width_voxels, num_height_voxels, num_length_voxels);

      fluid->mass = fluid->RHO_O*width*height*length/(num_height_points*num_width_points*num_length_points);
      // fluid -> RHO_O = RHO_O;
      // fluid->solver_iters = si;
      // fluid->fps = fps;
      // fluid->sf = sf;
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

  // glfwSetErrorCallback(error_callback);

  // createGLContexts();

  // Initialize the Fluid object
  fluid.buildGrid();

  // Initialize the FluidSimulator object
  app = new FluidSimulator(nullptr);
  app->loadFluid(&fluid);
  app->loadFluidParameters(&fp);
  app->loadCollisionObjects(&objects);
  app->init();

  // Call this after all the widgets have been defined

  // screen->setVisible(true);
  // screen->performLayout();

  // Attach callbacks to the GLFW window

  setGLFWCallbacks();

  while (true) {
    // glfwPollEvents();
    //
    // glClearColor(0.25f, 0.25f, 0.25f, 1.0f);
    // glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    app->drawContents();

    // Draw nanogui
    // screen->drawContents();
    // screen->drawWidgets();

    // glfwSwapBuffers(window);
    //
    // if (!app->isAlive()) {
    //   glfwSetWindowShouldClose(window, 1);
    // }
  }

  return 0;
}
