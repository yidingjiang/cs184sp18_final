#ifndef CGL_FLUID_SIMULATOR_H
#define CGL_FLUID_SIMULATOR_H

#include <nanogui/nanogui.h>

#include "camera.h"
#include "fluid.h"
#include "collision/collisionObject.h"

using namespace nanogui;

class FluidSimulator {
public:
  FluidSimulator(Screen *screen);
  ~FluidSimulator();

  void init();

  void loadFluid(Fluid *fluid);
  void loadFluidParameters(FluidParameters *fp);
  void loadCollisionObjects(vector<CollisionObject *> *objects);
  virtual bool isAlive();
  virtual void drawContents();

  // Screen events

  virtual bool cursorPosCallbackEvent(double x, double y);
  virtual bool mouseButtonCallbackEvent(int button, int action, int modifiers);
  virtual bool keyCallbackEvent(int key, int scancode, int action, int mods);
  virtual bool dropCallbackEvent(int count, const char **filenames);
  virtual bool scrollCallbackEvent(double x, double y);
  virtual bool resizeCallbackEvent(int width, int height);

private:
  virtual void initGUI(Screen *screen);
  // void drawParticle(GLShader &shader);
  // void drawWireframe(GLShader &shader);
  // void drawNormals(GLShader &shader);
  void drawPhong(GLShader &shader);

  // Camera methods

  virtual void resetCamera();
  virtual Matrix4f getProjectionMatrix();
  virtual Matrix4f getViewMatrix();

  // Default simulation values

  int frames_per_sec = 90;
  int simulation_steps = 2;

  CGL::Vector3D gravity = CGL::Vector3D(0, -9.8, 0);
  nanogui::Color color = nanogui::Color(1.0f, 0.0f, 0.0f, 1.0f);

  Fluid *fluid;
  FluidParameters *fp;
  vector<CollisionObject *> *collision_objects;

  // OpenGL attributes
  GLuint programID;
  GLfloat* g_vertex_buffer_data;
  GLuint positionsVAO;
  GLuint positionsVBO;

  enum e_shader { PHONG = 0 };
  e_shader activeShader = PHONG;

  vector<GLShader> shaders;

  GLShader phongShader;
  GLShader particleShader;

  // Camera attributes

  CGL::Camera camera;
  CGL::Camera canonicalCamera;

  double view_distance;
  double canonical_view_distance;
  double min_view_distance;
  double max_view_distance;

  double scroll_rate;

  // Screen methods

  Screen *screen;
  void mouseLeftDragged(double x, double y);
  void mouseRightDragged(double x, double y);
  void mouseMoved(double x, double y);

  // Mouse flags

  bool left_down = false;
  bool right_down = false;
  bool middle_down = false;

  // Keyboard flags

  bool ctrl_down = false;

  // Simulation flags

  bool is_paused = false;

  // Screen attributes

  int mouse_x;
  int mouse_y;

  int screen_w;
  int screen_h;

  bool is_alive = true;

  Vector2i default_window_size = Vector2i(1024, 800);
};

#endif // CGL_CLOTH_SIM_H
