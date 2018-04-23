#include <cmath>
#include <glad/glad.h>

#include <CGL/vector3D.h>
#include <nanogui/nanogui.h>

#include "fluidSimulator.h"

#include "camera.h"
#include "fluid.h"
#include "collision/plane.h"
#include "collision/particle.h"
#include "misc/camera_info.h"
#include "shader.hpp"

#include <glm/glm.hpp>

using namespace glm;
using namespace nanogui;
using namespace std;

FluidSimulator::FluidSimulator(Screen *screen) {
  this->screen = screen;

  // Initialize OpenGL buffers and shaders

  phongShader.initFromFiles("Phong", "../shaders/camera.vert",
                            "../shaders/phong.frag");
  // particleShader.initFromFiles("Particle", "../shaders/SimpleVertexShader.vertexshader")

  shaders.push_back(phongShader);
  // shaders.push_back(particleShader);

  glEnable(GL_PROGRAM_POINT_SIZE);
  glEnable(GL_DEPTH_TEST);
}

FluidSimulator::~FluidSimulator() {
  for (auto shader : shaders) {
    shader.free();
  }

  if (fluid) delete fluid;
  if (fp) delete fp;
  if (collision_objects) delete collision_objects;
}

void FluidSimulator::loadFluid(Fluid *fluid) {
  this->fluid = fluid;

  g_vertex_buffer_data = fluid->getBuffer();

  glGenVertexArrays(1, &positionsVAO);
  glGenBuffers(1, &positionsVBO);

  glBindVertexArray(positionsVAO);
  // Get a handle (ID) to the vertex buffer object (VBO), a buffer that holds the data that will be transferred to the GPU.
  // We bind the VBO to the global GL_ARRAY_BUFFER, and store our particles into the buffer. GL_STREAM_DRAW is selected
  // since we expect to update the vertex positions every frame.
  glBindBuffer(GL_ARRAY_BUFFER, positionsVBO);
  glBufferData(GL_ARRAY_BUFFER, sizeof(GLfloat) * fluid->particles.size() * 7, g_vertex_buffer_data, GL_STREAM_DRAW);

  // Enable gl_PointSize in the vertex shader to specify the size of a point
  glEnable(GL_VERTEX_PROGRAM_POINT_SIZE);
  // glEnable(GL_VERTEX_PROGRAM_POINT_SIZE);
  glEnable(GL_BLEND);
  glVertexAttribPointer(
    0,                  // vertex positions attribute specified at location 0 in the vertex shader.
    3,                  // position is a vec3
    GL_FLOAT,           // type
    GL_FALSE,           // do not normalize
    7 * sizeof(GLfloat),// Position is 3 floats, then color is 4 floats
    (GLvoid*)0          // first vector value starts at index 0 -> offset = 0
  );
  glEnableVertexAttribArray(0);

  glVertexAttribPointer(
    1,
    4,
    GL_FLOAT,
    GL_FALSE,
    7 * sizeof(GLfloat),
    (GLvoid*)(3* sizeof(GLfloat))
  );
  glEnableVertexAttribArray(1);

  // Unbind the VAO
  glBindVertexArray(0);
  glBindBuffer(GL_ARRAY_BUFFER, 0);
}

void FluidSimulator::loadFluidParameters(FluidParameters *fp) { this->fp = fp; }

void FluidSimulator::loadCollisionObjects(vector<CollisionObject *> *objects) { this->collision_objects = objects; }

/**
 * Initializes the fluid simulation and spawns a new thread to separate
 * rendering from simulation.
 */
void FluidSimulator::init() {
  // Initialize GUI
  initGUI(screen);
  screen->setSize(default_window_size);

  // Initialize camera

  CGL::Collada::CameraInfo camera_info;
  camera_info.hFov = 50;
  camera_info.vFov = 35;
  camera_info.nClip = 0.01;
  camera_info.fClip = 10000;

  // Try to intelligently figure out the camera target

  Vector3D avg_p_position(0, 0, 0);

  for (auto &p : fluid->particles) {
    avg_p_position += p.origin;
  }

  avg_p_position /= fluid->particles.size();
  std::cout << fluid->width << std::endl;
  // CGL::Vector3D target(avg_p_position.x, avg_p_position.y / 2,
  //                      avg_p_position.z);
  CGL::Vector3D target(0.5, 0.2, 0.5);
  CGL::Vector3D c_dir(0., 0., 0.);
  // canonical_view_distance = max(fluid->width, fluid->height) * 0.9;
  canonical_view_distance = 0.9;
  scroll_rate = canonical_view_distance / 10;

  view_distance = canonical_view_distance * 2;
  min_view_distance = canonical_view_distance / 10.0;
  max_view_distance = canonical_view_distance * 20.0;

  // canonicalCamera is a copy used for view resets

  camera.place(target, acos(c_dir.y), atan2(c_dir.x, c_dir.z), view_distance,
               min_view_distance, max_view_distance);
  canonicalCamera.place(target, acos(c_dir.y), atan2(c_dir.x, c_dir.z),
                        view_distance, min_view_distance, max_view_distance);

  screen_w = default_window_size(0);
  screen_h = default_window_size(1);

  camera.configure(camera_info, screen_w, screen_h);
  canonicalCamera.configure(camera_info, screen_w, screen_h);
}

bool FluidSimulator::isAlive() { return is_alive; }

mat4 convertToMat4(Matrix4f m) {
  return mat4(vec4(m(0,0), m(0,1), m(0,2), m(0,3)),
              vec4(m(1,0), m(1,1), m(1,2), m(1,3)),
              vec4(m(2,0), m(2,1), m(2,2), m(2,3)),
              vec4(m(3,0), m(3,1), m(3,2), m(3,3)));
}

void FluidSimulator::drawContents() {
  glEnable(GL_DEPTH_TEST);

  if (!is_paused) {
    vector<Vector3D> external_accelerations = {gravity};

    for (int i = 0; i < simulation_steps; i++) {
      fluid->simulate(frames_per_sec, simulation_steps, fp, external_accelerations, collision_objects);
    }
  }

  // instancing projection
  // GLint particleSizeLocation = glGetUniformLocation(programID, "particle_size");
  // glUniform1i(particleSizeLocation, 3*((float) camera.r)*5/45.0f);
  //
  // Matrix4f _view = getViewMatrix();
  // // convert CLG to GLSL
  // mat4 cameraTransform = convertToMat4(_view);
  //
  // Matrix4f _projection = getProjectionMatrix();
  // mat4 projection = convertToMat4(_projection);
  //
  // mat4 scaling = mat4(vec4(((float) camera.r)/45.0f, 0, 0, 0),
  //                     vec4(0, ((float) camera.r)/45.0f, 0, 0),
  //                     vec4(0, 0, ((float) camera.r)/45.0f, 0),
  //                     vec4(0, 0, 0, 1));
  //
  // GLuint transformLoc = glGetUniformLocation(programID, "transform");
  // GLint projLoc = glGetUniformLocation(programID, "projection");
  // GLint scalingLoc = glGetUniformLocation(programID, "scaling");
  // glUniformMatrix4fv(transformLoc, 1, GL_FALSE, &cameraTransform[0][0]);
  // glUniformMatrix4fv(projLoc, 1, GL_FALSE, &projection[0][0]);
  // glUniformMatrix4fv(scalingLoc, 1, GL_FALSE, &scaling[0][0]);
  //
  // glBindVertexArray(positionsVAO);
  // glBindBuffer(GL_ARRAY_BUFFER, positionsVBO);
  // g_vertex_buffer_data = fluid->getBuffer();
  // glBufferData(GL_ARRAY_BUFFER, sizeof(GLfloat) * fluid->particles.size() * 7, g_vertex_buffer_data, GL_STREAM_DRAW);
  // glDrawArrays(GL_POINTS, 0, fluid->particles.size());
  // glBindVertexArray(0);


  // new shit has come to light
  GLShader shader = shaders[activeShader];
  shader.bind();

  Matrix4f model;
  model.setIdentity();

  Matrix4f view = getViewMatrix();
  Matrix4f projection = getProjectionMatrix();

  Matrix4f viewProjection = projection * view;

  shader.setUniform("model", model);
  shader.setUniform("viewProjection", viewProjection);
  shader.setUniform("light", Vector3f(0.5, 2, 2));
  shader.setUniform("in_color", color);
  // shader.setUniform("particle_size", 6.*(7.-camera.r));
  shader.setUniform("particle_size", 80./camera.r);

  glBindVertexArray(positionsVAO);
  glBindBuffer(GL_ARRAY_BUFFER, positionsVBO);
  g_vertex_buffer_data = fluid->getBuffer();
  glBufferData(GL_ARRAY_BUFFER, sizeof(GLfloat) * fluid->particles.size() * 7, g_vertex_buffer_data, GL_STREAM_DRAW);
  glDrawArrays(GL_POINTS, 0, fluid->particles.size());
  glBindVertexArray(0);
  glBindBuffer(GL_ARRAY_BUFFER, 0);
  free(g_vertex_buffer_data);

  // is_paused = true;
}


// ----------------------------------------------------------------------------
// CAMERA CALCULATIONS
//
// OpenGL 3.1 deprecated the fixed pipeline, so we lose a lot of useful OpenGL
// functions that have to be recreated here.
// ----------------------------------------------------------------------------

void FluidSimulator::resetCamera() { camera.copy_placement(canonicalCamera); }

Matrix4f FluidSimulator::getProjectionMatrix() {
  Matrix4f perspective;
  perspective.setZero();

  double near = camera.near_clip();
  double far = camera.far_clip();

  double theta = camera.v_fov() * M_PI / 360;
  double range = far - near;
  double invtan = 1. / tanf(theta);

  perspective(0, 0) = invtan / camera.aspect_ratio();
  perspective(1, 1) = invtan;
  perspective(2, 2) = -(near + far) / range;
  perspective(3, 2) = -1;
  perspective(2, 3) = -2 * near * far / range;
  perspective(3, 3) = 0;

  return perspective;
}

Matrix4f FluidSimulator::getViewMatrix() {
  Matrix4f lookAt;
  Matrix3f R;

  lookAt.setZero();

  // Convert CGL vectors to Eigen vectors
  // TODO: Find a better way to do this!

  CGL::Vector3D c_pos = camera.position();
  CGL::Vector3D c_udir = camera.up_dir();
  CGL::Vector3D c_target = camera.view_point();

  Vector3f eye(c_pos.x, c_pos.y, c_pos.z);
  Vector3f up(c_udir.x, c_udir.y, c_udir.z);
  Vector3f target(c_target.x, c_target.y, c_target.z);

  R.col(2) = (eye - target).normalized();
  R.col(0) = up.cross(R.col(2)).normalized();
  R.col(1) = R.col(2).cross(R.col(0));

  lookAt.topLeftCorner<3, 3>() = R.transpose();
  lookAt.topRightCorner<3, 1>() = -R.transpose() * eye;
  lookAt(3, 3) = 1.0f;

  return lookAt;
}

// ----------------------------------------------------------------------------
// EVENT HANDLING
// ----------------------------------------------------------------------------

bool FluidSimulator::cursorPosCallbackEvent(double x, double y) {
  if (left_down && !middle_down && !right_down) {
    if (ctrl_down) {
      mouseRightDragged(x, y);
    } else {
      mouseLeftDragged(x, y);
    }
  } else if (!left_down && !middle_down && right_down) {
    mouseRightDragged(x, y);
  } else if (!left_down && !middle_down && !right_down) {
    mouseMoved(x, y);
  }

  mouse_x = x;
  mouse_y = y;

  return true;
}

bool FluidSimulator::mouseButtonCallbackEvent(int button, int action,
                                              int modifiers) {
  switch (action) {
  case GLFW_PRESS:
    switch (button) {
    case GLFW_MOUSE_BUTTON_LEFT:
      left_down = true;
      break;
    case GLFW_MOUSE_BUTTON_MIDDLE:
      middle_down = true;
      break;
    case GLFW_MOUSE_BUTTON_RIGHT:
      right_down = true;
      break;
    }
    return true;

  case GLFW_RELEASE:
    switch (button) {
    case GLFW_MOUSE_BUTTON_LEFT:
      left_down = false;
      break;
    case GLFW_MOUSE_BUTTON_MIDDLE:
      middle_down = false;
      break;
    case GLFW_MOUSE_BUTTON_RIGHT:
      right_down = false;
      break;
    }
    return true;
  }

  return false;
}

void FluidSimulator::mouseMoved(double x, double y) { y = screen_h - y; }

void FluidSimulator::mouseLeftDragged(double x, double y) {
  float dx = x - mouse_x;
  float dy = y - mouse_y;

  camera.rotate_by(-dy * (PI / screen_h), -dx * (PI / screen_w));
}

void FluidSimulator::mouseRightDragged(double x, double y) {
  camera.move_by(mouse_x - x, y - mouse_y, canonical_view_distance);
}

bool FluidSimulator::keyCallbackEvent(int key, int scancode, int action,
                                      int mods) {
  ctrl_down = (bool)(mods & GLFW_MOD_CONTROL);

  if (action == GLFW_PRESS) {
    switch (key) {
    case GLFW_KEY_ESCAPE:
      is_alive = false;
      break;
    case 'r':
    case 'R':
      fluid->reset();
      break;
    case ' ':
      resetCamera();
      break;
    case 'p':
    case 'P':
      is_paused = !is_paused;
      break;
    case 'n':
    case 'N':
      if (is_paused) {
        is_paused = false;
        drawContents();
        is_paused = true;
      }
      break;
    }
  }

  return true;
}

bool FluidSimulator::dropCallbackEvent(int count, const char **filenames) {
  return true;
}

bool FluidSimulator::scrollCallbackEvent(double x, double y) {
  camera.move_forward(y * scroll_rate);
  return true;
}

bool FluidSimulator::resizeCallbackEvent(int width, int height) {
  screen_w = width;
  screen_h = height;

  camera.set_screen_size(screen_w, screen_h);
  return true;
}

void FluidSimulator::initGUI(Screen *screen) {
  Window *window = new Window(screen, "Settings");
  window->setPosition(Vector2i(15, 15));
  window->setLayout(new GroupLayout(15, 6, 14, 5));

  // Simulation constants

  new Label(window, "Simulation", "sans-bold");

  {
    Widget *panel = new Widget(window);
    GridLayout *layout =
        new GridLayout(Orientation::Horizontal, 2, Alignment::Middle, 5, 5);
    layout->setColAlignment({Alignment::Maximum, Alignment::Fill});
    layout->setSpacing(0, 10);
    panel->setLayout(layout);

    new Label(panel, "frames/s :", "sans-bold");

    IntBox<int> *fsec = new IntBox<int>(panel);
    fsec->setEditable(true);
    fsec->setFixedSize(Vector2i(100, 20));
    fsec->setFontSize(14);
    fsec->setValue(frames_per_sec);
    fsec->setSpinnable(true);
    fsec->setCallback([this](int value) { frames_per_sec = value; });

    new Label(panel, "steps/frame :", "sans-bold");

    IntBox<int> *num_steps = new IntBox<int>(panel);
    num_steps->setEditable(true);
    num_steps->setFixedSize(Vector2i(100, 20));
    num_steps->setFontSize(14);
    num_steps->setValue(simulation_steps);
    num_steps->setSpinnable(true);
    num_steps->setMinValue(0);
    num_steps->setCallback([this](int value) { simulation_steps = value; });
  }

  // Damping slider and textbox

  new Label(window, "Damping", "sans-bold");

  {
    Widget *panel = new Widget(window);
    panel->setLayout(
        new BoxLayout(Orientation::Horizontal, Alignment::Middle, 0, 5));

    Slider *slider = new Slider(panel);
    slider->setValue(fp->damping);
    slider->setFixedWidth(105);

    TextBox *percentage = new TextBox(panel);
    percentage->setFixedWidth(75);
    percentage->setValue(to_string(fp->damping));
    percentage->setUnits("%");
    percentage->setFontSize(14);

    slider->setCallback([percentage](float value) {
      percentage->setValue(std::to_string(value));
    });
    slider->setFinalCallback([&](float value) {
      fp->damping = (double)value;
      // cout << "Final slider value: " << (int)(value * 100) << endl;
    });
  }

  // Gravity

  new Label(window, "Gravity", "sans-bold");

  {
    Widget *panel = new Widget(window);
    GridLayout *layout =
        new GridLayout(Orientation::Horizontal, 2, Alignment::Middle, 5, 5);
    layout->setColAlignment({Alignment::Maximum, Alignment::Fill});
    layout->setSpacing(0, 10);
    panel->setLayout(layout);

    new Label(panel, "x :", "sans-bold");

    FloatBox<double> *fb = new FloatBox<double>(panel);
    fb->setEditable(true);
    fb->setFixedSize(Vector2i(100, 20));
    fb->setFontSize(14);
    fb->setValue(gravity.x);
    fb->setUnits("m/s^2");
    fb->setSpinnable(true);
    fb->setCallback([this](float value) { gravity.x = value; });

    new Label(panel, "y :", "sans-bold");

    fb = new FloatBox<double>(panel);
    fb->setEditable(true);
    fb->setFixedSize(Vector2i(100, 20));
    fb->setFontSize(14);
    fb->setValue(gravity.y);
    fb->setUnits("m/s^2");
    fb->setSpinnable(true);
    fb->setCallback([this](float value) { gravity.y = value; });

    new Label(panel, "z :", "sans-bold");

    fb = new FloatBox<double>(panel);
    fb->setEditable(true);
    fb->setFixedSize(Vector2i(100, 20));
    fb->setFontSize(14);
    fb->setValue(gravity.z);
    fb->setUnits("m/s^2");
    fb->setSpinnable(true);
    fb->setCallback([this](float value) { gravity.z = value; });
  }

  // Appearance

  // new Label(window, "Appearance", "sans-bold");

  // {
  //   ComboBox *cb = new ComboBox(window, {"Wireframe", "Normals", "Shaded"});
  //   cb->setFontSize(14);
  //   cb->setCallback(
  //       [this, screen](int idx) { activeShader = static_cast<e_shader>(idx); });
  // }

  // new Label(window, "Color", "sans-bold");
  //
  // {
  //   ColorWheel *cw = new ColorWheel(window, color);
  //   cw->setCallback(
  //       [this](const nanogui::Color &color) { this->color = color; });
  // }
}
