#include "iostream"
#include <nanogui/nanogui.h>

#include "../fluidSimulator.h"
#include "plane.h"

using namespace std;
using namespace CGL;

#define SURFACE_OFFSET 0.0001

void Plane::collide_particle(Particle &pm) {
  // TODO (Part 3.2): Handle collisions with planes.
  float last_side = dot(pm.last_origin-point, normal);
  float current_side = dot(pm.origin-point, normal);
  //std::cout << "++" << std::endl;
  //std::cout << last_side << std::endl;
  //std::cout << current_side << std::endl;
  //if (abs(current_side) < pm.radius ) {
  if ( ((dot(point - pm.origin, normal) > 0.0) && !(dot(point - pm.last_origin, normal) > 0.0))){
    double t = dot((point - pm.origin), normal) / dot(normal, normal);
    Vector3D intersect_point = pm.origin + t * normal;
    Vector3D correction_vector = intersect_point - pm.last_origin + normal * SURFACE_OFFSET;
    pm.origin = pm.last_origin + correction_vector * (1 - friction);
  } else if ((dot(point - pm.origin, normal) < 0.0) && !(dot(point - pm.last_origin, normal) < 0.0)){
    // Reverse all the signs of normal vector.
    double t = dot((point - pm.origin), -normal) / dot(-normal, -normal);
    Vector3D intersect_point = pm.origin + t * -normal;
    Vector3D correction_vector = intersect_point - pm.last_origin + -normal * SURFACE_OFFSET;
    pm.origin = pm.last_origin + correction_vector * (1 - friction);
  }
}

void Plane::render(GLShader &shader) {
  /*nanogui::Color color(0.7f, 0.7f, 0.7f, 1.0f);

  Vector3f sPoint(point.x, point.y, point.z);
  Vector3f sNormal(normal.x, normal.y, normal.z);
  Vector3f sParallel(normal.y - normal.z, normal.z - normal.x,
                     normal.x - normal.y);
  sParallel.normalize();
  Vector3f sCross = sNormal.cross(sParallel);

  MatrixXf positions(3, 4);
  MatrixXf normals(3, 4);

  positions.col(0) << sPoint + 2 * (sCross + sParallel);
  positions.col(1) << sPoint + 2 * (sCross - sParallel);
  positions.col(2) << sPoint + 2 * (-sCross + sParallel);
  positions.col(3) << sPoint + 2 * (-sCross - sParallel);

  normals.col(0) << sNormal;
  normals.col(1) << sNormal;
  normals.col(2) << sNormal;
  normals.col(3) << sNormal;

  if (shader.uniform("in_color", false) != -1) {
    shader.setUniform("in_color", color);
  }
  shader.uploadAttrib("in_position", positions);
  shader.uploadAttrib("in_normal", normals);

  shader.drawArray(GL_TRIANGLE_STRIP, 0, 4);*/
}
