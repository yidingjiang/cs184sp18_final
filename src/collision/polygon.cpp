#include "iostream"
#include <nanogui/nanogui.h>

#include "../fluidSimulator.h"
#include "polygon.h"
#include "plane.cpp"

using namespace std;
using namespace CGL;



bool Polygon::particle_in_polygon(Particle &pm) {
  for (int i = 0; i < vertices.size(); i++) {
    int k = i-1;
    if (k < 0) k = vertices.size();
    Vector3D diff = vertices[i] - vertices[k];
    Vector3D center_diff = center - vertices[k];
    Vector3D particle_diff = pm.x_star - vertices[k];
    if (CGL::dot(CGL::cross(particle_diff, diff), CGL::cross(center_diff, diff)) < 0) return false;
  }
  return true;
}

void Polygon::collide_particle(Particle &pm) {
  // TODO (Part 3.2): Handle collisions with planes.
  float last_side = dot(pm.origin-point, normal);
  float current_side = dot(pm.x_star-point, normal);

  if ( ((dot(point - pm.x_star, normal) > 0.0) && !(dot(point - pm.origin, normal) > 0.0) && particle_in_polygon(pm)) ){
    double t = dot((point - pm.x_star), normal) / dot(normal, normal);
    Vector3D intersect_point = pm.x_star + t * normal;
    Vector3D correction_vector = intersect_point - pm.origin + normal * SURFACE_OFFSET;
    pm.x_star = pm.origin + correction_vector * (1 - friction);
  } else if ((dot(point - pm.x_star, normal) < 0.0) && !(dot(point - pm.origin, normal) < 0.0) && particle_in_polygon(pm) ){
    // Reverse all the signs of normal vector.
    double t = dot((point - pm.x_star), -normal) / dot(-normal, -normal);
    Vector3D intersect_point = pm.x_star + t * -normal;
    Vector3D correction_vector = intersect_point - pm.origin + -normal * SURFACE_OFFSET;
    pm.x_star = pm.origin + correction_vector * (1 - friction);
  }
}

void Polygon::render(GLShader &shader) {
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
