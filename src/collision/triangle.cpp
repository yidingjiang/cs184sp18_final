#include "iostream"
#include <nanogui/nanogui.h>

#include "../fluidSimulator.h"
#include "triangle.h"
#include "plane.cpp"

using namespace std;
using namespace CGL;

#define SURFACE_OFFSET 1e-6


// bool Triangle::point_in_triangle(Vector3D pt) {
//   std::vector<Vector3D> vertices = {a,b,c};
//   for (int i = 0; i < 3; i++) {
//     int k = (i == 0)? 2 : i-1;
//     Vector3D diff = vertices[i] - vertices[k];
//     Vector3D center_diff = plane->point - vertices[k];
//     Vector3D particle_diff = pt - vertices[k];
//     if (CGL::dot(CGL::cross(particle_diff, diff), CGL::cross(center_diff, diff)) < 0) return false;
//   }
//   return true;
// }

void Triangle::collide_particle(Particle &pm) {
//   // TODO (Part 3.2): Handle collisions with planes.
//   Vector3D point = this->plane->point;
//   Vector3D normal = this->plane->normal;
//   float last_side = dot(point - pm.origin, normal);
//   float current_side = dot(point - pm.x_star, normal);

//   if ((last_side > 0) == (current_side > 0)) return;

//   normal = (current_side > 0.0) ? normal : -normal;
//   double t = current_side/dot(normal, normal);
//   Vector3D intersect_point = pm.x_star + t * normal;

//   if (!point_in_triangle(intersect_point)) return;

//   Vector3D correction_vector = intersect_point - pm.origin + normal * SURFACE_OFFSET;
//   pm.x_star = pm.origin + correction_vector * (1 - friction);
}

void Triangle::render(GLShader &shader) {
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
