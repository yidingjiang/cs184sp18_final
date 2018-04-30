#ifndef COLLISIONOBJECT_TRI_H
#define COLLISIONOBJECT_TRI_H

#include <nanogui/nanogui.h>

#include "collisionObject.h"
#include "plane.h"
#include "particle.h"

using namespace nanogui;
using namespace CGL;
using namespace std;

#define SURFACE_OFFSET 1e-6

struct Triangle : public CollisionObject {
public:
  Plane* plane;
  Vector3D a,b,c;
  double friction;
  Triangle(const Vector3D &a, const Vector3D &b, const Vector3D &c, double friction){
    this->a = a;
    this->b = b;
    this->c = c;
    this->friction = friction;
    this->plane = new Plane( (a+b+c)/3.0 , CGL::cross(b-a, c-a), friction);
    }

  void render(GLShader &shader){};

  void collide_particle(Particle &pm){
    Vector3D point = this->plane->point;
    Vector3D normal = this->plane->normal;
    float last_side = dot(point - pm.origin, normal);
    float current_side = dot(point - pm.x_star, normal);

    if ((last_side > 0) == (current_side > 0)) return;

    normal = (current_side > 0.0) ? normal : -normal;
    double t = current_side/dot(normal, normal);
    Vector3D intersect_point = pm.x_star + t * normal;

    if (!point_in_triangle(intersect_point)) return;

    Vector3D correction_vector = intersect_point - pm.origin + normal * SURFACE_OFFSET;
    pm.x_star = pm.origin + correction_vector * (1 - friction);
  };

  bool point_in_triangle(Vector3D pt){
    std::vector<Vector3D> vertices = {a,b,c};
    for (int i = 0; i < 3; i++) {
      int k = (i == 0)? 2 : i-1;
      Vector3D diff = vertices[i] - vertices[k];
      Vector3D center_diff = plane->point - vertices[k];
      Vector3D particle_diff = pt - vertices[k];
      if (CGL::dot(CGL::cross(particle_diff, diff), CGL::cross(center_diff, diff)) < 0) return false;
    }
    return true;
  }
};

#endif /* COLLISIONOBJECT_TRI_H */
