#ifndef COLLISIONOBJECT_POLYGON_H
#define COLLISIONOBJECT_POLYGON_H

#include <nanogui/nanogui.h>

#include "plane.h"
#include "particle.h"

using namespace nanogui;
using namespace CGL;
using namespace std;

//Note: Only convex polygons allowed
struct Polygon : public Plane {
public:
  Polygon(Plane* plane, const std::vector<Vector3D> &vertices){
      this->plane = plane;
      this->vertices = vertices;
      for (auto vertex: vertices) center += vertex;
      center /= vertices.size();
  }

  void render(GLShader &shader);
  void collide_particle(Particle &pm);
  bool particle_in_polygon(Particle &pm);

  Plane* plane;
  std::vector<Vector3D> vertices;
  Vector3D center;

  double friction;
};

#endif /* COLLISIONOBJECT_PLANE_H */
