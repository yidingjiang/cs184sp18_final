#ifndef COLLISIONOBJECT_PARTICLE_H
#define COLLISIONOBJECT_PARTICLE_H

#include "CGL/CGL.h"
#include "CGL/misc.h"
#include "CGL/vector3D.h"
#include "collisionObject.h"

using namespace CGL;
using namespace std;

#define RADIUS_OFFSET 1e-8

struct Particle : public CollisionObject {
public:
  Particle(const Vector3D &origin, double radius, double friction)
      : origin(origin), last_origin(origin), start_origin(origin),
      	radius(radius), radius2(radius * radius),
        friction(friction) {}

  void render(nanogui::GLShader &shader);
  void collide_particle(Particle &pm);

  Vector3D origin;
  Vector3D last_origin;
  Vector3D start_origin;
  Vector3D forces;

  double radius;
  double radius2;

  double friction;
};

#endif /* COLLISIONOBJECT_PARTICLE_H */