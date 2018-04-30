#ifndef COLLISIONOBJECT
#define COLLISIONOBJECT

#include <nanogui/nanogui.h>
// #include "particle.h"

using namespace CGL;
using namespace std;
using namespace nanogui;

class Particle;
class CollisionObject {
public:
  virtual void render(GLShader &shader)=0;
  virtual void collide_particle(Particle &pm)=0;
  double friction;
};

#endif /* COLLISIONOBJECT */
