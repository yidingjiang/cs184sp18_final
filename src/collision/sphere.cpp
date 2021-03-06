#include <nanogui/nanogui.h>

#include "../misc/sphere_drawing.h"
#include "sphere.h"

using namespace nanogui;
using namespace CGL;

void Sphere::collide(Particle &pm) {
  // TODO (Part 3.1): Handle collisions with spheres.
	Vector3D d = pm.origin - origin;
	if (d.norm() <= radius) {
		Vector3D d_normal = d / d.norm();
		Vector3D tangent_p = origin + radius * d_normal;
		Vector3D correction_vec = tangent_p - pm.last_origin;
		pm.origin = pm.last_origin + (1.-friction)*correction_vec;
	}
}

void Sphere::render(GLShader &shader) {
  // We decrease the radius here so flat triangles don't behave strangely
  // and intersect with the sphere when rendered
  Misc::draw_sphere(shader, origin, radius * 0.92);
}
