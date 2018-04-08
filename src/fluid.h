#ifndef FLUID_H
#define FLUID_H

#include <unordered_set>
#include <unordered_map>
#include <vector>

#include "CGL/CGL.h"
#include "CGL/misc.h"
#include "collision/collisionObject.h"
#include "collision/particle.h"

using namespace CGL;
using namespace std;

enum e_orientation { HORIZONTAL = 0, VERTICAL = 1 };

struct FluidParameters {
  FluidParameters() {}
  FluidParameters(double damping,
                  double density, double ks)
      : damping(damping), density(density), ks(ks) {}
  ~FluidParameters() {}

  // Global simulation parameters

  double damping;

  // Mass-spring parameters
  double density;
  double ks;
};

struct Fluid {
  Fluid() {}
  Fluid(double width, double length, double height, double particle_radius,
        int num_particles, int num_height_points, 
        int num_width_points, int num_length_points);
  ~Fluid();

  void buildGrid();
  GLfloat* getBuffer();
  void simulate(double frames_per_sec, double simulation_steps, FluidParameters *fp,
                vector<Vector3D> external_accelerations,
                vector<CollisionObject *> *collision_objects);

  void reset();

  // Fluid properties
  double width;
  double length;
  double height;
  double particle_radius;
  int num_width_points;
  int num_length_points;
  int num_height_points;
  
  double radius;
  double friction;
  
  // Used to find neighboring particles
  double R=0.1;

  // Fluid components
  vector<Particle> particles;
  
  // Spatial hashing
  unordered_map<float, vector<Particle *> *> map;
  
  void build_spatial_map();
  float hash_position(Vector3D pos, int xOffset=0, int yOffset=0, int zOffset=0);
  
  std::vector<Particle *> getNeighbors(Vector3D pos);
};

#endif /* FLUID_H */
