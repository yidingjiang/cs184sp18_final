#ifndef FLUID_H
#define FLUID_H

#include <unordered_set>
#include <unordered_map>
#include <vector>

#include "CGL/CGL.h"
#include "CGL/misc.h"
#include "collision/collisionObject.h"
#include "collision/particle.h"
#include "nanoflann.hpp"
#include "utils.h"

using namespace CGL;
using namespace std;
using namespace nanoflann;

struct vertex {
  Vector3D p;
  Vector3D n;
  vertex(Vector3D pos, Vector3D norm){
    p = pos;
    n = norm;
  };
  vertex(double x, double y, double z){
    p = Vector3D(x,y,z);
    n = Vector3D(0,0,0);
  };
  vertex(Vector3D pos){
    p = pos;
    n = Vector3D(0,0,0);
  };
  vertex(){
    p = Vector3D(0,0,0);
    n = Vector3D(0,0,0);
  };
};

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

  void convertVoxelToFaces(Vector3D min, Vector3D sizeCell);
  void buildGrid();
  GLfloat* getBuffer();
  void simulate(double frames_per_sec, double simulation_steps, FluidParameters *fp,
                vector<Vector3D> external_accelerations,
                vector<CollisionObject *> *collision_objects, int step);

  void reset();
  void saveVoxelsToMitsuba(std::string fileName, Vector3D min, Vector3D max, bool orientation);
  Vector3D gradientNormal(Vector3D pos);

  // Fluid properties
  double width;
  double length;
  double height;
  double particle_radius;
  int num_width_points;
  int num_length_points;
  int num_height_points;
  int neighborhood_particle;
  int solver_iters = 3;

  Vector3D num_cells;
  bool firstFile = true;
  double viscosity;
  double vorticity;

  double radius;
  double friction;

  // Used to find neighboring particles
  double RHO_O = 25000;
  double mass = 1;
  int fps = 60;
  double sf = 1.0;
  
  int numberCube;

  double R=0.15;
  double W_CONSTANT =  315.0/(64.0*M_PI*R*R*R*R*R*R*R*R*R);
  double W_DEL_CONSTANT = 45.0/(M_PI*pow(R,6.0)); //TODO this maybe negated

  Vector3D maxBoundaries;
  Vector3D minBoundaries;

  // Fluid components
  vector<Particle> particles;

  // Spatial hashing
  unordered_map<string, vector<Particle *> *> map;

  // height, width, length
  std::vector<double> voxelGrid;
  std::vector<Vector3D> voxelOrientations;
  vector<vector<vertex>> triangles;

  void saveFacesToObjs(std::string fileName);
  void build_voxel_grid(int frameNum);
  string hash_position(Vector3D pos, int xOffset=0, int yOffset=0, int zOffset=0);
  std::vector<std::vector<Particle *>> generateNeighborArray();

  std::vector<Particle *> getNeighbors(Vector3D pos);

  double W(Vector3D r);
  Vector3D del_W(Vector3D r);
  double C_i(Particle p);
  double density(Particle p, std::vector<Particle *> neighbors);
  void update_delta_p(std::vector<std::vector<Particle *>> neighborArray);
  // double rho_i(Particle p);
  double lambda(Particle i, std::vector<std::vector<Particle *>>  neighborArray);
  void update_lambdas(std::vector<std::vector<Particle *>>  neighborArray);
  Vector3D del_ci_i(Particle i, std::vector<Particle *> neighbors);

  Vector3D del_ci_j(Particle i, Particle k);

  void apply_vorticity(std::vector<std::vector<Particle *>> neighborArray);
  void apply_viscosity(std::vector<std::vector<Particle *>> neighborArray);
  void update_omega(std::vector<std::vector<Particle *>> neighborArray);

  PointCloud cloud;
  typedef KDTreeSingleIndexAdaptor< L2_Simple_Adaptor<double, PointCloud> , PointCloud, 3 > kdtree;
  std::vector<std::vector<Particle *>> build_index();
  std::vector<std::vector<Particle *>> build_nearest_neighbors_index(int numNeighbors);

  double isotropic_kernel(Vector3D pos);

  kdtree *tree = NULL;

  void save_state_to_csv();
  void saveVoxelToCSV(std::string fileName, Vector3D min, Vector3D sizeCell);
};

#endif /* FLUID_H */
