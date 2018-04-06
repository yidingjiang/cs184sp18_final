#include <iostream>
#include <math.h>
#include <random>
#include <vector>

#include "fluid.h"
#include "collision/plane.h"
#include "collision/particle.h"

using namespace std;

Fluid::Fluid(double width, double height, int num_width_points,
             int num_height_points, float thickness) {
  this->width = width;
  this->height = height;
  this->num_width_points = num_width_points;
  this->num_height_points = num_height_points;
  this->thickness = thickness;

  buildGrid();
  buildFluidMesh();
}

Fluid::~Fluid() {
  point_masses.clear();
  springs.clear();

  if (clothMesh) {
    delete clothMesh;
  }
}

void Fluid::buildGrid() {
  // TODO (Part 1.1): Build a grid of masses.
  double w_offset = width / ((double) num_width_points);
  double l_offset = length / ((double) num_length_points);
  double h_offset = height / ((double) num_height_points);
  for (int i = 0; i < num_width_points; i++) {
    for (int j = 0; j < num_length_points; j++) {
      for (int k = 0; k < num_height_points; k++) {
        Vector3D pos = Vector3D(i * w_offset, j * l_offset, k * h_offset);
        Particle p = Particle(pos, radius, 0.0);
        particles.emplace_back(p);
      }
    }
  }

}

void Fluid::simulate(double frames_per_sec, double simulation_steps, FluidParameters *fp,
                     vector<Vector3D> external_accelerations,
                     vector<CollisionObject *> *collision_objects) {
  double mass = width * height * fp->density / num_width_points / num_height_points;
  double delta_t = 1.0f / frames_per_sec / simulation_steps;

  // TODO (Part 2.1): Compute total force acting on each point mass.
  for (auto &p : particles) {
    // reseting all forces
    p.forces = Vector3D(0,0,0);
  }

  for (auto ea: external_accelerations){
    for (auto &p: particles) {
      p.forces += mass * ea;
    }
  }

  // TODO (Part 2.2): Use Verlet integration to compute new point mass positions
  for (auto &m : point_masses) {
    if (!m.pinned) {
      Vector3D temp = m.position;
      // std::cout << m.forces << std::endl;
      // std::cout << (1.-cp->damping/100.) * (m.position - m.last_position) + pow(delta_t, 2) * m.forces/mass << std::endl;
      m.position += (1.-cp->damping/100.) * (m.position - m.last_position) + pow(delta_t, 2) * m.forces/mass;
      m.last_position = temp;
      // std::cout << m.position << std::endl;
    }
  }

  // This won't do anything until you complete Part 4.
  // build_spatial_map();
  // for (Particle &pm : point_masses) {
  //   self_collide(pm, simulation_steps);
  // }

  // This won't do anything until you complete Part 3.
  // for (Particle &pm : point_masses) {
  //   for (CollisionObject *co : *collision_objects) {
  //     co->collide(pm);
  //   }
  // }
}

// void Fluid::build_spatial_map() {
//   for (const auto &entry : map) {
//     delete(entry.second);
//   }
//   map.clear();
//   // TODO (Part 4.2): Build a spatial map out of all of the point masses.
//   for (Particle &pm : point_masses) {

//     float k = hash_position(pm.position);

//     if (map.find(k) == map.end()) {
//       // std::cout << "not found" << std::endl;
//       // std::cout << count << std::endl;
//       map[k] = new vector<Particle *>;
//       map[k]->push_back(&pm);
//     } else {
//       // std::cout << "found" << std::endl;
//       map[k]->push_back(&pm);
//     }
//   }

// }

// float Fluid::hash_position(Vector3D pos) {
//   // TODO (Part 4.1): Hash a 3D position into a unique float identifier that represents
//   // membership in some uniquely identified 3D box volume.
//   float w = 3. * width / ((float) num_width_points);
//   float h = 3. * height / ((float) num_height_points);
//   float t = std::max(w, h);
//   float x = pos.x - fmod(pos.x, w);
//   float y = pos.y - fmod(pos.y, h);
//   float z = pos.z - fmod(pos.z, t);
//   return (x * 31 + y) * 31 + z;
//   // return 0.f;
// }

///////////////////////////////////////////////////////
/// YOU DO NOT NEED TO REFER TO ANY CODE BELOW THIS ///
///////////////////////////////////////////////////////

void Fluid::reset() {
  Particle *p = &particles[0];
  for (int i = 0; i < point_masses.size(); i++) {
    pm->position = pm->start_position;
    pm->last_position = pm->start_position;
    pm++;
  }
}
