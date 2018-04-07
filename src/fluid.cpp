#include <iostream>
#include <math.h>
#include <random>
#include <vector>

#include "fluid.h"
#include "collision/plane.h"
#include "collision/particle.h"

using namespace std;

Fluid::Fluid(double width, double length, double height, double particle_radius,
             int num_particles, int num_height_points, 
             int num_width_points, int num_length_points) {
  this->width = width;
  this->height = height;
  this->num_width_points = num_width_points;
  this->num_height_points = num_height_points;

  buildGrid();
}

Fluid::~Fluid() {
  particles.clear();
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
        Particle p = Particle(pos, radius, friction);
        particles.emplace_back(p);
      }
    }
  }

}

void Fluid::simulate(double frames_per_sec, double simulation_steps, FluidParameters *fp,
                     vector<Vector3D> external_accelerations,
                     vector<CollisionObject *> *collision_objects) {
  double mass = 0.1;
  double delta_t = 1.0f / frames_per_sec / simulation_steps;

  // TODO (Part 2.1): Compute total force acting on each point mass.
  for (auto &p : this->particles) {
    // reseting all forces
    p.forces = Vector3D(0,0,0);
  }

  for (auto ea: external_accelerations){
    for (auto &p: particles) {
      p.forces += mass * ea;
    }
  }
  
  for (auto &m : this->particles) {
    Vector3D temp = m.origin;
    m.origin += m.origin-m.last_origin + pow(delta_t, 2) * m.forces/mass;
    m.last_origin = temp;
  }

  // TODO (Part 2.2): Use Verlet integration to compute new point mass origins
  /*for (auto &m : particles) {
    Vector3D temp = m.origin;
    m.origin += (1.-fp->damping/100.) * (m.origin - m.last_origin) + pow(delta_t, 2) * m.forces/mass;
    m.last_origin = temp;
  }*/

  // This won't do anything until you complete Part 4.
  // build_spatial_map();
  // for (Particle &pm : point_masses) {
  //   self_collide(pm, simulation_steps);
  // }

  // This won't do anything until you complete Part 3.
  for (auto &m : this->particles) {
    for (CollisionObject *co : *collision_objects) {
      co->collide_particle(m);
    }
  }
}

///////////////////////////////////////////////////////
/// YOU DO NOT NEED TO REFER TO ANY CODE BELOW THIS ///
///////////////////////////////////////////////////////

void Fluid::reset() {
  Particle *pm = &particles[0];
  for (int i = 0; i < particles.size(); i++) {
    pm->origin = pm->start_origin;
    pm->last_origin = pm->start_origin;
    pm++;
  }
}
