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
        Vector3D pos = Vector3D(i * w_offset + ((double)(rand()%10 - 5))/100.,
                                j * l_offset + ((double)(rand()%10 - 5))/100., 
                                k * h_offset + ((double)(rand()%10 - 5))/100.);
        Particle p = Particle(pos, radius, friction);
        particles.emplace_back(p);
      }
    }
  }

}

void Fluid::build_spatial_map() {
  for (const auto &entry : map) {
    delete(entry.second);
  }
  map.clear();

  for (Particle &particle : this->particles){
    float key = hash_position(particle.origin);
    if (map.find(key) == map.end()){
      map[key] = new std::vector<Particle *>();
    }
    map[key]->push_back(&particle);
  }
}

float Fluid::hash_position(Vector3D pos, int xOffset, int yOffset, int zOffset) {
  // TODO (Part 4.1): Hash a 3D position into a unique float identifier that represents
  // membership in some uniquely identified 3D box volume.
  double w = 3 * width / num_width_points;
  double h = 3 * height / num_height_points;
  double t = max(w,h);
  float xVol = pos.x - fmod(pos.x, w);
  float yVol = pos.y - fmod(pos.y, h);
  float zVol = pos.z - fmod(pos.z, t);
  
  xVol += xOffset;
  yVol += yOffset;
  zVol += zOffset;
  
  return (xVol * 31 + yVol) * 31 + zVol;
}

std::vector<Particle *> Fluid::getNeighbors(Vector3D pos){
  std::vector<Particle *> neighbors = std::vector<Particle *>();
  
  
  if (map.find(hash_position(pos)) == map.end()){
    vector<Particle *> currCell = *map[hash_position(pos)];
    neighbors.insert(neighbors.end(), currCell.begin(), currCell.end());
  }
  
  if (map.find(hash_position(pos, 1, 0, 0)) == map.end()){
    vector<Particle *> currCell = *map[hash_position(pos, 1, 0, 0)];
    neighbors.insert(neighbors.end(), currCell.begin(), currCell.end());
  }
  
  if (map.find(hash_position(pos, -1, 0, 0)) == map.end()){
    vector<Particle *> currCell = *map[hash_position(pos, -1, 0, 0)];
    neighbors.insert(neighbors.end(), currCell.begin(), currCell.end());
  }
  
  if (map.find(hash_position(pos, 0, 1, 0)) == map.end()){
    vector<Particle *> currCell = *map[hash_position(pos, 0, 1, 0)];
    neighbors.insert(neighbors.end(), currCell.begin(), currCell.end());
  }
  
  if (map.find(hash_position(pos, 0, -1, 0)) == map.end()){
    vector<Particle *> currCell = *map[hash_position(pos, 0, -1, 0)];
    neighbors.insert(neighbors.end(), currCell.begin(), currCell.end());
  }
  
  if (map.find(hash_position(pos, 0, 0, 1)) == map.end()){
    vector<Particle *> currCell = *map[hash_position(pos, 0, 0, 1)];
    neighbors.insert(neighbors.end(), currCell.begin(), currCell.end());
  }
  
  if (map.find(hash_position(pos, 0, 0, -1)) == map.end()){
    vector<Particle *> currCell = *map[hash_position(pos, 0, 0, -1)];
    neighbors.insert(neighbors.end(), currCell.begin(), currCell.end());
  }
  
  return neighbors;
}

void Fluid::simulate(double frames_per_sec, double simulation_steps, FluidParameters *fp,
                     vector<Vector3D> external_accelerations,
                     vector<CollisionObject *> *collision_objects) {
  double mass = 0.01;
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
    m.origin += 1.0 * (m.origin-m.last_origin) + pow(delta_t, 2) * m.forces/mass;
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
  
  build_spatial_map();

  // This won't do anything until you complete Part 3.
  for (auto &m : this->particles) {
    for (CollisionObject *co : *collision_objects) {
      co->collide_particle(m);
    }
    for (auto &p : this->particles) {
      if ((p.origin-m.origin).norm() > 1e-10) p.collide_particle(m);
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
