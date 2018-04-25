#include <iostream>
#include <math.h>
#include <random>
#include <vector>

#include "fluid.h"
#include "collision/plane.h"
#include "collision/particle.h"

using namespace std;
// TODO instantiate particles with the correct mass, size, and distances.

#define EPSILON 1e-6

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
        Vector3D pos = Vector3D(i * w_offset,
                                j * l_offset,
                                k * h_offset);
        Particle p = Particle(pos, radius, friction);
        particles.emplace_back(p);
      }
    }
  }

}

GLfloat* Fluid::getBuffer() {
    GLfloat* data = (GLfloat*) malloc(sizeof(GLfloat)*this->particles.size()*7);
    int count = 0;
    for (auto particle : particles) {
        data[count * 7] = particle.origin.x;
        data[count * 7+1] = particle.origin.y;
        data[count * 7+2] = particle.origin.z;
        data[count * 7+3] = particle.color.x;
        data[count * 7+4] = particle.color.y;
        data[count * 7+5] = particle.color.z;
        // data[count * 7+3] = particle.origin.x * particle.origin.x;
        // data[count * 7+4] = 1.0f;
        // data[count * 7+5] = particle.origin.z * particle.origin.z;
        data[count * 7+6] = 1.0f;
        count += 1;
    }
    return data;
}

void Fluid::build_spatial_map() {
  for (const auto &entry : map) {
    delete(entry.second);
  }
  map.clear();

  for (Particle &particle : this->particles){
    string key = hash_position(particle.origin);
    if (map.find(key) == map.end()){
      map[key] = new std::vector<Particle *>();
    }
    map[key]->push_back(&particle);
  }
}

string Fluid::hash_position(Vector3D pos, int xOffset, int yOffset, int zOffset) {
  // TODO (Part 4.1): Hash a 3D position into a unique float identifier that represents
  // membership in some uniquely identified 3D box volume.
  double threeR = 3*R;
  int xVol = floor(pos.x / threeR);
  int yVol = floor(pos.y / threeR);
  int zVol = floor(pos.z / threeR);

  if (pos.x < 0){
    xVol = ceil(pos.x / threeR);
  }
  if (pos.y < 0){
    yVol = ceil(pos.y / threeR);
  }
  if (pos.z < 0){
    zVol = ceil(pos.z / threeR);
  }

  xVol += xOffset;
  yVol += yOffset;
  zVol += zOffset;

  return std::to_string(xVol) + "/" + std::to_string(yVol) + "/" +  std::to_string(zVol);
}



std::vector<Particle *> Fluid::getNeighbors(Vector3D pos){
  std::vector<Particle *> neighbors = std::vector<Particle *>();
  // Get the location of all neighboring cells in the hashmap.
  std::vector<string> neighborCellsHashes = std::vector<string>();
  neighborCellsHashes.push_back(hash_position(pos));
  neighborCellsHashes.push_back(hash_position(pos, 1, 0, 0));
  neighborCellsHashes.push_back(hash_position(pos, -1, 0, 0));
  neighborCellsHashes.push_back(hash_position(pos, 0, 1, 0));
  neighborCellsHashes.push_back(hash_position(pos, 0, -1, 0));
  neighborCellsHashes.push_back(hash_position(pos, 0, 0, 1));
  neighborCellsHashes.push_back(hash_position(pos, 0, 0, -1));


  // Iterate through the neighbor cell and check if within R distance.
  for (string neighborCellsHash : neighborCellsHashes){
    if (map[neighborCellsHash] != NULL){
      vector<Particle *> currCell = *map[neighborCellsHash];
      for (Particle* particle : currCell){
        if ((pos-particle->origin).norm() < R) {
          neighbors.push_back(particle);
        }
      }
    }
  }
  return neighbors;
}

void Fluid::simulate(double frames_per_sec, double simulation_steps, FluidParameters *fp,
                     vector<Vector3D> external_accelerations,
                      vector<CollisionObject *> *collision_objects) {
  double delta_t = 1.0f / frames_per_sec / simulation_steps;
  for (auto &p: particles) {
    p.last_origin = p.origin;
    for (auto ea: external_accelerations){
      p.velocity += delta_t*ea;
    }
    p.x_star = p.origin + delta_t*p.velocity;
  }
  int i = 0;
  build_spatial_map();
  std::vector<std::vector<Particle *>>  neighborArray = generateNeighborArray();
  for(int iter=0; iter<solver_iters; iter++) {
    this->update_density(neighborArray);
    this->update_lambdas(neighborArray);
    this->update_delta_p(neighborArray);
    i = 0;
    //apply delta_p and perform collision detection
    for (Particle &p: this->particles) {
      p.x_star += p.delta_p*1e-3;
      // cout << i << endl;
      // cout << "c_i" << p.density/RHO_O - 1 << endl; //<< " delta_p_norm calc: " << p.delta_p.norm2() <<  " " << " Lambda: " << p.lambda << endl;
      // cout << p.delta_p.norm2() << endl;
      // p.delta_p *= 0.0;
      for (CollisionObject *co : *collision_objects) co->collide_particle(p);
      i++;
    }
  }
  i  = 0;
  for (Particle &p: this->particles) {
    p.velocity = (p.x_star-p.origin)/delta_t;
     //vorticity confinement
    // Vector3D corrective_force = f_vorticity(p);
    //TODO don't know what to do with the corrective force.
    // p.velocity += corrective_force*delta_t/mass;
    // viscosity
    // this->apply_viscosity(p);
    p.last_origin = p.origin;
    p.origin = p.x_star;

    //update color
    p.color = Vector3D( (p.origin.x < 0.51) && (p.origin.x > 0), (p.origin.y < 0.51) && (p.origin.y > 0), (p.origin.z < 0.51) && (p.origin.z > 0));
    i++;
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

std::vector<std::vector<Particle *>>  Fluid::generateNeighborArray(){
  std::vector<std::vector<Particle *>> to_return = std::vector<std::vector<Particle *>>();
  for(Particle &p: this->particles) to_return.emplace_back(getNeighbors(p.origin));
  return to_return;
}

/*
* Simulation Physics Code
*/ 

double Fluid::W(Vector3D r) { //density kernel
  double r_norm = r.norm();
  if (r_norm > R) return 0;
  return pow(R*R-r_norm*r_norm, 3.0)*W_CONSTANT;
}

Vector3D Fluid::del_W(Vector3D r) { // gradient of density kernel
  double r_norm = r.norm();
  if (r_norm > R) return Vector3D(0.0,0.0,0.0); //TODO possibly return 0 is r_norm is small
  return -r*pow(R-r_norm,2.0)*W_DEL_CONSTANT/(r_norm + 1e-6); //TODO this could be negated.
}


double Fluid::C_i(Particle p){
  return p.density/RHO_O - 1;
}


void Fluid::update_density(std::vector<std::vector<Particle *>>  neighborArray) {
  // double accum = 0;
  int i = 0;
  for (Particle &p: this-> particles) {
    Vector3D pi = p.origin;
    double r = 0;
    for (Particle * &pj: neighborArray[i]) {
      // cout << neighborArray[i].size() << endl;
      r += W(pi-pj->origin);
    }
    p.density = r; //TODO maybe include mass?
    // accum += particles[i].density;
    i ++;
  }
  // cout << accum << endl;
  // cout << RHO_O << " " << accum/(particles.size()*RHO_O) -1  << endl;
}




Vector3D Fluid::del_ci_j(Particle i, Particle k) {
  Vector3D to_return =  -del_W(i.x_star-k.x_star)/RHO_O; //TODO: could be negated
  return to_return;
}

Vector3D Fluid::del_ci_i(Particle i, std::vector<Particle *> neighbors) {
  Vector3D accum;
  for (Particle * &j : neighbors) {
    accum += del_W(i.x_star-j->x_star); //TODO: might be negated
  }
  Vector3D to_return = accum/RHO_O;
  // cout << "ii:" << to_return << endl;
  return to_return;
}



void Fluid::update_lambdas(std::vector<std::vector<Particle *>>  neighborArray) {
  int i = 0;
  for (Particle &p: this-> particles) {
    double ci = p.density/RHO_O - 1.0;
    // cout << ci << endl;
    std::vector<Particle *> neighbors  = neighborArray[i];
    double sum_sq_norm = del_ci_i(p, neighbors).norm2();
    for (Particle * &j : neighbors) {
      sum_sq_norm += del_ci_j(p, *j).norm2();
    }

    p.lambda = -ci/(sum_sq_norm+EPSILON);
    // cout << "Lambda: " << p.lambda << endl;
    i++;
  }
}


void Fluid::update_delta_p(std::vector<std::vector<Particle *>> neighborArray){
  int i = 0;
  for (Particle &p: this-> particles) {
    Vector3D pi = p.x_star;
    std::vector<Particle *> neighbors = neighborArray[i];
    p.delta_p = Vector3D(0.0,0.0,0.0);
    // //TODO add scorr term
    if (neighbors.size() == 1) continue;
    for (Particle * &pj: neighbors) {
      double scorr = -0.001*pow(W(pi-pj->x_star)/W(Vector3D(0.0,0.0,0.0)), 4.0);
      p.delta_p += (p.lambda+pj->lambda+scorr)*del_W(pi-pj->x_star);
    }
    p.delta_p /= RHO_O;
    i++;
  }
}

// Vector3D Fluid::f_vorticity(Particle p){
//   Vector3D pi = p.origin;
//   Vector3D omega_i = p.omega;
//   double epsilon = EPSILON;
//   Vector3D N;
//   for (Particle * &j: getNeighbors(pi)) {
//     N += (mass/j->density)*(j->omega).norm()*del_W(pi-j->origin);
//   }
//   if (N.norm() > 1e-6) N.normalize();
//   return epsilon*CGL::cross(N, p.omega);
// }

// void Fluid::update_omega(){
//   for (Particle &p: this->particles){
//       Vector3D accum;
//       Vector3D pi = p.origin;
//       for (Particle * &j: getNeighbors(pi)) {
//         Vector3D vij = p.velocity - j->velocity;
//         accum += CGL::cross(vij,del_W(pi-j->origin)); //TODO delW may be negative
//       }
//       p.omega = accum;
//   }
// }


// void Fluid::apply_viscosity(Particle p) {
//   Vector3D accum;
//   Vector3D pi = p.origin;
//   for (Particle * &j: getNeighbors(pi)) {
//     Vector3D vij = p.velocity - j->velocity;
//     accum += W(pi-j->origin)*vij;
//   }
//   p.velocity += 0.01*accum;
// }
