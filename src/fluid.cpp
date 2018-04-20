#include <iostream>
#include <math.h>
#include <random>
#include <vector>

#include "fluid.h"
#include "collision/plane.h"
#include "collision/particle.h"

using namespace std;
// TODO instantiate particles with the correct mass, size, and distances.

#define EPSILON 5000.0

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
        Vector3D pos = Vector3D((i-num_width_points/2.0) * w_offset + width/2.0,
                                j * l_offset,
                                (k - num_height_points/2.0) * h_offset + height/2.0);
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
  std::vector<std::vector<Particle *>>  neighborArray = build_index();

  for(int iter=0; iter<solver_iters; iter++) {
    this->update_lambdas(neighborArray);
    this->update_delta_p(neighborArray);
    // i = 0;
    //apply delta_p and perform collision detection

    for (Particle &p: this->particles) {
      p.x_star += p.delta_p*sf;
      // for (CollisionObject *co : *collision_objects) co->collide_particle(p);
      // i++;
    }
    for (Particle &p: this->particles) {
      for (CollisionObject *co : *collision_objects) co->collide_particle(p);
      // i++;
    }
  }
  i  = 0;
  for (Particle &p: this->particles) {
    p.velocity = (p.x_star-p.origin)/delta_t;
     //vorticity confinement
    //TODO don't know what to do with the corrective force.
    // viscosity

    // this->apply_viscosity(p);
    // p.last_origin = p.origin;
    p.origin = p.x_star;

    //update color
    // p.color = Vector3D( neighborArray[i].size()/10 , 1, 1);
    p.color = Vector3D( p.density/RHO_O ,0,(int)(p.origin.y <= -0.19));
    i++;
  }
  // int nidx = (int) rand()*1000.0/RAND_MAX;
  // for(int k = 0; k < neighborArray[nidx].size(); k++) neighborArray[nidx][k]->color = Vector3D( 1, 1, 1);
  // cout << nidx << endl;

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
/*
* Simulation Physics Code
*/

double Fluid::W(Vector3D r) { //density kernel
  double r_norm = r.norm();
  // std::cout << "norm " << r_norm << '\n';
  if (r_norm > R) return 0;
  // std::cout << "W CONSTANT " << W_CONSTANT << '\n';
  // std::cout << "R " << R << '\n';
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


// double Fluid::density(Particle p, std::vector<Particle *>  neighbors) {
//   Vector3D pi = p.x_star;
//   double r = 0;
//   for (Particle * &pj: neighbors) r += W(pi-pj->x_star);
//   p.density = r*mass; //TODO maybe include mass?
//   return r*mass;
// }

double Fluid::density(Particle p, std::vector<Particle *>  neighbors) {
  Vector3D pi = p.x_star;
  double r = 0;
  // cout << neighbors.size() << endl;
  for (Particle * &pj: neighbors) {
    // std::cout << "W " << W(pi-pj->x_star) << '\n';
    r += W(pi-pj->x_star);
  }
  p.density = r; //TODO maybe include mass?
  return r;
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



// void Fluid::update_lambdas(std::vector<std::vector<Particle *>>  neighborArray) {
//   int i = 0;
//   double accum;
//   for (Particle &p: this-> particles) {
//
//     std::vector<Particle *> neighbors  = neighborArray[i];
//     double ci = density(p, neighbors)/RHO_O - 1.0;
//     cout << ci << endl;
//     double sum_sq_norm;
//
//     double del_j=0;
//     Vector3D del_i(0.0,0.0,0.0);
//     Vector3D del_temp(0.0,0.0,0.0);
//
//     for (Particle * &j: neighbors){
//       del_temp = del_W(p.x_star-j->x_star)/RHO_O;
//       del_i += del_temp;
//       del_j += del_temp.norm2();
//     }
//     sum_sq_norm = del_j + del_i.norm2();
//     // cout << "deli " << del_i.norm2() << " delj " << del_j << endl;
//     p.lambda = -ci/(sum_sq_norm+EPSILON);
//     accum += ci;
//     // cout << "Lambda: " << p.lambda << endl;
//     i++;
//   }
//   // cout << accum/particles.size() << endl;
// }

void Fluid::update_lambdas(std::vector<std::vector<Particle *>>  neighborArray) {
  int i = 0;
  for (Particle &p: this-> particles) {

    std::vector<Particle *> neighbors  = neighborArray[i];
    // std::cout << neighbors.size() << '\n';
    // std::cout << density(p, neighbors) << '\n';
    double ci = density(p, neighbors)/RHO_O - 1.0;
    // std::cout << ci << '\n';
    double sum_sq_norm;

    double del_j=0;
    Vector3D del_i(0.0,0.0,0.0);
    Vector3D del_temp(0.0,0.0,0.0);

    for (Particle * &j: neighbors){
      del_temp = del_W(p.x_star-j->x_star)/RHO_O;
      del_i += del_temp;
      del_j += del_temp.norm2();
    }
    sum_sq_norm = del_j + del_i.norm2();
    // std::cout << "divisor " << sum_sq_norm << '\n';
    p.lambda = -ci/(sum_sq_norm+EPSILON);
    // std::cout << "lambda " << p.lambda << '\n';
    i++;
  }
}


// void Fluid::update_delta_p(std::vector<std::vector<Particle *>> neighborArray){
//   int i = 0;
//   for (Particle &p: this-> particles) {
//     Vector3D pi = p.x_star;
//     std::vector<Particle *> neighbors = neighborArray[i];
//     p.delta_p = Vector3D(0.0,0.0,0.0);
//     // //TODO add scorr term
//     // if (neighbors.size() == 1) continue;
//     for (Particle * &pj: neighbors) {
//       double scorr_denom = W(Vector3D(0.1, 0, 0)*R);
//       double scorr = -0.01*pow(W(pi-pj->x_star)/scorr_denom, 4.0);
//       // std::cout << "scorr " << scorr << '\n';
//       // std::cout << "scorr " << p.lambda+pj->lambda+scorr << '\n';
//       p.delta_p += (p.lambda+pj->lambda+scorr)*del_W(pi-pj->x_star);
//     }
//     // std::cout << "scorr " << p.delta_p << '\n';
//     p.delta_p /= RHO_O;
//     // std::cout << "scorr " << p.delta_p << '\n';
//     i++;
//   }
// }

void Fluid::update_delta_p(std::vector<std::vector<Particle *>> neighborArray){
  int i = 0;
  for (Particle &p: this->particles) {
    Vector3D p_pred = p.x_star;
    std::vector<Particle *> neighbors = neighborArray[i];
    Vector3D delta = Vector3D(0.0,0.0,0.0);
    float l = p.lambda;
    for (Particle * &pj: neighbors) {
      float scorr = -0.001 * pow(W(p_pred-pj->x_star)/W(Vector3D(0.01, 0.01, 0.01)*R), 4.0);
      Vector3D gradient = del_W(p_pred-pj->x_star);
      delta += (l + pj->lambda+scorr) * gradient;
    }
    p.delta_p = delta / RHO_O;
    i++;
  }
}



// vec3 delta = vec3(0.0f);
// float l = p->lambda;
// for (int i=0; i<num_neighbors; i++) {
//     float s_corr = -PRESSURE_K * pow(poly6Kernel(p->pred_position, neighbors[i]->pred_position) / poly6Kernel(p->pred_position, vec3(DELTA_Q) + p->pred_position), PRESSURE_N);
//     vec3 gradient = spikyKernelGradient(p->pred_position, neighbors[i]->pred_position);
//
//     delta += (l + neighbors[i]->lambda + s_corr) * gradient;
// }
// p->delta_position = 1.0f / RO_REST * delta;



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


std::vector<std::vector<Particle *>> Fluid::build_index(){
    // Populate cloud
    this->cloud.pts = this->particles;

    kdtree index(3, cloud, KDTreeSingleIndexAdaptorParams(3));
    index.buildIndex();

    std::vector<std::pair<size_t,double> > ret_matches;
    SearchParams params;
    params.sorted = false;

    double query_pt[3] = {0.5,0.5,0.5};

    std::vector<std::vector<Particle *>> to_return;
    for  (int k =0; k < particles.size(); k++){
      std::vector<Particle *> to_append;

      query_pt[0] = particles[k].x_star.x;
      query_pt[1] = particles[k].x_star.y;
      query_pt[2] = particles[k].x_star.z;

      double nMatches = index.radiusSearch(&query_pt[0], R*R, ret_matches, params); //TODO should be R or squared?
      for (auto &pair: ret_matches) {
        to_append.emplace_back(&(this->particles[pair.first]));
      }
      to_return.emplace_back(to_append);
      // cout << to_append.size() << endl;
    }

    return to_return;


}
