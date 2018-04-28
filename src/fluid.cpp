#include <fstream>
#include <iostream>
#include <string>
#include <math.h>
#include <random>
#include <vector>

#include "fluid.h"
#include "collision/plane.h"
#include "collision/particle.h"
#include "float.h"
#include "cube.cpp"

using namespace std;
// TODO instantiate particles with the correct mass, size, and distances.

#define EPSILON 5000.0
// #define C_VISCOSITY 0.0001
// #define C_VORTICITY 0.005


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
        Vector3D pos = Vector3D((i- (num_width_points/2.0)) * w_offset,
                                j * l_offset,
                                (k - (num_height_points/2.0)) * h_offset);
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
        data[count * 7+6] = 0.5f;
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

void Fluid::build_voxel_grid(int frameNum) {
  vector<double> voxels((num_cells.x) * (num_cells.y) * (num_cells.z), 0);

  this->voxelGrid = voxels;
  Vector3D min = Vector3D(-2, -2, -2);
  Vector3D max = Vector3D(2, 2, 2); //TODO: In future, adjust thee based on scene params


  Vector3D sizeGrid = Vector3D(max.x - min.x, max.y - min.y, max.z - min.z);
  Vector3D sizeCell = Vector3D(sizeGrid.x / num_cells.x, sizeGrid.y / num_cells.y, sizeGrid.z / num_cells.z);


  for (Particle &particle : this->particles){
    Vector3D newWithMinCoord = particle.origin - min;

    double positionArrayX = newWithMinCoord.x / sizeCell.x;
    double positionArrayY = newWithMinCoord.y / sizeCell.y;
    double positionArrayZ = newWithMinCoord.z / sizeCell.z;

    positionArrayX = floor(positionArrayX);
    positionArrayY = floor(positionArrayY);
    positionArrayZ = floor(positionArrayZ);

    Vector3D cellNum = Vector3D(positionArrayX, positionArrayY, positionArrayZ);

    this->voxelGrid[cellNum.x + (num_cells.x) * (cellNum.y + (num_cells.y) * (cellNum.z))] = isotropic_kernel(particle.origin);
    
  }
  
  convertVoxelToFaces(min, sizeCell);
  saveFacesToObjs(std::to_string(frameNum));
}

Vector3D Fluid::gradientNormal(Vector3D pos){
  Vector3D normal = Vector3D(isotropic_kernel(Vector3D(pos.x + 1, pos.y, pos.z)) - isotropic_kernel(Vector3D(pos.x - 1, pos.y, pos.z)), 
                     isotropic_kernel(Vector3D(pos.x, pos.y + 1, pos.z)) - isotropic_kernel(Vector3D(pos.x, pos.y - 1, pos.z)), 
                     isotropic_kernel(Vector3D(pos.x, pos.y, pos.z + 1)) - isotropic_kernel(Vector3D(pos.x, pos.y, pos.z - 1))
                   );
  normal.normalize();
  return normal;
}

void Fluid::convertVoxelToFaces(Vector3D min, Vector3D sizeCell){
  //mcubes(	double start_x, double start_y, double start_z, double end_x, double end_y, double end_z,
	//		double step_x, double step_y, double step_z);
  //cube * test = new cube();

  triangles = vector<vector<vertex>>();
  // DO ORIENTATION STUFF

  //std::cout << "KKK" << std::endl;
  for (int xpos = 0; xpos < num_cells.x; ++xpos) {
    for (int ypos = 0; ypos < num_cells.y; ++ypos) {
      for (int zpos = 0; zpos < num_cells.z; ++zpos) {
        
        if ((xpos < num_cells.x-1) && (ypos < num_cells.y-1) && (zpos < num_cells.z-1)) {
          vector<double> grid = vector<double>();
          grid.push_back(this->voxelGrid[xpos + (num_cells.x) * (ypos + (num_cells.y) * zpos)]);
          grid.push_back(this->voxelGrid[xpos+1 + (num_cells.x) * (ypos + (num_cells.y) * zpos)]);
          grid.push_back(this->voxelGrid[xpos+1 + (num_cells.x) * (ypos + (num_cells.y) * (zpos+1))]);
          grid.push_back(this->voxelGrid[xpos + (num_cells.x) * (ypos + (num_cells.y) * (zpos+1))]);
          grid.push_back(this->voxelGrid[xpos + (num_cells.x) * (ypos+1 + (num_cells.y) * zpos)]);
          grid.push_back(this->voxelGrid[xpos+1 + (num_cells.x) * (ypos+1 + (num_cells.y) * zpos)]);
          grid.push_back(this->voxelGrid[xpos+1 + (num_cells.x) * (ypos+1 + (num_cells.y) * (zpos+1))]);
          grid.push_back(this->voxelGrid[xpos + (num_cells.x) * (ypos+1 + (num_cells.y) * (zpos+1))]);

          vector<vertex> positions = vector<vertex>();
          positions.push_back(vertex(Vector3D(xpos,ypos,zpos), -gradientNormal(Vector3D(xpos,ypos,zpos))));
          positions.push_back(vertex(Vector3D(xpos+1,ypos,zpos), -gradientNormal(Vector3D(xpos+1,ypos,zpos))));
          positions.push_back(vertex(Vector3D(xpos+1,ypos,zpos+1), -gradientNormal(Vector3D(xpos+1,ypos,zpos+1))));
          positions.push_back(vertex(Vector3D(xpos,ypos,zpos+1), -gradientNormal(Vector3D(xpos,ypos,zpos+1))));
          positions.push_back(vertex(Vector3D(xpos,ypos+1,zpos), -gradientNormal(Vector3D(xpos,ypos+1,zpos))));
          positions.push_back(vertex(Vector3D(xpos+1,ypos+1,zpos), -gradientNormal(Vector3D(xpos+1,ypos+1,zpos))));
          positions.push_back(vertex(Vector3D(xpos+1,ypos+1,zpos+1), -gradientNormal(Vector3D(xpos+1,ypos+1,zpos+1))));
          positions.push_back(vertex(Vector3D(xpos,ypos+1,zpos+1), -gradientNormal(Vector3D(xpos,ypos+1,zpos+1))));

          vector<vector<vertex>> currTriangles = Polygonise(grid,0.1, positions);
          
          for (vector<vertex> &triangle : currTriangles){
            for (vertex &point : triangle){
               
               point.p.x *= sizeCell.x;
               point.p.y *= sizeCell.y;
               point.p.z *= sizeCell.z;

               point.p += min;
               
               //point.n = gradientNormal(point.p);
            }
          }

          triangles.insert(std::end(triangles), std::begin(currTriangles), std::end(currTriangles));
        }
      }
    }
  }

  // Convert triangles to correct coordinates:
  /*for (vector<vertex> &triangle : triangles){
    for (vertex &point : triangle){
       point.n.normalize();
       point.p.x *= sizeCell.x;
       point.p.y *= sizeCell.y;
       point.p.z *= sizeCell.z;

       point.p += min;
    }
  }*/

  //std::cout << triangles.size() << std::endl;
  /*
  vector<double> grid = vector<double>();
  grid.push_back(1);
  grid.push_back(1);
  grid.push_back(0);
  grid.push_back(0);
  grid.push_back(0);
  grid.push_back(0);
  grid.push_back(0);
  grid.push_back(0);

  vector<Vector3D> positions = vector<Vector3D>();
  positions.push_back(Vector3D(0,0,0));
  positions.push_back(Vector3D(1,0,0));
  positions.push_back(Vector3D(1,0,1));
  positions.push_back(Vector3D(0,0,1));
  positions.push_back(Vector3D(0,1,0));
  positions.push_back(Vector3D(1,1,0));
  positions.push_back(Vector3D(1,1,1));
  positions.push_back(Vector3D(0,1,1));
  std::cout << "HELLO" << std::endl;
  std::cout << Polygonise(grid,0.5, positions).size() << std::endl;

  std::cout << "BYE" << std::endl;*/
}



void Fluid::saveFacesToObjs(std::string fileName){
  ofstream myfile;
  myfile.open ("../mitsuba/input/face" + fileName + ".obj");

  std::unordered_map<std::string,int> mymap;


  int numVertex = 1;
  for (vector<vertex> face : triangles){
    vector<int> indicedVertices = vector<int>();
    for (vertex point : face){
      if (mymap.find(std::to_string(point.p.x) + "/" + std::to_string(point.p.y) + "/" +  std::to_string(point.p.z)) == mymap.end()){
        myfile << "v " + std::to_string(point.p.x) + " " + std::to_string(point.p.y) + " " + std::to_string(point.p.z) + "\n";
        myfile << "vn " + std::to_string(point.n.x) + " " + std::to_string(point.n.y) + " " + std::to_string(point.n.z) + "\n";

        mymap[std::to_string(point.p.x) + "/" + std::to_string(point.p.y) + "/" +  std::to_string(point.p.z)] = numVertex;
        numVertex++;
      }
      indicedVertices.push_back(mymap[std::to_string(point.p.x) + "/" + std::to_string(point.p.y) + "/" +  std::to_string(point.p.z)]);
    }
    myfile << "f " + std::to_string(indicedVertices[2]) + "//"  + std::to_string(indicedVertices[2]) + " " + std::to_string(indicedVertices[1]) + "//"  + std::to_string(indicedVertices[1]) + " " + std::to_string(indicedVertices[0])+ "//"  + std::to_string(indicedVertices[0])  + "\n";
    indicedVertices.clear();
  }


  myfile.close();
}

void Fluid::saveVoxelsToMitsuba(std::string fileName, Vector3D min, Vector3D max, bool orientation){
  ofstream fout;
  fout.open(fileName, ios::binary | ios::out);

  char a[4] = {'V', 'O', 'L', (char) 3};
  fout.write((char*) &a, sizeof(a));

  uint32_t encodingId = 1;
  fout.write((char*)&encodingId,sizeof(encodingId));

  uint32_t X = num_cells.x;
  fout.write((char*)&X,sizeof(X));
  uint32_t Y = num_cells.y;
  fout.write((char*)&Y,sizeof(Y));
  uint32_t Z = num_cells.z;
  fout.write((char*)&Z,sizeof(Z));

  if (!orientation){
    uint32_t numChannel = 1;
    fout.write((char*)&numChannel,sizeof(numChannel));
  } else {
    uint32_t numChannel = 3;
    fout.write((char*)&numChannel,sizeof(numChannel));
  }

  float xmin = min.x;
  float ymin = min.y;
  float zmin = min.z;
  float xmax = max.x;
  float ymax = max.y;
  float zmax = max.z;
  fout.write((char*)&xmin,sizeof(xmin));
  fout.write((char*)&ymin,sizeof(ymin));
  fout.write((char*)&zmin,sizeof(zmin));
  fout.write((char*)&xmax,sizeof(xmax));
  fout.write((char*)&ymax,sizeof(ymax));
  fout.write((char*)&zmax,sizeof(zmax));

  if (!orientation){
    for (int xpos = 0; xpos < num_cells.x; ++xpos) {
      for (int ypos = 0; ypos < num_cells.y; ++ypos) {
        for (int zpos = 0; zpos < num_cells.z; ++zpos) {
          float currVal = this->voxelGrid[xpos + num_cells.x * (ypos + num_cells.y * zpos)];
          //std::cout << currVal << std::endl;
          fout.write((char*)&currVal,sizeof(currVal));
        }
      }
    }
  } else{
    for (int xpos = 0; xpos < num_cells.x; ++xpos) {
      for (int ypos = 0; ypos < num_cells.y; ++ypos) {
        for (int zpos = 0; zpos < num_cells.z; ++zpos) {
          Vector3D currVal = this->voxelOrientations[xpos + num_cells.x * (ypos + num_cells.y * zpos)];
          float x = currVal.x;
          float y = currVal.y;
          float z = currVal.z;
          fout.write((char*)&x,sizeof(x));
          fout.write((char*)&y,sizeof(y));
          fout.write((char*)&z,sizeof(z));
        }
      }
    }
  }


  fout.close();
}

void Fluid::simulate(double frames_per_sec, double simulation_steps, FluidParameters *fp,
                     vector<Vector3D> external_accelerations,
                      vector<CollisionObject *> *collision_objects, int step) {
  double delta_t = 1.0f / fps / simulation_steps;
  for (auto &p: particles) {
    p.last_origin = p.origin;
    for (auto ea: external_accelerations){
      p.velocity += delta_t*(ea+p.forces);
    }
    p.x_star = p.origin + delta_t*p.velocity;
  }
  

  std::vector<std::vector<Particle *>>  neighborArray = build_index();
  
  // for surfacing only
  //build_spatial_map();
  build_voxel_grid(step);

  for(int iter=0; iter<solver_iters; iter++) {
    this->update_lambdas(neighborArray);
    this->update_delta_p(neighborArray);
    //apply delta_p and perform collision detection
    for (Particle &p: this->particles) {
      p.x_star += p.delta_p*sf;
    }
    for (Particle &p: this->particles) {
      for (CollisionObject *co : *collision_objects) co->collide_particle(p);
    }
  }
  for (Particle &p: this->particles) {
    p.velocity = (p.x_star-p.origin)/delta_t;
  }

  this->update_omega(neighborArray);
  this->apply_vorticity(neighborArray);
  this->apply_viscosity(neighborArray);

for (Particle &p: this->particles) {
    p.origin = p.x_star;
    //update color
    // p.color = Vector3D( neighborArray[i].size()/10 , 1, 1);
    // p.color = Vector3D( p.density/RHO_O ,0,(int)(p.origin.y <= -0.19));
    p.color = Vector3D(0.05,5*(p.origin.y+0.2)*(p.origin.y+0.2), 1.0);
    // i++;
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
    pm->forces *= 0;
    pm->x_star = pm->origin;
    pm->omega *= 0;
    pm->delta_p *= 0;
    pm->velocity *= 0;
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


void Fluid::update_delta_p(std::vector<std::vector<Particle *>> neighborArray){
  int i = 0;
  for (Particle &p: this->particles) {
    Vector3D p_pred = p.x_star;
    std::vector<Particle *> neighbors = neighborArray[i];
    Vector3D delta = Vector3D(0.0,0.0,0.0);
    float l = p.lambda;
    for (Particle * &pj: neighbors) {
      float scorr = -0.001 * pow(W(p_pred-pj->x_star)/W(Vector3D(0.02, 0.02, 0.02)*R), 4.0);
      Vector3D gradient = del_W(p_pred-pj->x_star);
      delta += (l + pj->lambda+scorr) * gradient;
    }
    p.delta_p = delta / RHO_O;
    i++;
  }
}

void Fluid::apply_vorticity(std::vector<std::vector<Particle *>> neighborArray){
  for (int i = 0; i < particles.size(); i++){
    Particle &p = particles[i];
    Vector3D N;

    p.forces *= 0;
    continue;

    if (p.omega.norm() > 1e-8) {
      std::vector<Particle *> neighbors = neighborArray[i];
      Vector3D pi = p.x_star;
      Vector3D omega_i = p.omega;
      for (Particle * &j: neighbors) {
        N += (j->omega).norm()*del_W(pi-j->x_star); //TODO density might not be included.
      }
      if (N.norm() > 1e-8) N.normalize();
    }
    p.forces = vorticity*CGL::cross(N, p.omega);
  }
}

void Fluid::update_omega(std::vector<std::vector<Particle *>> neighborArray){
  for (int i = 0; i < particles.size(); i++){
    Particle &p = particles[i];
    std::vector<Particle *> neighbors = neighborArray[i];

    Vector3D accum;
    Vector3D pi = p.x_star;
    for (Particle * &j: neighbors) {
      Vector3D vij =  j->velocity - p.velocity;
      accum += CGL::cross(vij,del_W(pi-j->origin));
    }
    p.omega = accum;
  }
}

void Fluid::apply_viscosity(std::vector<std::vector<Particle *>> neighborArray) {
  for (int i = 0; i < particles.size(); i++){
    Particle &p = particles[i];
    std::vector<Particle *> neighbors = neighborArray[i];

    Vector3D accum;
    Vector3D pi = p.x_star;
    for (Particle * &j: neighbors) {
      Vector3D vij =  j->velocity - p.velocity;
      accum += W(pi-j->x_star)*vij;
    }

    p.velocity += viscosity*accum; //TODO viscosity works, but need to use smaller hyperparams than paper
  }

}

// ==================================neighbors=======================================

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

std::vector<std::vector<Particle *>> Fluid::build_index(){
    // Populate cloud
    this->cloud.pts = this->particles;

    if (this->tree != NULL) free(this->tree);
    this->tree =  new kdtree(3, cloud, KDTreeSingleIndexAdaptorParams(3));
    this->tree->buildIndex();

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

      double nMatches = this->tree->radiusSearch(&query_pt[0], R*R, ret_matches, params); //TODO should be R or squared?

      for (auto &pair: ret_matches) {
        to_append.emplace_back(&(this->particles[pair.first]));
      }
      to_return.emplace_back(to_append);
      // cout << to_append.size() << endl;
    }

    return to_return;
}

std::vector<std::vector<Particle *>> Fluid::build_nearest_neighbors_index(int numNeighbors){
    //Finds numNeighbors nearest neighbors for each particle
    // Populate cloud
    this->cloud.pts = this->particles;

    kdtree index(3, cloud, KDTreeSingleIndexAdaptorParams(3));
    index.buildIndex();

    double query_pt[3] = {0.5,0.5,0.5};

    std::vector<std::vector<Particle *>> to_return;
    for  (int k =0; k < particles.size(); k++){
      std::vector<Particle *> to_append;

      query_pt[0] = particles[k].x_star.x;
      query_pt[1] = particles[k].x_star.y;
      query_pt[2] = particles[k].x_star.z;

      std::vector<size_t> ret_index(numNeighbors);
      std::vector<double> dist_sq(numNeighbors);

      int num_results = index.knnSearch(&query_pt[0], numNeighbors, &ret_index[0], &dist_sq[0]);
      for (int i = 0; i < numNeighbors; i++) to_append.emplace_back(&(particles[ret_index[i]]));
      to_return.emplace_back(to_append);
    }
    return to_return;
}

double Fluid::isotropic_kernel(Vector3D pos){
  // Computes the scalar density field interpolated at a point pos
 
  // find particles in radius R

  // sum W(x-x_j)/rho_j


  std::vector<std::pair<size_t,double> > ret_matches;
  SearchParams params;
  params.sorted = false;

  double query_pt[3] = {pos.x,pos.y,pos.z};

  double to_return;

  double nMatches = this->tree->radiusSearch(&query_pt[0], R*R, ret_matches, params);
  for (auto &pair: ret_matches) {
    Particle p = this->particles[pair.first];
    to_return += W(pos-p.x_star)/p.density; //TODO density is from previous time step
  }

  
  return (0.3 - (pos - Vector3D(0.3,0.3,0.3)).norm()); //threshold this at 0.

}

// ==================================state saving=======================================

void Fluid::save_state_to_csv() {
  ofstream fs;
  std::string filename = "state.csv";
  fs.open(filename);

  std::vector<std::vector<Particle *>>  neighborArray = build_index();
  for (int i = 0; i < particles.size(); i++) {
    // write self state
    Particle *p = &particles[i];
    fs << p->origin.x << "," << p->origin.y << "," << p->origin.z << ","
       << p->velocity.x << "," << p->velocity.y << "," << p->velocity.z << std::endl;
    std::cout << neighborArray[i].size() << '\n';
    for (int j = 0; j < 3; j++) {
      p = neighborArray[i][j];
      fs << p->origin.x << "," << p->origin.y << "," << p->origin.z << ","
         << p->velocity.x << "," << p->velocity.y << "," << p->velocity.z << std::endl;
    }
  }

  fs.close();
}
