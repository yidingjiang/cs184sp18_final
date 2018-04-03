## CS 184: Computer Graphics and Imaging, Spring 2018

### 3D Position Based Fluid Simulation and Surfacing

Gauthier Dieppedalle (24362705), Utkarsh Singhal (25674685), YiDing Jiang (25551317)

**Summary**

In this project, we plan to implement the **[position based fluids](http://mmacklin.com/pbf_sig_preprint.pdf)** from Nvidia. The particle simulation will include various properties such as viscosity, surface tension and we will also implement surfacing (surface reconstruction) with OpenGL/Mitsuba. If time permits, we want to try machine learning techniques for accelerating the particle simulation.

**Problem Description**

Fluid simulation is an interesting topic that has many important applications in computer graphics. We are interested in producing photo realistic simulation of water and understanding the underlying physics of water; therefore, we think particle based simulation is an ideal candidate. It is important to have physically realistic simulation of water because many engineering appplications require accurate simulations that reflect the actual physical interaction. The main challenges of this task are: 

1. The algorithm is mathematically involved
2. The simulation requires efficient data strucutre for neighbor finding
3. The simulation must be numerially stable
4. Machine learning would be novelly hard if we do use it

We will first try to replicate the paper and then build on top of the paper if we have sufficient time.

**Goal and Deliverables**

1. **What we plan to deliver**

![oa](goal.png)

This is out goal for replicating the paper. The bottom shows the particle simulation and the top shows the surfacing result plus or minus the bunny. We plan to measure the quality of the simulation by visually inspecting the realism of the simulation and also comparing the result to the paper. Then, we plan to measure the performance of the system by measuring the runtime as a function of the number of particles. We also want to explore different parameters for both simlation and surfacing.

We want to explore the tradeoff curve for number of particles vs. runtime vs. render quality. We are also interested in how to better parametrize the system for machine learning.

2. **What we hope to deliver**

We want to investigate the possible use of machine learning in fluid simulation task. In particular, we are interested in accelerating the simulation by estimating the acceleration of each particle directly via either a regression tree or a neural network.

**Schedule**

1. 04/08/2018 
   - Setting up the infrastructures for the project
   - Finish collision detection and nearest neighbor data strucuture
   - Finish basic visualization for debugging
2. 04/15/2018
   - Make progress on surfacing
   - Finish implementing various helper functions $$C$$ for enforcing the imcompressibility constraints and kernel function $$W$$
   - Create helper function for tensile instability and vorticity confinement
3. 04/22/2018
   - Finish the update loop
   - Finish surfacing
   - Debug
4. 04/29/2018
   - Add machine learning for acceleration if all goes well
   - Performance profiling
   - Create rendering and prepare for presentation

**Resources**

[Position Based Fluids](http://mmacklin.com/pbf_sig_preprint.pdf)

[Surfacing](https://www.cc.gatech.edu/~turk/my_papers/sph_surfaces.pdf)

[Machine learning for fluid](https://www.inf.ethz.ch/personal/ladickyl/fluid_sigasia15.pdf)

[Other resources](http://blog.mmacklin.com/position-based-fluids/)

For software, we plan to used OpenGL and Mitsuba. For hardware, we want to use the hive and a laptop with 970M if needed. Tensorflow also has native C++ API.