#include <stdlib.h>     // rand
#include <iostream>
#include <vector>
#include <string>
#include <algorithm>    // std::max

#include <Eigen/Dense>

#include "../simulator/io.h"
#include "../simulator/sampling.h"
#include "../simulator/projectiveDynamics.h"


int main()
{
	int w = 40; // number of points horizontally
    int h = 30; // number of points vertically
    double s = 0.1; // half distance between the particles

	std::vector<Eigen::Vector3d> cloth_particles; // vectors to store particle positions

    // the square cloth has left up position at (0,0,0), length = 20m, and parallel to the z=0 plane
    for (int i = 0; i < h; ++i)
        for (int j = 0; j < w; ++j)
            cloth_particles.push_back(Eigen::Vector3d(-j*2*s, 0.0, -i*2*s));
    
    std::vector<Eigen::Vector3d> solids_velocities;

    solids_velocities.resize(cloth_particles.size(), Eigen::Vector3d::Zero());

    // initialize all the constraints;
    std::vector<Eigen::Vector3d> distance_constraints;
    for (int i = 0; i < h; ++i)
    {
        for (int j = 0; j < w - 1; ++j)
        {
            distance_constraints.push_back(Eigen::Vector3d(i * w + j, i * w + j + 1, 2 * s)); // right
        }
    }

	for (int i = 0; i < h - 1; ++i)
    {
        for (int j = 0; j < w; ++j)
        {
            distance_constraints.push_back(Eigen::Vector3d(i * w + j, (i + 1) * w + j, 2 * s)); // down
        }
    }

    for (int i = 0; i < h - 1; ++i)
    {
        for (int j = 0; j < w - 1; ++j)
        {
            distance_constraints.push_back(Eigen::Vector3d(i * w + j, (i + 1) * w + j + 1, sqrt(2) * 2 * s)); // down right
            distance_constraints.push_back(Eigen::Vector3d(i * w + j + 1, (i + 1) * w + j, sqrt(2) * 2 * s)); // down left
        }
    }

	std::vector<double> solids_mass;
	solids_mass.resize(cloth_particles.size(), 0.5);

	// initialize the fixed points
    std::unordered_set<int> fixed_points;
    fixed_points.insert(0);
    fixed_points.insert(w-1);

	// initialize the PD solver
	pd::ProjectiveDynamics pd_solver;
	pd_solver.initialize(
		cloth_particles, 
		solids_mass, 
		solids_velocities, 
		1000, 
		10, 
		0.01, 
		30.0, 
		1e3, 
		distance_constraints, 
		fixed_points);

	// run the simulation
	pd_solver.run();
	return 0;
}