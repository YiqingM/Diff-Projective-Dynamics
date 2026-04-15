#include <stdlib.h>
#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <cmath>

#include <Eigen/Dense>

#include "../simulator/io.h"
#include "../simulator/sampling.h"
#include "../simulator/projectiveDynamics.h"

// Helper: create initial cloth particles, constraints, etc.
// (extracted so we can reuse for multiple runs)
struct ClothSetup {
    std::vector<Eigen::Vector3d> particles;
    std::vector<Eigen::Vector3d> velocities;
    std::vector<Eigen::Vector3d> constraints;
    std::vector<double> mass;
    std::unordered_set<int> fixed_points;
    int sim_iterations;
    int solver_iterations;
    double dt;
    double fps;
};

ClothSetup createClothSetup() {
    ClothSetup setup;
    int w = 40;
    int h = 30;
    double s = 0.1;

    // Particle positions
    for (int i = 0; i < h; ++i)
        for (int j = 0; j < w; ++j)
            setup.particles.push_back(Eigen::Vector3d(-j*2*s, 0.0, -i*2*s));

    // Velocities
    setup.velocities.resize(setup.particles.size(), Eigen::Vector3d::Zero());

    // Distance constraints
    for (int i = 0; i < h; ++i)
        for (int j = 0; j < w - 1; ++j)
            setup.constraints.push_back(Eigen::Vector3d(i*w+j, i*w+j+1, 2*s));

    for (int i = 0; i < h - 1; ++i)
        for (int j = 0; j < w; ++j)
            setup.constraints.push_back(Eigen::Vector3d(i*w+j, (i+1)*w+j, 2*s));

    for (int i = 0; i < h - 1; ++i)
        for (int j = 0; j < w - 1; ++j) {
            setup.constraints.push_back(Eigen::Vector3d(i*w+j, (i+1)*w+j+1, sqrt(2)*2*s));
            setup.constraints.push_back(Eigen::Vector3d(i*w+j+1, (i+1)*w+j, sqrt(2)*2*s));
        }

    // Mass
    setup.mass.resize(setup.particles.size(), 0.5);

    // Fixed points
    setup.fixed_points.insert(0);
    setup.fixed_points.insert(w - 1);

    // Simulation parameters
    setup.sim_iterations = 200;   // reduced for faster optimization
    setup.solver_iterations = 10;
    setup.dt = 0.01;
    setup.fps = 30.0;

    return setup;
}

int main()
{
    ClothSetup setup = createClothSetup();

    // ============ Step 1: Generate ground truth ============
    double gt_stiffness = 1000.0;

    pd::ProjectiveDynamics gt_solver;
    gt_solver.initialize(
        setup.particles, 
        setup.mass, 
        setup.velocities,
        setup.sim_iterations, 
        setup.solver_iterations,
        setup.dt, setup.fps, 
        gt_stiffness,
        setup.constraints, 
        setup.fixed_points
    );
    gt_solver.run();

    std::vector<Eigen::Vector3d> target = gt_solver.getPositions();
    std::cout << "Ground truth generated with stiffness = " << gt_stiffness << std::endl;

    // ============ Step 2: Optimization loop ============
    double w_guess = 1500.0;          // wrong initial guess
    double learning_rate = 3.0;      // may need tuning
    int max_opt_iters = 200;

    for (int opt = 0; opt < max_opt_iters; opt++) {
        // Forward simulation with current guess
        pd::ProjectiveDynamics solver;
        solver.initialize(
            setup.particles, 
            setup.mass, 
            setup.velocities,
            setup.sim_iterations, 
            setup.solver_iterations,
            setup.dt, 
            setup.fps, 
            w_guess,
            setup.constraints, 
            setup.fixed_points
        );
        solver.run();

        // Compute loss and gradient
        double dLdw = 0.0;
        double loss = solver.computeGradient(target, dLdw);

        std::cout << "iter=" << opt
                  << "  loss=" << loss
                  << "  stiffness=" << w_guess
                  << "  dLdw=" << dLdw
                  << std::endl;

        // Check convergence
        if (loss < 1e-8) {
            std::cout << "Converged!" << std::endl;
            break;
        }

        // Gradient descent update
        w_guess -= learning_rate * dLdw;

        // Clamp to positive values
        if (w_guess < 1.0) w_guess = 1.0;
    }

    std::cout << "Final stiffness: " << w_guess
              << " (ground truth: " << gt_stiffness << ")" << std::endl;

    return 0;
}