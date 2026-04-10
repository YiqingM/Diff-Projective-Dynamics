#include "projectiveDynamics.h"
#include "io.h"
#include <iostream>
#include <omp.h>

namespace pd
{
    void ProjectiveDynamics::initialize(
        std::vector<Eigen::Vector3d>& positions,
        std::vector<double>& mass,
        std::vector<Eigen::Vector3d>& velocities,
        int sim_iterations,
        int solver_iterations,
        double dt,
        double fps,
        double stiffness,
        std::vector<Eigen::Vector3d>& edge_constraints,
        std::unordered_set<int> fixed_points
    ) {
        m_positions = positions;
        m_velocities = velocities;
        m_mass = mass;
        m_edge_constraints = edge_constraints;
        m_fixed_points = fixed_points;

        m_sim_iterations = sim_iterations;
        m_solver_iterations = solver_iterations;
        m_dt = dt;
        m_fps = fps;
        m_stiffness = stiffness;
        m_ct = 0.0;
        m_ct_frame = 0.0;
        m_frame_time = 1.0 / fps;
        m_c_frame = 0;

        // Initialize PD solver variables
        m_s.resize(positions.size(), Eigen::Vector3d::Zero());
        m_q.resize(positions.size(), Eigen::Vector3d::Zero());
        m_projections.resize(edge_constraints.size(), Eigen::Vector3d::Zero());

        // Precompute system matrix and its factorization
        prefactorize();
    }

    void ProjectiveDynamics::prefactorize() {
        int N = (int)m_positions.size();
        int constraint_count = (int)m_edge_constraints.size();

        // Use triplets to build sparse matrix efficiently
        std::vector<Eigen::Triplet<double>> triplets;
        triplets.reserve(N + 4 * constraint_count);

        // Mass term: M / h^2 on the diagonal
        for (int i = 0; i < N; i++) {
            triplets.push_back(Eigen::Triplet<double>(i, i, m_mass[i] / (m_dt * m_dt)));
        }

        // Constraint term: w_i * S_i^T * A_i^T * A_i * S_i for each edge
        // A_i * S_i = [1, -1], so A_i^T * A_i contributes the 2x2 block
        for (int i = 0; i < constraint_count; i++) {
            int idx0 = (int)m_edge_constraints[i][0];
            int idx1 = (int)m_edge_constraints[i][1];

            triplets.push_back(Eigen::Triplet<double>(idx0, idx0,  m_stiffness));
            triplets.push_back(Eigen::Triplet<double>(idx1, idx1,  m_stiffness));
            triplets.push_back(Eigen::Triplet<double>(idx0, idx1, -m_stiffness));
            triplets.push_back(Eigen::Triplet<double>(idx1, idx0, -m_stiffness));
        }

        m_systemMatrix.resize(N, N);
        m_systemMatrix.setFromTriplets(triplets.begin(), triplets.end());

        // Prefactorize: Cholesky decomposition (only done once)
        m_solver.compute(m_systemMatrix);
        if (m_solver.info() != Eigen::Success) {
            std::cerr << "Cholesky factorization failed!" << std::endl;
        }
    }

    void ProjectiveDynamics::run() {
        for (int iter = 0; iter < m_sim_iterations; iter++) {
            m_ct += m_dt;
            step(m_dt);
            if (m_ct >= m_ct_frame)
            {
                // todo: write the particles to file
                learnSPH::write_particles_to_vtk("../res/example_cloth_set" + std::to_string(m_c_frame) + ".vtk", m_positions, m_mass);
                m_ct_frame += m_frame_time;
                m_c_frame += 1;
            }
        }
    }

    void ProjectiveDynamics::step(const double dt) {
        computeMomentumEstimate();
        for (int iter = 0; iter < m_solver_iterations; iter++) {
            localProjection();
            globalSolve();
        }
        updateVelocities();
    }

    void ProjectiveDynamics::computeMomentumEstimate() {
        // s_n = q_n + h*v_n + h^2 * M^{-1} * f_ext
        #pragma omp parallel for
        for (int i = 0; i < (int)m_positions.size(); i++) {
            m_s[i] = m_positions[i] + m_dt * m_velocities[i] + (m_dt * m_dt) * gravity;
        }
        m_q = m_s;
        for (int idx : m_fixed_points) {
            m_s[idx] = m_positions[idx];
            m_q[idx] = m_positions[idx];
        }
    }

    void ProjectiveDynamics::localProjection() {
        #pragma omp parallel for
        for (int i = 0; i < (int)m_edge_constraints.size(); i++) {
            int idx0 = (int)m_edge_constraints[i][0];
            int idx1 = (int)m_edge_constraints[i][1];
            double rest_length = m_edge_constraints[i][2];

            Eigen::Vector3d edge = m_q[idx0] - m_q[idx1];
            double length = edge.norm();

            if (length > 1e-7) {
                m_projections[i] = rest_length * edge / length;
            } else {
                m_projections[i] = Eigen::Vector3d::Zero();
            }
        }
    }

    void ProjectiveDynamics::globalSolve() {
        int N = (int)m_positions.size();

        // We solve A * q_coord = rhs_coord for each coordinate (x, y, z) separately
        // because A is the same N x N matrix for all three
        for (int coord = 0; coord < 3; coord++) {

            // Step 1: Build the right-hand side vector
            Eigen::VectorXd rhs(N);

            // Mass term: M / h^2 * s_n
            for (int i = 0; i < N; i++) {
                rhs(i) = (m_mass[i] / (m_dt * m_dt)) * m_s[i][coord];
            }

            // Constraint term: sum of w_i * S_i^T * A_i^T * B_i * p_i
            // For edge constraint with A_i = B_i = [1, -1]:
            //   idx0 gets +w_i * p_i
            //   idx1 gets -w_i * p_i
            for (int i = 0; i < (int)m_edge_constraints.size(); i++) {
                int idx0 = (int)m_edge_constraints[i][0];
                int idx1 = (int)m_edge_constraints[i][1];

                rhs(idx0) += m_stiffness * m_projections[i][coord];
                rhs(idx1) -= m_stiffness * m_projections[i][coord];
            }

            // Step 2: Solve A * q = rhs using prefactored Cholesky
            Eigen::VectorXd result = m_solver.solve(rhs);

            // Step 3: Write result back to m_q
            for (int i = 0; i < N; i++) {
                m_q[i][coord] = result(i);
            }
        }

        // Enforce fixed points: reset their positions
        for (int idx : m_fixed_points) {
            m_q[idx] = m_positions[idx];
        }
    }

    void ProjectiveDynamics::updateVelocities() {
        #pragma omp parallel for
        for (int i = 0; i < (int)m_positions.size(); i++) {
            m_velocities[i] = (m_q[i] - m_positions[i]) / m_dt;
            m_positions[i] = m_q[i]; // update positions for the next iteration
        }
    }
}