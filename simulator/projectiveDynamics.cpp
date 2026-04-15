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

        m_x_history.clear();
        m_y_history.clear();

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

        m_y_history.push_back(m_s);

        for (int iter = 0; iter < m_solver_iterations; iter++) {
            localProjection();
            globalSolve();
        }
        updateVelocities();

        m_x_history.push_back(m_positions);
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
    
    // ========================================================================
    // DiffPD backward pass
    // ========================================================================
 
    std::vector<Eigen::Vector3d> ProjectiveDynamics::backwardStep(
        const std::vector<Eigen::Vector3d>& x,
        const std::vector<Eigen::Vector3d>& dLdx,
        std::vector<Eigen::Vector3d>& z_out
    ) {
        // Implements Algorithm 2 from DiffPD (Du et al. 2021).
        //
        // Solves the adjoint system: nabla^2 g(x) * z = (dL/dx)^T
        // using the splitting nabla^2 g(x) = A - deltaA, where:
        //   A     = constant system matrix (prefactorized)
        //   deltaA = sum_c w_c * G_c^T * (dp_c/dx) -- local nonlinear part
        //
        // Iterative scheme (Eq. 23): A * z_{k+1} = deltaA * z_k + (dL/dx)^T
 
        int N = (int)x.size();
        int C = (int)m_edge_constraints.size();
 
        // Initialize z = 0
        std::vector<Eigen::Vector3d> z(N, Eigen::Vector3d::Zero());
 
        for (int iter = 0; iter < m_solver_iterations; iter++) {
 
            // --- Local step: compute deltaA * z ---
            //
            // For each edge constraint c with p_c = d * edge / ||edge||:
            //   dp_c/d(edge) = (d / ||edge||) * (I - e_hat * e_hat^T), e_hat = edge / ||edge||
            //
            // We need: w_c * G_c^T * (dp_c/d(edge)) * (G_c * z), second G_c is needed for d(edege)/d(x)
            // where G_c * z = z[idx0] - z[idx1]
 
            std::vector<Eigen::Vector3d> deltaAz(N, Eigen::Vector3d::Zero());
 
            for (int c = 0; c < C; c++) {
                int idx0 = (int)m_edge_constraints[c][0];
                int idx1 = (int)m_edge_constraints[c][1];
                double rest_length = m_edge_constraints[c][2];
 
                Eigen::Vector3d edge = x[idx0] - x[idx1];
                double length = edge.norm();
 
                if (length > 1e-7) {
                    Eigen::Vector3d edge_hat = edge / length;
 
                    // G_c * z
                    Eigen::Vector3d Gz = z[idx0] - z[idx1];
 
                    // (dp_c/d(edge)) * Gz = (d/||e||) * (Gz - (e_hat . Gz) * e_hat)
                    // This is the Jacobian of normalization applied to Gz
                    double proj = edge_hat.dot(Gz);
                    Eigen::Vector3d dpGz = (rest_length / length) * (Gz - proj * edge_hat);
 
                    // Scatter to global: G_c^T * dpGz
                    deltaAz[idx0] += m_stiffness * dpGz;
                    deltaAz[idx1] -= m_stiffness * dpGz;
                }
            }
 
            // --- Global step: z = A^{-1} * (deltaA * z + dL/dx) ---
            // Reuses the prefactored Cholesky decomposition of A
            // solves A * z = deltaA * z + dL/dx for each coordinate separately.
            for (int coord = 0; coord < 3; coord++) {
                Eigen::VectorXd rhs(N);
                for (int i = 0; i < N; i++) {
                    rhs(i) = deltaAz[i][coord] + dLdx[i][coord];
                }
 
                Eigen::VectorXd result = m_solver.solve(rhs);
 
                for (int i = 0; i < N; i++) {
                    z[i][coord] = result(i);
                }
            }
 
            // Fixed points have zero adjoint (their positions don't change)
            for (int idx : m_fixed_points) {
                z[idx] = Eigen::Vector3d::Zero();
            }
        }
 
        // Output z for dL/dw accumulation
        z_out = z;
 
        // Compute dL/dy = (1/h^2) * M * z  (Eq. 13 from DiffPD)
        std::vector<Eigen::Vector3d> dLdy(N);
        for (int i = 0; i < N; i++) {
            dLdy[i] = (1.0 / (m_dt * m_dt)) * m_mass[i] * z[i];
        }
 
        return dLdy;
    }
 
    double ProjectiveDynamics::computeGradient(
        const std::vector<Eigen::Vector3d>& target_positions,
        double& dLdw
    ) {
        int N = (int)m_positions.size();
        int C = (int)m_edge_constraints.size();

        // Loss is on final positions (which is also q of the last state)
        double loss = 0.0;
        
        // Track dl_dx_next and dl_dv_next separately, x_next is the particle position from last timestep, v_next is the particle velocity.
        std::vector<Eigen::Vector3d> dl_dx_next(N, Eigen::Vector3d::Zero());
        std::vector<Eigen::Vector3d> dl_dv_next(N, Eigen::Vector3d::Zero());
        
        // L: = sum_i ||x_i - target_i||^2
        for (int i = 0; i < N; i++) {
            dl_dx_next[i] = 2.0 * (m_positions[i] - target_positions[i]);
            loss += (m_positions[i] - target_positions[i]).squaredNorm();
        }

        dLdw = 0.0;

        // backprop through time
        for (int t = (int)m_x_history.size() - 1; t >= 0; t--) 
        {
            const auto& x_t = m_x_history[t];

            // ============ Step 3 backward: v_t = (x_t - x_prev)/h ============
            // dl_dv_next: gradient of loss w.r.t. v_t, propagated from future steps.
            // dl_dx_next: gradient of loss w.r.t. x_t, propagated from future steps.
            // v_next depends on two variables:
            //   x_t    with coefficient +1/h  =>  x_t    receives dl_dv_next / h
            //   x_prev with coefficient -1/h  =>  x_prev receives -dl_dv_next / h
            
            // Total gradient for x_t (input to adjoint solve):
            // dl_dx_total = dl_dx_next + dl_dv_next / h
            //
            // Partial gradient for x_prev (will accumulate dl_dy in Step 1 backward):
            //   dl_dx_prev = -dl_dv_next / h
            std::vector<Eigen::Vector3d> dl_dx_total(N);
            std::vector<Eigen::Vector3d> dl_dx_prev(N);
            for (int i = 0; i < N; i++) {
                dl_dx_total[i] = dl_dx_next[i] + dl_dv_next[i] / m_dt;
                dl_dx_prev[i] = -dl_dv_next[i] / m_dt;
            }
            // ============ Step 2 Backward: adjoint solve (Algorithm 2 from DiffPD) ============
            // Input:  dl_dx_total (total gradient w.r.t. x_t, used as dL/dx in Algorithm 2)
            // Output: dl_dy (gradient w.r.t. momentum estimate y_t)
            //         z     (adjoint vector, reused for dL/dw accumulation)
            std::vector<Eigen::Vector3d> z;
            auto dl_dy = backwardStep(x_t, dl_dx_total, z);

            // Fixed points: x = x_prev directly (not from PD solve), so gradient passes through
            for (int idx : m_fixed_points) {
                dl_dy[idx] = dl_dx_total[idx];
            }

            // Accumulate dL/dw at this timestep:
            // Second step of PD here is to solve grad g(x) = 0
            // d(grad g)/d(x) * d(x)/d(w) + d(grad g)/d(w) = 0
            // d(x)/d(w) = -(\grad^2 g)^{-1} * d(grad g)/d(w)
            //   dL/dw = dL/dx * d(x)/d(w) = -dL/dx * (\grad^2 g)^{-1} * d(grad g)/d(w)
            //          = z^T * (d(nabla_g)/dw)
            //          = -sum_c (G_c * z)^T * (G_c * x_t - p_c)
            for (int c = 0; c < C; c++) {
                int idx0 = (int)m_edge_constraints[c][0];
                int idx1 = (int)m_edge_constraints[c][1];
                double rest_length = m_edge_constraints[c][2];

                Eigen::Vector3d edge = x_t[idx0] - x_t[idx1];
                double length = edge.norm();

                Eigen::Vector3d p_c = Eigen::Vector3d::Zero();
                if (length > 1e-7) {
                    p_c = rest_length * edge / length;
                }

                Eigen::Vector3d residual = edge - p_c;
                Eigen::Vector3d Gz = z[idx0] - z[idx1];
                dLdw -= Gz.dot(residual);
            }

            // ============ Step 1 (backward): y = x_prev + h*v_prev + h^2*g ============
            // dy/dx_prev = I   =>  dl_dx_prev += dl_dy
            // dy/dv_prev = hI  =>  dl_dv_prev  = h * dl_dy  (fresh assignment, v_prev only appears in y)
            for (int i = 0; i < N; i++) {
                dl_dx_prev[i] += dl_dy[i];
                dl_dv_next[i] = m_dt * dl_dy[i];  // this becomes dl_dv_next for the next iteration
            }

            // Pass to next iteration: x_prev of this step is x_{t-1}
            dl_dx_next = dl_dx_prev;
        }

        return loss;
    }
}