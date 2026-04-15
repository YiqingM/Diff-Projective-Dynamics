#pragma once
#include <vector>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>
#include <unordered_set>

namespace pd
{
    class ProjectiveDynamics
    {
    public:
        ProjectiveDynamics() = default;
        ~ProjectiveDynamics() = default;

        void initialize(
            std::vector<Eigen::Vector3d>& positions,
            std::vector<double>& mass,          // per-vertex mass (not inverse)
            std::vector<Eigen::Vector3d>& velocities,
            int sim_iterations,
            int solver_iterations,
            double dt,
            double fps,
            double stiffness,
            std::vector<Eigen::Vector3d>& edge_constraints, // (idx0, idx1)
            std::unordered_set<int> fixed_points
        );

        void run(); // main simulation loop

        std::vector<Eigen::Vector3d> getPositions() const { return m_positions; }
        const std::vector<std::vector<Eigen::Vector3d>>& getXHistory() const { return m_x_history; }
        double computeGradient(
            const std::vector<Eigen::Vector3d>& target_positions,
            double& dLdw
        );

    private:
        // --- Core PD pipeline (called each timestep) ---
        void step(const double dt);
        void computeMomentumEstimate();  // s_n = q_n + h*v_n + h^2 * M^{-1} * f_ext
        void localProjection();          // project each constraint independently -> p_i
        void globalSolve();              // solve the global linear system A * q = rhs
        void updateVelocities();         // v_{n+1} = (q_{n+1} - q_n) / h

        // --- Precomputation (called once in initialize) ---
        // Assemble system matrix A = M/h^2 + sum_i w_i * S_i^T * A_i^T * A_i * S_i
        // and compute its Cholesky factorization.
        // For edge constraints with differential coordinates, A_i maps an edge
        // to the difference of its endpoints. The matrix A is constant as long as
        // constraints don't change, so we only factorize once.
        void prefactorize();

        // --- DiffPD backward (Algorithm 2) ---
        // Solve the adjoint system: nabla^2 g(x) * z = (dL/dx)^T
        // using the local-global iterative solver that reuses A's Cholesky factorization.
        // Returns dL/dy for this timestep and outputs z for dL/dw accumulation.
        std::vector<Eigen::Vector3d> backwardStep(
            const std::vector<Eigen::Vector3d>& x,
            const std::vector<Eigen::Vector3d>& dLdx,
            std::vector<Eigen::Vector3d>& z_out
        );

        // --- Particle state ---
        std::vector<Eigen::Vector3d> m_positions;    // current positions q_n
        std::vector<Eigen::Vector3d> m_velocities;   // current velocities v_n
        std::vector<double> m_mass;                   // per-vertex mass

        // --- PD solver variables ---
        std::vector<Eigen::Vector3d> m_s;             // momentum estimate s_n
        std::vector<Eigen::Vector3d> m_q;             // working positions during local/global iteration
        std::vector<Eigen::Vector3d> m_projections;   // local projection results p_i, one per constraint

        // --- Constraints ---
        // Each entry stores (vertex_idx_0, vertex_idx_1)
        std::vector<Eigen::Vector3d> m_edge_constraints;

        // --- Global system matrix and solver ---
        // Size: N x N where N = number of vertices.
        // We solve the same N x N system three times (once for x, y, z)
        // since the matrix is identical for all three coordinates.
        Eigen::SparseMatrix<double> m_systemMatrix;
        Eigen::SimplicialLLT<Eigen::SparseMatrix<double>> m_solver;

        // --- Simulation parameters ---
        int m_sim_iterations;       // total number of timesteps
        int m_solver_iterations;    // number of local/global iterations per timestep
        double m_dt;                // timestep size
        double m_ct;                // accumulated simulation time
        double m_fps;               // output frame rate
        double m_ct_frame;          // time threshold for next frame output
        double m_frame_time;        // time interval between frames (1/fps)
        int m_c_frame;              // current frame counter
        double m_stiffness;         // constraint weight w_i in the energy
        Eigen::Vector3d gravity = Eigen::Vector3d(0, 0, -9.81); // external force (gravity)
        // --- Fixed particles ---
        // Indices of vertices with Dirichlet boundary conditions
        // (positions remain unchanged during simulation)
        std::unordered_set<int> m_fixed_points;

        // --- Forward history for backward pass ---
        std::vector<std::vector<Eigen::Vector3d>> m_x_history; // positions after each timestep
        std::vector<std::vector<Eigen::Vector3d>> m_y_history; // momentum estimates each timestep
    };
}