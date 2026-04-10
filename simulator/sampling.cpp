#include "sampling.h"

void learnSPH::sampling::fluid_box(std::vector<Eigen::Vector3d>& particles, const Eigen::Vector3d& bottom, const Eigen::Vector3d& top, const double sampling_distance)
{
	// Center the sampling
	Eigen::Vector3d b = bottom + 0.5*sampling_distance*Eigen::Vector3d::Ones();
	Eigen::Vector3d t = top - 0.5*sampling_distance*Eigen::Vector3d::Ones();
	for (int i = 0; i < 3; i++) {
		const double l = t[i] - b[i];
		const int n = (int)(l/sampling_distance) + 1;
		const double sampled = n*sampling_distance;
		const double remainder = l - sampled;
		b[i] += 0.5 * remainder;
		t[i] -= 0.5 * remainder;
	}

	// Sample
	const double eps = 1e-10;
	for (double x = b[0]; x < t[0] + eps; x += sampling_distance) {
		for (double y = b[1]; y < t[1] + eps; y += sampling_distance) {
			for (double z = b[2]; z < t[2] + eps; z += sampling_distance) {
				particles.push_back({x, y, z});
			}
		}
	}
}

void learnSPH::sampling::triangle_mesh(std::vector<Eigen::Vector3d>& particles, const std::vector<Eigen::Vector3d>& vertices, const std::vector<std::array<int, 3>>& triangles, const double sampling_distance)
{
	std::vector<std::array<int, 2>> edges = _find_edges(triangles);

	// Sample triangles
	for (const std::array<int, 3> &triangle : triangles) {
		const Eigen::Vector3d& a = vertices[triangle[0]];
		const Eigen::Vector3d& b = vertices[triangle[1]];
		const Eigen::Vector3d& c = vertices[triangle[2]];
		_sample_triangle(particles, a, b, c, sampling_distance);
	}
}

std::vector<std::array<int, 2>> learnSPH::sampling::_find_edges(const std::vector<std::array<int, 3>>& triangles)
{
	// Find the number of nodes
	int n_nodes = 0;
	for (const std::array<int, 3>& triangle : triangles) {
		n_nodes = std::max(n_nodes, triangle[0]);
		n_nodes = std::max(n_nodes, triangle[1]);
		n_nodes = std::max(n_nodes, triangle[2]);
	}

	// Keep track of the already found edges
	std::vector<std::vector<int>> node_edge_map(n_nodes);
	int edge_count = 0;
	for (int elem_i = 0; elem_i < (int)triangles.size(); elem_i++) {
		for (int loc_node_i = 0; loc_node_i < 3; loc_node_i++) {
			for (int loc_node_j = loc_node_i + 1; loc_node_j < 3; loc_node_j++) {
				std::array<int, 2> edge = { std::min(triangles[elem_i][loc_node_i], triangles[elem_i][loc_node_j]),
											std::max(triangles[elem_i][loc_node_i], triangles[elem_i][loc_node_j]) };

				// Check if this edge already exist
				bool found = false;
				for (int glob_node_j : node_edge_map[edge[0]]) {  // Loop through the nodes connected to glob_node_i
					if (glob_node_j == edge[1]) {
						found = true;
						break;
					}
				}

				if (!found) {
					node_edge_map[edge[0]].push_back(edge[1]);
					edge_count++;
				}
			}
		}
	}

	// List unique edges
	std::vector<std::array<int, 2>> edges(edge_count);
	int counter = 0;
	for (int glob_node_i = 0; glob_node_i < n_nodes; glob_node_i++) {
		for (int glob_node_j : node_edge_map[glob_node_i]) {
			edges[counter][0] = glob_node_i;
			edges[counter][1] = glob_node_j;
			counter++;
		}
	}

	return edges;
}

void learnSPH::sampling::_sample_triangle(std::vector<Eigen::Vector3d>& particles, const Eigen::Vector3d& a, const Eigen::Vector3d& b, const Eigen::Vector3d& c, const double sampling_distance)
{
	// Check whether the triangle is too small and just place a particle on the centroid
	const Eigen::Vector3d triangle_centroid = (a + b + c) / 3.0;
	const double triangle_radius = std::max((a - triangle_centroid).norm(), std::max((b - triangle_centroid).norm(), (c - triangle_centroid).norm()));
	if (triangle_radius < sampling_distance) {
		particles.push_back(triangle_centroid);
	}

	// Otherwise use a grid to place the particles and discard the ones outside of the extended triangle
	else {
		// If obtuse angle, pick another configuration
		if ((b - a).dot(c - a) < 0.0) {
			_sample_triangle(particles, b, c, a, sampling_distance);
			return;
		}

		const int particles_begin = (int)particles.size();
		const double particle_radius = 0.5*sampling_distance;
		const double particle_diameter = 2.0 * particle_radius;
		const double length_u = (b - a).norm();
		const double length_v = (c - a).norm();
		const Eigen::Vector3d u = (b - a) / length_u;
		const Eigen::Vector3d normal = u.cross(c - a).normalized();
		const Eigen::Vector3d v = normal.cross(u);

		// Enlarged triangle for pointInTriangle
		const Eigen::Vector3d l_a = a + ((a - b).normalized() + (a - c).normalized()) * 2.0 * particle_radius;
		const Eigen::Vector3d l_b = b + ((b - a).normalized() + (b - c).normalized()) * 2.0 * particle_radius;
		const Eigen::Vector3d l_c = c + ((c - a).normalized() + (c - b).normalized()) * 2.0 * particle_radius;

		const double enlarged_sampling_length_u = std::max(length_u, (c - a).dot(u)) + particle_radius;
		const double enlarged_sampling_length_v = (c - a).dot(v) + particle_radius;

		Eigen::Vector3d particles_centroid = Eigen::Vector3d::Zero();
		// Hexagonal grid
		const double hexa_step_u = 0.8 * particle_diameter;
		const int nu = (int)(enlarged_sampling_length_u / hexa_step_u);
		const int nv = (int)(enlarged_sampling_length_v / particle_diameter);
		for (int i = 0; i < nu; i++) {
			const double offset = (i % 2 != 0) ? -particle_radius : 0.0;
			for (int j = 0; j < nv; j++) {
				const Eigen::Vector3d particle = a + u * i * hexa_step_u + v * (j * particle_diameter + offset);
				if (_point_on_triangle(l_a, l_b, l_c, normal, particle)) {
					particles_centroid += particle;
					particles.push_back(particle);
				}
			}
		}
	}
}

bool learnSPH::sampling::_point_on_triangle(const Eigen::Vector3d& a, const Eigen::Vector3d& b, const Eigen::Vector3d& c, const Eigen::Vector3d& n, const Eigen::Vector3d& p)
{
	const Eigen::Vector3d pab = (b - a).cross(p - a);
	if (pab.dot(n) < 0.0) { return false; }

	const Eigen::Vector3d pbc = (c - b).cross(p - b);
	if (pbc.dot(n) < 0.0) { return false; }

	const Eigen::Vector3d pca = (a - c).cross(p - c);
	if (pca.dot(n) < 0.0) { return false; }

	return true;
}
