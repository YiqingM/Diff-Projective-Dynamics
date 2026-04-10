#pragma once
#include <vector>
#include <array>
#include <fstream>
#include <iostream>

#include <Eigen/Dense>

namespace learnSPH
{
	// .OBJ
	struct TriMesh
	{
		std::vector<Eigen::Vector3d> vertices;
		std::vector<std::array<int, 3>> triangles;
	};
	std::vector<TriMesh> read_tri_meshes_from_obj(std::string filename);

	// VTK
	void write_tri_mesh_to_vtk(std::string path, const std::vector<Eigen::Vector3d>& vertices, const std::vector<std::array<int, 3>>& triangles, const std::vector<Eigen::Vector3d>& normals = std::vector<Eigen::Vector3d>());
	void write_particles_to_vtk(std::string path, const std::vector<Eigen::Vector3d>& particles, const std::vector<double>& particle_scalar_data, const std::vector<Eigen::Vector3d>& particle_vector_data);
	void write_particles_to_vtk(std::string path, const std::vector<Eigen::Vector3d>& particles, const std::vector<double>& particle_scalar_data);
	void write_particles_to_vtk(std::string path, const std::vector<Eigen::Vector3d>& particles);
	void write_empty(std::string path);

	// Aux
	template<class T>
	void _swap_bytes_inplace(T* arr, const int size)
	{
		constexpr unsigned int N = sizeof(T);
		constexpr unsigned int HALF_N = N / 2;
		char* char_arr = reinterpret_cast<char*>(arr);
		for (int e = 0; e < size; e++) {
			for (int w = 0; w < HALF_N; w++) {
				std::swap(char_arr[e * N + w], char_arr[(e + 1) * N - 1 - w]);
			}
		}
	}
}


