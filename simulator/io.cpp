#include "io.h"

#define TINYOBJLOADER_IMPLEMENTATION // define this in only *one* .cc
#include "tiny_obj_loader.h"

void learnSPH::write_particles_to_vtk(std::string path, const std::vector<Eigen::Vector3d>& particles)
{
	write_particles_to_vtk(path, particles, std::vector<double>(), std::vector<Eigen::Vector3d>());
}

void learnSPH::write_empty(std::string path)
{
	std::vector<Eigen::Vector3d> empty;
	write_particles_to_vtk(path, empty);
}

void learnSPH::write_tri_mesh_to_vtk(std::string path, const std::vector<Eigen::Vector3d>& vertices, const std::vector<std::array<int, 3>>& triangles, const std::vector<Eigen::Vector3d>& normals)
{
	// Open the file
	std::ofstream outfile(path, std::ios::binary);
	if (!outfile) {
		std::cout << "learnSPH error in saveTriMeshToVTK: Cannot open the file " << path << std::endl;
		abort();
	}

	// Parameters
	int n_particles = (int)vertices.size();
	int n_triangles = (int)triangles.size();

	// Header
	outfile << "# vtk DataFile Version 4.2\n";
	outfile << "\n";
	outfile << "BINARY\n";
	outfile << "DATASET UNSTRUCTURED_GRID\n";

	// Vertices
	{
		outfile << "POINTS " << n_particles << " double\n";
		std::vector<double> particles_to_write;
		particles_to_write.reserve(3 * n_particles);
		for (const Eigen::Vector3d& vertex : vertices) {
			particles_to_write.push_back(vertex[0]);
			particles_to_write.push_back(vertex[1]);
			particles_to_write.push_back(vertex[2]);
		}
		_swap_bytes_inplace<double>(&particles_to_write[0], (int)particles_to_write.size());
		outfile.write(reinterpret_cast<char*>(&particles_to_write[0]), particles_to_write.size() * sizeof(double));
		outfile << "\n";
	}

	// Connectivity
	{
		outfile << "CELLS " << n_triangles << " " << 4 * n_triangles << "\n";
		std::vector<int> connectivity_to_write;
		connectivity_to_write.reserve(4 * n_triangles);
		for (int tri_i = 0; tri_i < n_triangles; tri_i++) {
			connectivity_to_write.push_back(3);
			connectivity_to_write.push_back(triangles[tri_i][0]);
			connectivity_to_write.push_back(triangles[tri_i][1]);
			connectivity_to_write.push_back(triangles[tri_i][2]);
		}
		_swap_bytes_inplace<int>(&connectivity_to_write[0], (int)connectivity_to_write.size());
		outfile.write(reinterpret_cast<char*>(&connectivity_to_write[0]), connectivity_to_write.size() * sizeof(int));
		outfile << "\n";
	}

	// Cell types
	{
		outfile << "CELL_TYPES " << n_triangles << "\n";
		int cell_type_swapped = 5;
		_swap_bytes_inplace<int>(&cell_type_swapped, 1);
		std::vector<int> cell_type_arr(n_triangles, cell_type_swapped);
		outfile.write(reinterpret_cast<char*>(&cell_type_arr[0]), cell_type_arr.size() * sizeof(int));
		outfile << "\n";
	}

	// Normals
	if (normals.size() > 0) {
		outfile << "POINT_DATA " << n_particles << "\n";
		outfile << "FIELD FieldData 1\n";
		outfile << "vector 3 " << n_particles << " double\n";
		std::vector<double> normals_to_write;
		normals_to_write.reserve(3 * n_particles);
		for (const Eigen::Vector3d& normal : normals) {
			normals_to_write.push_back(normal[0]);
			normals_to_write.push_back(normal[1]);
			normals_to_write.push_back(normal[2]);
		}
		_swap_bytes_inplace<double>(&normals_to_write[0], (int)normals_to_write.size());
		outfile.write(reinterpret_cast<char*>(&normals_to_write[0]), normals_to_write.size() * sizeof(double));
		outfile << "\n";
	}
}

void learnSPH::write_particles_to_vtk(std::string path, const std::vector<Eigen::Vector3d>& particles, const std::vector<double>& particle_scalar_data, const std::vector<Eigen::Vector3d>& particle_vector_data)
{
	// Input checking
	if ((path.size() < 5) || (path.substr(path.size() - 4) != ".vtk")) {
		std::cout << "learnSPH error in saveParticlesToVTK: Filename does not end with '.vtk'" << std::endl;
		abort();
	}

	if ((particle_scalar_data.size() != 0) && (particles.size() != particle_scalar_data.size())) {
		std::cout << "learnSPH error in saveParticlesToVTK: Number of particles not equal to particle scalar data." << std::endl;
		abort(); 
	}
	if ((particle_vector_data.size() != 0) && (particles.size() != particle_vector_data.size())) {
		std::cout << "learnSPH error in saveParticlesToVTK: Number of particles not equal to particle vector data." << std::endl;
		abort(); 
	}

	// Open the file
	std::ofstream outfile(path, std::ios::binary);
	if (!outfile) {
		std::cout << "learnSPH error in saveParticlesToVTK: Cannot open the file " << path << std::endl;
		abort();
	}

	// Parameters
	int n_particles = (int)particles.size();

	if (n_particles == 0) {
		// Header
		outfile << "# vtk DataFile Version 4.2\n";
		outfile << "\n";
		outfile << "ASCII\n";
		outfile << "DATASET UNSTRUCTURED_GRID\n";
		outfile << "POINTS 0 double\n";
		outfile << "CELLS 0 0\n";
		outfile << "CELL_TYPES 0\n";
		outfile.close();
		return;
	}

	// Header
	outfile << "# vtk DataFile Version 4.2\n";
	outfile << "\n";
	outfile << "BINARY\n";
	outfile << "DATASET UNSTRUCTURED_GRID\n";

	// Vertices
	{
		outfile << "POINTS " << n_particles << " double\n";
		std::vector<double> particles_to_write;
		particles_to_write.reserve(3 * n_particles);
		for (const Eigen::Vector3d& vertex : particles) {
			particles_to_write.push_back(vertex[0]);
			particles_to_write.push_back(vertex[1]);
			particles_to_write.push_back(vertex[2]);
		}
		_swap_bytes_inplace<double>(&particles_to_write[0], (int)particles_to_write.size());
		outfile.write(reinterpret_cast<char*>(&particles_to_write[0]), particles_to_write.size() * sizeof(double));
		outfile << "\n";
	}

	// Connectivity
	{
		outfile << "CELLS " << n_particles << " " << 2 * n_particles << "\n";
		std::vector<int> connectivity_to_write;
		connectivity_to_write.reserve(2 * n_particles);
		for (int particle_i = 0; particle_i < (int)particles.size(); particle_i++) {
			connectivity_to_write.push_back(1);
			connectivity_to_write.push_back(particle_i);
		}
		_swap_bytes_inplace<int>(&connectivity_to_write[0], (int)connectivity_to_write.size());
		outfile.write(reinterpret_cast<char*>(&connectivity_to_write[0]), connectivity_to_write.size() * sizeof(int));
		outfile << "\n";
	}

	// Cell types
	{
		outfile << "CELL_TYPES " << n_particles << "\n";
		int cell_type_swapped = 1;
		_swap_bytes_inplace<int>(&cell_type_swapped, 1);
		std::vector<int> cell_type_arr(n_particles, cell_type_swapped);
		outfile.write(reinterpret_cast<char*>(&cell_type_arr[0]), cell_type_arr.size() * sizeof(int));
		outfile << "\n";
	}

	// Point data
	{
		int num_fields = 0;
		if (particle_scalar_data.size() > 0) { num_fields++; }
		if (particle_vector_data.size() > 0) { num_fields++; }

		outfile << "POINT_DATA " << n_particles << "\n";
		outfile << "FIELD FieldData " << std::to_string(num_fields) << "\n";


		if (particle_scalar_data.size() > 0) {
			outfile << "scalar" << " 1 " << n_particles << " double\n";
			std::vector<double> scalar_to_write;
			scalar_to_write.insert(scalar_to_write.end(), particle_scalar_data.begin(), particle_scalar_data.end());
			_swap_bytes_inplace<double>(&scalar_to_write[0], (int)scalar_to_write.size());
			outfile.write(reinterpret_cast<char*>(&scalar_to_write[0]), scalar_to_write.size() * sizeof(double));
			outfile << "\n";
		}

		if (particle_vector_data.size() > 0) {
			outfile << "vector" << " 3 " << n_particles << " double\n";
			std::vector<double> vector_to_write;
			vector_to_write.reserve(3 * n_particles);
			for (const Eigen::Vector3d& vector : particle_vector_data) {
				vector_to_write.push_back(vector[0]);
				vector_to_write.push_back(vector[1]);
				vector_to_write.push_back(vector[2]);
			}
			_swap_bytes_inplace<double>(&vector_to_write[0], (int)vector_to_write.size());
			outfile.write(reinterpret_cast<char*>(&vector_to_write[0]), vector_to_write.size() * sizeof(double));
			outfile << "\n";
		}
	}

	outfile.close();
}

void learnSPH::write_particles_to_vtk(std::string path, const std::vector<Eigen::Vector3d>& particles, const std::vector<double>& particle_scalar_data)
{
	write_particles_to_vtk(path, particles, particle_scalar_data, std::vector<Eigen::Vector3d>());
}

std::vector<learnSPH::TriMesh> learnSPH::read_tri_meshes_from_obj(std::string filename)
{
	std::vector<TriMesh> tri_meshes;

	// Initialise tinyobjloader objects and read the file
	tinyobj::attrib_t attrib;
	std::vector<tinyobj::shape_t> shapes;
	std::vector<tinyobj::material_t> materials;
	std::string warn;
	std::string err;
	bool ret = tinyobj::LoadObj(&attrib, &shapes, &materials, &warn, &err, filename.c_str());

	// Input checking
	if (!warn.empty()) {
		std::cout << warn << std::endl;
	}

	if (!err.empty()) {
		std::cerr << err << std::endl;
	}

	if (!ret) {
		exit(1);
	}

	// Global information
	const int total_n_vertices = (int)attrib.vertices.size() / 3;

	// Write the geometric information into individual triangular meshes
	// Loop over meshes
	for (int shape_i = 0; shape_i < shapes.size(); shape_i++) {

		// Initialize individual triangular mesh
		TriMesh tri_mesh;
		std::vector<bool> global_nodes_present(total_n_vertices, false);

		// Loop over triangles
		int index_offset = 0;
		for (int tri_i = 0; tri_i < shapes[shape_i].mesh.num_face_vertices.size(); tri_i++) {
			if (shapes[shape_i].mesh.num_face_vertices[tri_i] != 3) {
				std::cout << "learnSPH error: readTriMeshesFromObj can only read triangle meshes." << std::endl;
			}

			// Gather triangle global indices
			std::array<int, 3> triangle_global_indices;
			for (int vertex_i = 0; vertex_i < 3; vertex_i++) {
				tinyobj::index_t idx = shapes[shape_i].mesh.indices[(int)(3 * tri_i + vertex_i)];
				const int global_vertex_index = idx.vertex_index;
				triangle_global_indices[vertex_i] = global_vertex_index;
				global_nodes_present[global_vertex_index] = true;
			}
			tri_mesh.triangles.push_back(triangle_global_indices);
		}

		// Reduce global indexes to local indexes
		std::vector<int> global_to_local_vertex_idx(total_n_vertices, -1);
		int local_vertices_count = 0;
		for (int global_vertex_i = 0; global_vertex_i < total_n_vertices; global_vertex_i++) {
			if (global_nodes_present[global_vertex_i]) {
				// Map global -> local
				global_to_local_vertex_idx[global_vertex_i] = local_vertices_count;
				local_vertices_count++;

				// Add vertex to the local mesh vertex vector
				tinyobj::real_t vx = attrib.vertices[(int)(3 * global_vertex_i + 0)];
				tinyobj::real_t vy = attrib.vertices[(int)(3 * global_vertex_i + 1)];
				tinyobj::real_t vz = attrib.vertices[(int)(3 * global_vertex_i + 2)];
				tri_mesh.vertices.push_back({ vx, vy, vz });
			}
		}

		// Change triangle indices
		for (int tri_i = 0; tri_i < tri_mesh.triangles.size(); tri_i++) {
			for (int vertex_i = 0; vertex_i < 3; vertex_i++) {
				tri_mesh.triangles[tri_i][vertex_i] = global_to_local_vertex_idx[tri_mesh.triangles[tri_i][vertex_i]];
			}
		}

		tri_meshes.push_back(tri_mesh);
	}

	return tri_meshes;
}
