#pragma once
#include "../pch.h"
#include <assimp/Importer.hpp>      // C++ importer interface
#include <assimp/scene.h>           // Output data structure
#include <assimp/postprocess.h>     // Post processing flags


/*
the mesh holds a bunch of vertex data and index for each. The mesh can have N vertices then data for each and indexing into it

a mesh holds the data for every vertex, then it holds index values for the triangles
then a model could have multiple meshes

*/
class Triangle;

class MeshCache
{
public:
	struct Mesh
	{
		void reserve_space(size_t _size)
		{
			positions.reserve(_size);
			normals.reserve(_size);
			texcoords.reserve(_size);
			tangents.reserve(_size);
			bitangents.reserve(_size);
			indices.reserve(_size);
		}

		//vertex data
		std::vector<glm::vec3> positions;
		std::vector<glm::vec3> normals;
		std::vector<glm::vec2> texcoords;
		std::vector<glm::vec3> tangents;
		std::vector<glm::vec3> bitangents;

		//triangle indices
		std::vector<unsigned int> indices;
	};
	//models can have multiple meshes 
	struct Model
	{
		std::vector<Mesh> meshes;
		std::string mesh_name;
	};

public:
	MeshCache() = default;
	~MeshCache() = default;

	//no copying/moving should be allowed from this class
	// disallow copy and assignment
	MeshCache(MeshCache const&) = delete;
	MeshCache(MeshCache&&) = delete;
	MeshCache& operator=(MeshCache const&) = delete;
	MeshCache& operator=(MeshCache&&) = delete;

	static void LoadMeshFromFile(std::string filename);
	static void Print();
	//static std::vector<Triangle> ConvertToTriangles(std::string name, glm::mat4 M);


	static std::unordered_map<std::string, Model> modelCache;
};


/*ASSIMP NOTES
	https://assimp-docs.readthedocs.io/en/latest/index.html

	 has a logger if you want to use it
	 has a coordinate system you can change, texture coordinate system, and winding order
	 you have nodes with have children, and they just hold references to the meshes, scene has all the meshes
	 *if im working with heirarchy docs have an example of creating nodes with the translation matrices
	 meshes only have 1 material
	 AI_SUCCESS
	*/

//could binary export data so I don't have to load again
class ASSIMPLoader
{
public:
	static bool Load(std::string path, MeshCache::Model& modeldata);

private:
	static void Process_Node(aiNode* node, const aiScene* scene, MeshCache::Model& modeldata);
	static void Process_Mesh(aiMesh* mesh, const aiScene* scene, MeshCache::Mesh& modeldata);
};