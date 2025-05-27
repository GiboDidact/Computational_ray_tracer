#include "../pch.h"
#include "AssetManager.h"
#include "Shapes.h"

//MESHCACHE
std::unordered_map<std::string, MeshCache::Model> MeshCache::modelCache;

void MeshCache::LoadMeshFromFile(std::string filename)
{
	if (modelCache.count(filename) != 0)
	{
		std::cout<<"Loading mesh ", filename, " twice\n";
		return;
	}

	//load information, can contain multiple meshes
	Model& model_info = modelCache[filename];
	model_info.mesh_name = filename;

	ASSIMPLoader::Load(filename, model_info);
	
	std::cout << "model: " << filename << " " << " number of meshes: " << model_info.meshes.size()
		<< " number of vertex: " << model_info.meshes[0].positions.size()
		<< " number of triangles: " << model_info.meshes[0].indices.size() / 3 << std::endl;
}

void MeshCache::Print()
{
	for (auto model : modelCache)
	{
		std::cout << "model name: "<<model.second.mesh_name << std::endl;
		for (int i = 0; i < model.second.meshes.size(); i++)
		{
			std::cout << "vertex count: " << model.second.meshes[i].positions.size() << std::endl;
		}
	}
}

/*
std::vector<Triangle> MeshCache::ConvertToTriangles(std::string name, glm::mat4 M)
{
	std::vector<Triangle> tris;
	if (modelCache.count(name) != 0)
	{
		tris.reserve(10000);
		Model& model = modelCache[name];

		for (int i = 0; i < model.meshes.size(); i++)
		{
			Mesh& mesh = model.meshes[i];
			for (int j = 0; j < mesh.indices.size(); j += 3)
			{
				std::vector<unsigned int>& indices = mesh.indices;
				Triangle t1("tri", M, mesh.positions[indices[j]], mesh.positions[indices[j + 1]], mesh.positions[indices[j + 2]]);
				tris.push_back(t1);
				//std::cout << "index: " << indices[j] << " "<< indices[j+1]<<" "<< indices[j+2]<<std::endl;
				//std::cout << mesh.positions[indices[j]].x << " " << mesh.positions[indices[j+1]].x << " " << mesh.positions[indices[j+2]].x << " \n";
			}
		}
	}

	return tris;
}
*/

//ASSIMP LOADER
bool ASSIMPLoader::Load(std::string path, MeshCache::Model& modeldata)
{
	Assimp::Importer importer;
	//importer.GetErrorString();
	//importer.ReadFile(path, aiProcess_Triangulate | aiProcess_CalcTangentSpace);

	const aiScene* scene = importer.ReadFile(path, aiProcess_Triangulate | aiProcess_CalcTangentSpace | aiProcess_GenNormals);

	if (!scene) {
		std::cout<<importer.GetErrorString(), "\n";
		return false;
	}
	aiNode* node = scene->mRootNode;
	Process_Node(scene->mRootNode, scene, modeldata);
	//scene holds materials and you can loop through materials
	//scene->mMaterials[0]->GetTexture(aiTextureType::aiTextureType_NORMALS, 0);
	//std::cout << "Loaded	" + path + " Successfully\n";

	return true;
}

void ASSIMPLoader::Process_Node(aiNode* node, const aiScene* scene, MeshCache::Model& modeldata)
{
	modeldata.meshes.resize(node->mNumMeshes);
	for (int x = 0; x < node->mNumMeshes; x++)
	{
		modeldata.meshes[x].reserve_space(10000);

		aiMesh* mesh = scene->mMeshes[node->mMeshes[x]];
		Process_Mesh(mesh, scene, modeldata.meshes[x]);
	}

	for (int i = 0; i < node->mNumChildren; i++) {
		Process_Node(node->mChildren[i], scene, modeldata);
	}
}

void ASSIMPLoader::Process_Mesh(aiMesh* mesh, const aiScene* scene, MeshCache::Mesh& modeldata)
{

	//go through every vertex and get: position, normal, and texcoords	
	//position, normals, texcoords, tangent, bitangent 
	bool positions = true;
	bool normals = true;
	bool texcoords = true;
	bool tangentbi = true;
	for (int i = 0; i < mesh->mNumVertices; i++)
	{
		if (mesh->HasPositions())
		{
			glm::vec3 vertex_pos(mesh->mVertices[i].x, mesh->mVertices[i].y, mesh->mVertices[i].z);
			modeldata.positions.push_back(vertex_pos);
		}
		else
		{
			positions = false;
		}

		if (mesh->HasNormals())
		{
			glm::vec3 normal(mesh->mNormals[i].x, mesh->mNormals[i].y, mesh->mNormals[i].z);
			modeldata.normals.push_back(normal);
		}
		else
		{
			modeldata.normals.push_back(glm::vec3(0, 0, 0));
			normals = false;
		}

		if (mesh->mTextureCoords[0])
		{
			glm::vec2 texcoords(mesh->mTextureCoords[0][i].x, mesh->mTextureCoords[0][i].y);
			modeldata.texcoords.push_back(texcoords);
		}
		else
		{
			modeldata.texcoords.push_back(glm::vec2(0,0));
			texcoords = false;
		}

		if (mesh->HasTangentsAndBitangents())
		{
			glm::vec3 tangents(mesh->mTangents[i].x, mesh->mTangents[i].y, mesh->mTangents[i].z);
			modeldata.tangents.push_back(tangents);

			glm::vec3 bitangents(mesh->mBitangents[i].x, mesh->mBitangents[i].y, mesh->mBitangents[i].z);
			modeldata.bitangents.push_back(tangents);
		}
		else
		{
			modeldata.tangents.push_back(glm::vec3(0,0,0));
			modeldata.bitangents.push_back(glm::vec3(0, 0, 0));
			tangentbi = false;
		}

	}

	if (!tangentbi)
	{
		std::cout << "mesh doesn't have tangents and bitangents\n";
	}
	if (!normals)
	{
		std::cout << "mesh doesn't have normals\n";
	}
	if (!texcoords)
	{
		std::cout << "mesh doesn't have uv\n";
	}
	if (!positions)
	{
		std::cout << "mesh doesn't have positions\n";
	}

	//go through every face and get indices
	for (int i = 0; i < mesh->mNumFaces; i++)
	{
		aiFace face = mesh->mFaces[i];
		for (int j = 0; j < face.mNumIndices; j++)
		{
			modeldata.indices.push_back(face.mIndices[j]);
		}
	}

}