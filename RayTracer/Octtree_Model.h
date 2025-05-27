#pragma once
#include "../pch.h"
#include "AssetManager.h"
#include "Shapes.h"
#include "../ThirdParty/AABB_triangle_Moller.h"

//octtree for a triangle mesh, so it will be associated with a mesh, contain only triangles
//this octree should exist in world space, really the bounds are in world space, its input is in worldspace ray
class Octtree_Model
{
public:
	struct index_info
	{
		int mesh_id;
		int tri_id;
	};
	struct node
	{
		node() = default;
		node(Bounds3 bound) : bounds(bound) {}

		Bounds3 bounds;
		std::vector<index_info> triangle_info;
		bool leaf = true;
		int parent_id = -1;
		std::vector<int> child_id;
	};
public:
	Octtree_Model(TriModel& _model) : model(_model) {}
	~Octtree_Model() = default;

	//creation
	void CreateOcttree()
	{
		//each node is a bounding aabb, and the leaf holds a group of triangles
		//you just keep adding each triangle it will find its leaf
		//then if you hold more than x amount a leaf, theres the splitting method where you split it up into 8 more and hopefuly the average is 4 or less
		//(ideally 1)
		//start with the precomputedbounds of the model
		octtree.clear();
		octtree.reserve(10000);
		node root;
		root.bounds = model.Bounds();
		root.leaf = true;
		root.parent_id = 0;
		//std::cout << root.bounds.pmin.x << " " << root.bounds.pmin.y << " " << root.bounds.pmin.z << std::endl;
		//std::cout << root.bounds.pmax.x << " " << root.bounds.pmax.y << " " << root.bounds.pmax.z << std::endl;
		//root.bounds.Transform(glm::scale(glm::mat4(1.0f), glm::vec3(20, 20, 20)));
		std::cout << root.bounds.pmax.x << " " << root.bounds.pmax.y << " " << root.bounds.pmax.z << std::endl;
		std::cout << root.bounds.pmin.x << " " << root.bounds.pmin.y << " " << root.bounds.pmin.z << std::endl;
		octtree.push_back(root);
		std::cout<<octtree[0].bounds.pmax.x<<std::endl;
		
		for (int mesh_id = 0; mesh_id < model.triangles.size(); mesh_id++)
		{
			const std::vector<Triangle>& tris = model.triangles[mesh_id];
			for (int tri_id = 0; tri_id < tris.size(); tri_id++)
			{
				AddTriangle({ mesh_id, tri_id });
				//std::cout << octtree.size()<< " triangle: " << tri_id << std::endl;
			}
		}
	}

	//traversal
	std::optional<LocalSurfaceInfo> Traverse(Ray& ray)
	{
		//do a breath search through the whole tree, but if ray isnt in a parent you can skip all of that branch, keep going
		//if you hit a leaf then check all the triangles in that, and you can change the tmax to early out other ones farther away
		//ultimately return the triangle that hit and is the closest

		index_info info;
		Triangle::TriangleIntersect tri_isect;
		bool found = false;
		float tMax = std::numeric_limits<float>::max();
		//float tMin = std::numeric_limits<float>::max();
		std::queue<int> queue;
		queue.push(0);
		while (!queue.empty())
		{
			int current_id = queue.front();
			queue.pop();

			if (octtree[current_id].bounds.IntersectP(ray, tMax)) {
				if (octtree[current_id].leaf) {
					std::vector<index_info>& leaf_triangle_info = octtree[current_id].triangle_info;
					for (int t = 0; t < leaf_triangle_info.size(); t++)
					{
						Triangle& tri = model.triangles[leaf_triangle_info[t].mesh_id][leaf_triangle_info[t].tri_id];
						if (model.enable_cull_back_face && !model.back_facing.empty() &&
							model.back_facing[leaf_triangle_info[t].mesh_id][leaf_triangle_info[t].tri_id])
						{
							//std::cout << "culled out\n";
							continue;
						}

						std::optional<Triangle::TriangleIntersect> isect = tri.BasicIntersect(ray, tMax);
						if (isect.has_value())
						{
							//tMax = std::min(tMax, isect.value().t);
							if (isect.value().t < tMax)
							{
								found = true;
								tMax = isect.value().t;
								tri_isect = isect.value();
								info = leaf_triangle_info[t];
							}
						}
					}
				}
				else {
					for (int i = 0; i < 8; i++)
					{
						queue.push(octtree[current_id].child_id[i]);
					}
				}
			}
		}
		if (!found) {
			//std::cout << "traverse oct tree: no hits\n";
			return {};
		}
		
		std::optional<LocalSurfaceInfo> local_info = model.triangles[info.mesh_id][info.tri_id].CalculateLocalSurface(tri_isect);

		return { local_info };
	}

	int getTreeSize()
	{
		return octtree.size();
	}

	void PrintInfo()
	{
		int leaf_count = 0;
		float avg_triangle_leaf = 0;
		int empty_leafs = 0;
		int sum = 0;
		int max_triangle_leaf = std::numeric_limits<float>::min();
		int real_nodes = 0;
		std::queue<int> queue;
		queue.push(0);
		while (!queue.empty())
		{
			int current_id = queue.front();
			queue.pop();
			real_nodes++;
			//std::cout << octtree[current_id].triangle_info.size() << std::endl;
			if (!octtree[current_id].leaf)
			{
				for (int i = 0; i < 8; i++)
				{
					queue.push(octtree[current_id].child_id[i]);
				}
			}
			else
			{
				sum += octtree[current_id].triangle_info.size();
				max_triangle_leaf = std::max(max_triangle_leaf, (int)octtree[current_id].triangle_info.size());
				leaf_count++;
				if (octtree[current_id].triangle_info.empty())
					empty_leafs++;
			}
		}
		avg_triangle_leaf = sum / (float)leaf_count;

		//int depth;
		//int last_index = octtree.size() - 1;
		//if (last_index == 0) depth = 1;
		//else if (last_index == 1) depth = 2;
		//else depth = (std::log(last_index - 1) / std::log(8)) + 2;

		std::cout << "node count: " << octtree.size()<<" real node count: "<<real_nodes << " leafs: " << leaf_count << " empty leafs: " << empty_leafs
			<< " avg_leaftriangle: " << avg_triangle_leaf << " maxleaf: " << max_triangle_leaf << std::endl;
	}

	Octtree_Model::node GetNode(int i) { return octtree[i]; }
private:
	void AddTriangle(index_info info)
	{
		//find what leaf it belongs too
		//then place it in there, or split and do that splitting
		//aabb triangle intersection to see if its in the bounds or not
		//can be in multiple ones

		//go through every one starting from root in order and quit if parent isn't good
		MeshCache::Mesh& mesh = modelcache[model.getModelName()].meshes[info.mesh_id];
		glm::vec3 p0 = mesh.positions[mesh.indices[3 * info.tri_id]];
		glm::vec3 p1 = mesh.positions[mesh.indices[3 * info.tri_id + 1]];
		glm::vec3 p2 = mesh.positions[mesh.indices[3 * info.tri_id + 2]];
		if (!model.mesh_has_precomputed_worldposition)
		{
			p0 = model.ObjectToRender * glm::vec4(p0.x, p0.y, p0.z, 1);
			p1 = model.ObjectToRender * glm::vec4(p1.x, p1.y, p1.z, 1);
			p2 = model.ObjectToRender * glm::vec4(p2.x, p2.y, p2.z, 1);
		}

		std::unordered_map<int, bool> parent_invalid_list;
		//what boxes is this triangle inside?
		int current_id = 0;
		int max_id = octtree.size() - 1; //octtree might split but we shouldn't check those ones
		bool found = false;
		std::queue<int> queue;
		//if you hit a parent, push all its children on a stack we know their indices are valid. if leaf don't push anything on
		queue.push(0);
		while (!queue.empty()) {
			current_id = queue.front();
			queue.pop();
			//std::cout << "current id: "<<current_id << std::endl;
			if (tri_boundsIntersection(p0, p1, p2, octtree[current_id].bounds)) {
				//std::cout << "bounds hit " << current_id << " triangles held: " << octtree[current_id].triangle_info.size() << "\n";
				//if (!octtree[current_id].leaf)
					//std::cout << octtree[Child(current_id, 1)].triangle_info.size() << std::endl;
				//if its a parent do nothing, if its a leaf add and see if we need to split
				if (octtree[current_id].leaf) {
					//std::cout << "is a leaf " << current_id << " size: " << octtree[current_id].triangle_info.size() << "\n";
					octtree[current_id].triangle_info.push_back(info);
					found = true;
					if (octtree[current_id].triangle_info.size() >= TRIANGLE_CAPACITY) {
						//std::cout << "splitting\n";
						Split(current_id);
					}
				}
				else {
					//push children in
					//std::cout << "pushing children\n";
					for (int i = 0; i < 8; i++) {
						//std::cout << Child(current_id, i) << std::endl;
						queue.push(octtree[current_id].child_id[i]);
					}
				}
			}
		}

		/*
		while (current_id <= max_id)
		{
			//std::cout << current_id << std::endl;
			if (!parent_invalid_list.contains(Parent(current_id)))
			{
				std::cout << "foudn parent: "<< current_id<<"\n";
				if (tri_boundsIntersection(p0, p1, p2, octtree[current_id].bounds))
				{
					std::cout << "bounds intersection "<< current_id<<" "<< octtree[current_id].triangle_info.size()<< "\n";
					if(!octtree[current_id].leaf)
						std::cout<<octtree[Child(current_id, 1)].triangle_info.size() << std::endl;
					//if its a parent do nothing, if its a leaf add and see if we need to split
					if (octtree[current_id].leaf)
					{
						std::cout << "leaf " << current_id <<" size: "<< octtree[current_id].triangle_info.size()<< "\n";
						octtree[current_id].triangle_info.push_back(info);
						found = true;
						if (octtree[current_id].triangle_info.size() >= TRIANGLE_CAPACITY)
						{
							std::cout << "splitting\n";
							Split(current_id);
						}
					}
				}
				else
				{
					if (parent_invalid_list.contains(current_id)) std::cout << "oct:addtriangle parent already added1\n";
					parent_invalid_list[current_id] = true;
				}
			}
			else
			{
				if (parent_invalid_list.contains(current_id)) std::cout << "oct:addtriangle parent already added2\n";
				parent_invalid_list[current_id] = true;
			}
			current_id++;
		}
		*/
		if (!found)
			std::cout << "triangle didn't find any octbox\n";
	}

	void Split(int split_id)
	{
		Bounds3 B = octtree[split_id].bounds;
		float padding = 0.01f;
		glm::vec3 half_d = glm::vec3(B.pmax.x - B.pmin.x, B.pmax.y - B.pmin.y, B.pmax.z - B.pmin.z) / 2.0f;
		glm::vec3 C = B.pmin + half_d;
		//create 8 new nodes 
		half_d += glm::vec3(padding, padding, padding);
		Bounds3 bound_topfl(C + glm::vec3(-half_d.x, 0, -half_d.z), C + glm::vec3(0, half_d.y, 0));
		Bounds3 bound_topfr(C + glm::vec3(0, 0, -half_d.z), C + glm::vec3(half_d.x, half_d.y, 0));
		Bounds3 bound_topbl(C + glm::vec3(-half_d.x, 0, 0), C + glm::vec3(0, half_d.y, half_d.z));
		Bounds3 bound_topbr(C + glm::vec3(0, 0, 0), C + glm::vec3(half_d.x, half_d.y, half_d.z));

		Bounds3 bound_botfl(C + glm::vec3(-half_d.x, -half_d.y, -half_d.z), C + glm::vec3(0, 0, 0));
		Bounds3 bound_botfr(C + glm::vec3(0, -half_d.y, -half_d.z), C + glm::vec3(half_d.x, 0, 0));
		Bounds3 bound_botbl(C + glm::vec3(-half_d.x, -half_d.y, 0), C + glm::vec3(0, 0, half_d.z));
		Bounds3 bound_botbr(C + glm::vec3(0, -half_d.y, 0), C + glm::vec3(half_d.x, 0, half_d.z));

		node node_topfl(bound_topfl); node node_topfr(bound_topfr); node node_topbl(bound_topbl); node node_topbr(bound_topbr);
		node node_botfl(bound_botfl); node node_botfr(bound_botfr); node node_botbl(bound_botbl); node node_botbr(bound_botbr);
		std::vector<node> nodes = { node_topfl, node_topfr, node_topbl, node_topbr,
									node_botfl, node_botfr, node_botbl, node_botbr };

		//now put each of the triangles into their new boxes (can be multiple)
		for (int i = 0; i < octtree[split_id].triangle_info.size(); i++)
		{
			index_info current_info = octtree[split_id].triangle_info[i];
			MeshCache::Mesh& mesh = modelcache[model.getModelName()].meshes[current_info.mesh_id];
			glm::vec3 p0 = mesh.positions[mesh.indices[3 * current_info.tri_id]];
			glm::vec3 p1 = mesh.positions[mesh.indices[3 * current_info.tri_id + 1]];
			glm::vec3 p2 = mesh.positions[mesh.indices[3 * current_info.tri_id + 2]];
			if (!model.mesh_has_precomputed_worldposition)
			{
				p0 = model.ObjectToRender * glm::vec4(p0.x, p0.y, p0.z, 1);
				p1 = model.ObjectToRender * glm::vec4(p1.x, p1.y, p1.z, 1);
				p2 = model.ObjectToRender * glm::vec4(p2.x, p2.y, p2.z, 1);
			}

			bool found_one = false;
			for (int n = 0; n < nodes.size(); n++)
			{
				if (tri_boundsIntersection(p0, p1, p2, nodes[n].bounds))
				{
					found_one = true;
					nodes[n].triangle_info.push_back(current_info);
				}
			}
			if (!found_one)
			{
				std::cout << "error triangle didnt find new aabb\n";
			}
		}
		//see if should abort the split, like if all of them just go into another bin theres no point in splitting
		int triangle_count = octtree[split_id].triangle_info.size();
		bool leave = false;
		for (int n = 0; n < 8; n++)
		{
			if (nodes[n].triangle_info.size() == triangle_count)
				leave = true;
		}
		if (leave)
			return;
		

		//octtree.reserve(octtree.size() + 8);
		octtree[split_id].child_id.clear();
		octtree[split_id].child_id.reserve(8);

		for (int n = 0; n < nodes.size(); n++)
		{
			nodes[n].parent_id = split_id;
			//std::cout << "node: " << n << " size: " << nodes[n].triangle_info.size() << std::endl;
			octtree.push_back(nodes[n]);
			octtree[split_id].child_id.push_back(octtree.size() - 1);
		}

		//remove triangles from parent
		octtree[split_id].triangle_info.clear();
		octtree[split_id].leaf = false;
	}

	//tri/bounds intersection
	bool tri_boundsIntersection(const glm::vec3& p0, const glm::vec3& p1, const glm::vec3& p2, const Bounds3& bounds)
	{
		glm::vec3 half_d = glm::vec3(bounds.pmax.x - bounds.pmin.x, bounds.pmax.y - bounds.pmin.y, bounds.pmax.z - bounds.pmin.z) / 2.0f;
		glm::vec3 C = bounds.pmin + half_d;

		return Moller::triBoxOverlap(C, half_d, { p0, p1, p2 });

		/*
		//std::cout << p0.x << " " << p0.y << " " << p0.z << std::endl;
		//std::cout << p1.x << " " << p1.y << " " << p1.z << std::endl;
		//std::cout << p2.x << " " << p2.y << " " << p2.z << std::endl;

		if (bounds.pmin.x <= p0.x && p0.x <= bounds.pmax.x && bounds.pmin.y <= p0.y && p0.y <= bounds.pmax.y && bounds.pmin.z <= p0.z && p0.z <= bounds.pmax.z)
			return true;
		
		if (bounds.pmin.x <= p1.x && p1.x <= bounds.pmax.x && bounds.pmin.y <= p1.y && p1.y <= bounds.pmax.y && bounds.pmin.z <= p1.z && p1.z <= bounds.pmax.z)
			return true;

		if (bounds.pmin.x <= p2.x && p2.x <= bounds.pmax.x && bounds.pmin.y <= p2.y && p2.y <= bounds.pmax.y && bounds.pmin.z <= p2.z && p2.z <= bounds.pmax.z)
			return true;
	
		return false;
		*/
	}

public:
	//each level has full amount 8^(level-1). 1,8,64,512,4096,32768
	static const int TRIANGLE_CAPACITY = 40;

	//k starts at (1,8)
	/*
	inline int Child(int parent, int k)
	{
		return parent*8 + k;
	}

	inline int Parent(int i)
	{
		if (i == 0) return 0;
		return std::floor((i - 1) / 8.0f);
	}
	*/
/*
	inline bool isLeaf(int i)
	{
		if (Child(i, 1) > octtree.size() - 1)
			return true;
		return false;
	}
	*/
	//the tree is stored in one big array index is simple
	//kth child of index i (kCi) = i*8 + k
	//parent of index i (Pi) = floor(i-1 / 8), if i=0 then NA
private:
	std::vector<node> octtree;
	TriModel& model;

	std::unordered_map<std::string, MeshCache::Model>& modelcache = MeshCache::modelCache;
};



