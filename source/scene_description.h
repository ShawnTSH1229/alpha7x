#pragma once
#include <embree4/rtcore.h>
#include <map>
#include <string>
#include <filesystem>
#include <glm\mat4x4.hpp>

#include "bxdf.h"
#include "pbrt/sampling.h"

class CRay
{
public:
	CRay() = default;
	CRay(glm::vec3 origin, glm::vec3 direction) : origin(origin), direction(direction) {}

	glm::vec3 origin;
	glm::vec3 direction;
};

class CInteraction
{
public:
	glm::vec3 position;
	glm::vec3 norm;
	glm::vec3 wo;
};

struct SShapeSample
{
	CInteraction inter;
	float pdf;
};

struct SLightSample
{
public:
	glm::vec3 wi;
	glm::vec3 L;
	float pdf;

	CInteraction iteraction;
};

/*************************************
* Scene Mesh
**************************************/

struct SA7XGeometry
{
	RTCGeometry geometry;
	int material_idx;
};

struct STriangleMesh
{
	std::vector<int> indices;
	std::vector<glm::vec3> points;
	std::vector<glm::vec3> normals;
	std::vector<glm::vec2> uvs;
};

/*************************************
* Scene Mateials
**************************************/

class CMaterial
{
public:
	inline std::shared_ptr<CBxDF> getBxdf() { return bxdf; }
protected:
	std::shared_ptr<CBxDF> bxdf;
};

class CDiffuseMaterial : public CMaterial
{
public:
	CDiffuseMaterial(glm::vec3 reflectance) :reflectance(reflectance)
	{
		bxdf = std::make_shared<CDiffuseBxDF>(reflectance);
	};
private:
	glm::vec3 reflectance;
};

class CSurfaceInterraction :public CInteraction
{
public:
	bool is_hit = false;
	CMaterial* material = nullptr;

	inline CRay spawnRay(glm::vec3 direction) const { return CRay(position + direction * 1e-4f, direction); };
	inline glm::vec3 OffsetRayOrigin()const { return norm * 1e-4f; }
	inline CBSDF getBSDF() const { return CBSDF(norm, material->getBxdf()); }
};

class CDielectricMaterial : public CMaterial
{
public:
	CDielectricMaterial(float eta) :eta(eta)
	{
		bxdf = std::make_shared<CDielectricBxDF>(eta);
	};
private:
	float eta;
};

/*************************************
* Scene Ray Tracing Accelerator
**************************************/

class SShapeSceneEntity;
class CAccelerator
{
public:
	CAccelerator();
	~CAccelerator();

	CSurfaceInterraction intersection(CRay ray);

	bool traceVisibilityRay(CRay ray, float max_t);

	RTCGeometry createRTCGeometry(SShapeSceneEntity* shape_entity, int ID, const std::filesystem::path& file_path, glm::mat4x4 transform = glm::mat4x4());
	void finalizeRtSceneCreate();

private:
	RTCGeometry readPLY(const std::string& file_name, int ID, glm::mat4x4 transform = glm::mat4x4());

	friend class CAlpa7XScene;

	RTCScene rt_scene;
	RTCDevice rt_device;

	std::map<std::string, int> mat_name_idx_map;
	std::vector<CMaterial*> scene_materials;
	std::vector<SA7XGeometry> scene_geometries;
	std::vector<STriangleMesh> lights_triangles;
};

/*************************************
* Mesh Triangle
**************************************/

class CTriangle
{
public:

	inline float area()
	{
		const std::shared_ptr<STriangleMesh> tri_mesh = triangle_mesh;
		const std::vector<int>& indices = tri_mesh->indices;
		const std::vector<glm::vec3> points = tri_mesh->points;
		glm::i32vec3 vtx_indices(indices[tri_index * 3 + 0], indices[tri_index * 3 + 1], indices[tri_index * 3 + 2]);
		glm::vec3 positions[3] = { points[vtx_indices.x],points[vtx_indices.y], points[vtx_indices.z] };
		return 0.5f * glm::length(glm::cross(positions[1] - positions[0], positions[2] - positions[0]));
	}

	inline SShapeSample sample(glm::vec2 u)
	{
		const std::shared_ptr<STriangleMesh> tri_mesh = triangle_mesh;
		const std::vector<int>& indices = tri_mesh->indices;
		const std::vector<glm::vec3> points = tri_mesh->points;
		const std::vector<glm::vec3> normals = tri_mesh->normals;

		glm::i32vec3 vtx_indices(indices[tri_index * 3 + 0], indices[tri_index * 3 + 1], indices[tri_index * 3 + 2]);
		glm::vec3 positions[3] = { points[vtx_indices.x],points[vtx_indices.y], points[vtx_indices.z] };
		glm::vec3 normal[3] = { normals[vtx_indices.x],normals[vtx_indices.y], normals[vtx_indices.z] };
		glm::vec3 barycentric_coords = sampleUniformTriangle(u); //!

		glm::vec3 sampled_pos = positions[0] * barycentric_coords.x + positions[1] * barycentric_coords.y + positions[2] * barycentric_coords.z;
		glm::vec3 sampled_normal = normal[0] * barycentric_coords.x + normal[1] * barycentric_coords.y + normal[2] * barycentric_coords.z;
		sampled_normal = glm::normalize(sampled_normal);

		return SShapeSample{ CInteraction {sampled_pos,sampled_normal},1.0f / area()};
	}

	inline SShapeSample sample(glm::vec3 sample_position, glm::vec2 u)
	{
		SShapeSample shape_sample = sample(u);

		glm::vec3 wi = shape_sample.inter.position - sample_position;
		float squared_length = (wi.x * wi.x + wi.y * wi.y + wi.z * wi.z);
		if (squared_length == 0)
		{
			return SShapeSample{ CInteraction(),0 };
		}

		wi = glm::normalize(wi);

		float distance = glm::distance(shape_sample.inter.position, sample_position);

		float tri_area = area();
		float cos_theta = std::abs(glm::dot(shape_sample.inter.norm, -wi));
		float sample_pdf = (distance * distance) / (tri_area * cos_theta);

		if (std::isinf(sample_pdf))
		{
			return SShapeSample{ CInteraction(),0 };
		}

		return SShapeSample{ CInteraction {shape_sample.inter.position,shape_sample.inter.norm},sample_pdf };
	}

	std::shared_ptr<STriangleMesh> triangle_mesh;
	int tri_index = -1;
};

/*************************************
* Scene Lights
**************************************/

class CLight
{
public:
	CLight(glm::vec3 ipt_l_emit)
		:l_emit(ipt_l_emit) {};

	virtual glm::vec3 sampleLe(const glm::vec2& u1, const glm::vec2& u2, CRay& ray, glm::vec3& normal_light, float& pdf_pos, float& pdf_dir) = 0;
	virtual SLightSample SampleLi(glm::vec3 position, glm::vec3 normal, glm::vec2 u) = 0;
	glm::vec3 l_emit;
};

class CDiffuseAreaLight : public CLight
{
public:
	CDiffuseAreaLight(CTriangle ipt_triangle, glm::vec3 ipt_l_emit)
		:CLight(ipt_l_emit)
		, triangle(ipt_triangle) {};

	SLightSample SampleLi(glm::vec3 position, glm::vec3 normal, glm::vec2 u);
	glm::vec3 sampleLe(const glm::vec2& u1, const glm::vec2& u2, CRay& ray, glm::vec3& normal_light, float& pdf_pos, float& pdf_dir)override;

	inline glm::vec3 L(glm::vec3 normal, glm::vec3 w)
	{
		if (glm::dot(normal, w) < 0)
		{
			return glm::vec3(0, 0, 0);
		}
		return l_emit;
	}

	CTriangle triangle;
};


