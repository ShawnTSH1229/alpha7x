#pragma once
#define GLM_ENABLE_EXPERIMENTAL
#include <glm/gtx/component_wise.hpp>
#include "glm-aabb/AABB.hpp"
#include "samplers.h"
#include "cameras.h"
#include "scene_description.h"

class CIntegrator
{
public:
	CIntegrator(CAccelerator* ipt_accelerator) :accelerator(ipt_accelerator) {};

	virtual void render() = 0;
	CSurfaceInterraction intersect(CRay ray)const;
	bool traceVisibilityRay(CRay ray, float max_t) { return accelerator->traceVisibilityRay(ray, max_t); }
private:
	CAccelerator* accelerator;
};

struct SCameraPixelPath
{
	struct SVisiblePoint
	{
		glm::vec3 position;
		glm::vec3 wo;
		CBSDF bsdf;
		glm::vec3 beta;
	};

	glm::vec3 phi;

	float radius = 0.0;
	SVisiblePoint visible_point;
	glm::vec3 l_d;

	glm::vec3 tau = glm::vec3(0, 0, 0);

	int m = 0;
	float n = 0.0;
};

struct SGridVisiblePointsList
{
	SCameraPixelPath* cam_pixel_path;
	SGridVisiblePointsList* next_point = nullptr;
};

class CSPPMGrid
{
public:
	CSPPMGrid(int ipt_image_area, const std::vector<SCameraPixelPath>& camera_paths);
	void collectGridVisiblePointsList(std::vector<SCameraPixelPath>& camera_paths, const glm::u32vec2 image_size);
	void addPhoton(glm::vec3 photon_pos, glm::vec3 photon_ray_dir, glm::vec3 beta);

private:
	bool getGridOffest(glm::vec3 p, const glm::AABB& bounds, const int grid_res[3], glm::ivec3& out_point);
	uint32_t hashVisPoint(glm::ivec3 photon_grid_index, int hash_size);

	glm::AABB grid_bound;
	std::allocator<SGridVisiblePointsList> node_allocator;
	std::vector<SGridVisiblePointsList*> grids;
	int grid_res[3];
	int image_area;
};

class CSPPMIntegrator : public CIntegrator
{
public:
	CSPPMIntegrator(int max_depth, float initial_radius, CPerspectiveCamera* camera, CSampler* sampler, CAccelerator* ipt_accelerator, std::vector<std::shared_ptr<CLight>> lights);
	void render();

private:
	void traceCameraPath(glm::u32vec2 image_size, int iter_idx);
	void generatePhotons(int photons_per_iteration, int iter_idx, CSPPMGrid& sppm_grid);
	void updateCameraVisiblePoints(glm::u32vec2 image_size);

	glm::vec3 SampleLd(const CSurfaceInterraction& sf_interaction, const CBSDF* bsdf, CSampler* sampler);

	std::vector<SCameraPixelPath>camera_paths;

	float initial_radius;
	std::shared_ptr<CLightSampler> light_sampler;
	int max_depth;
	CPerspectiveCamera* camera;
	int iteration_num = 0;
	CSampler* sampler_prototype;
};

class CPathIntegrator : public CIntegrator
{
public:
	CPathIntegrator(int max_depth, CPerspectiveCamera* camera, CSampler* sampler, CAccelerator* ipt_accelerator, std::vector<std::shared_ptr<CLight>> lights);

	void render();
private:

	glm::vec3 evaluatePixelSample(glm::vec2 pixel_pos, CSampler* sampler);
	glm::vec3 Li(CRay ray, CSampler* sampler);

	glm::vec3 SampleLd(const CSurfaceInterraction& sf_interaction, const CBSDF* bsdf, CSampler* sampler);

	int max_depth;
	std::shared_ptr<CLightSampler> light_sampler;
	CPerspectiveCamera* camera;
	CSampler* sampler_prototype;
};

