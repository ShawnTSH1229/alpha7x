#pragma once
#include <glm/vec2.hpp>
#include <iostream>
#include "pbrt/sampling.h"
#include "scene_description.h"

class CLightSampler
{
public:
	CLightSampler(std::vector<std::shared_ptr<CLight>> input_lights);
	void Sample(float u, std::shared_ptr<CLight>& light, float& pmf);
private:
	std::vector<std::shared_ptr<CLight>>lights;
	std::vector<float> light_cdf;
};

class CSampler
{
public:
	CSampler(int ipt_spp) :samplers_per_pixel(ipt_spp) {};

	virtual void initPixelSample(glm::u32vec2 pos, int sample_index, int dim = 0) = 0;
	virtual float get1D() = 0;
	virtual glm::vec2 get2D() = 0;
	virtual glm::vec2 getPixel2D() = 0;

	inline int getSamplersPerPixel()const { return samplers_per_pixel; }
private:
	int samplers_per_pixel;
};


class CSobelSampler : public CSampler
{
public:
	CSobelSampler(int ipt_spp,glm::ivec2 full_resolution) :CSampler(ipt_spp) 
	{
		scale = roundUpPow2(std::max(full_resolution.x, full_resolution.y));
	};

	inline void initPixelSample(glm::u32vec2 pos, int sample_index, int dim = 0)
	{
		pixel = pos;
		dimension = std::max(2, dim);
		sobol_index = sobolIntervalToIndex(log2Int(scale), sample_index, pixel);
	};

	virtual float get1D() // 0 - 1
	{
		if (dimension >= NSobolDimensions)
		{
			dimension = 2;
		}
		return sampleDimension(dimension++);
	};

	virtual glm::vec2 get2D() // 0 - 1
	{ 
		if (dimension + 1 >= NSobolDimensions)
		{
			dimension = 2;
		}
		glm::vec2 u(sampleDimension(dimension), sampleDimension(dimension + 1));
		dimension += 2;
		return u;
	};

	virtual glm::vec2 getPixel2D()
	{
		glm::vec2 u(sobolSample(sobol_index, 0, SNoRandomizer()), sobolSample(sobol_index, 1, SNoRandomizer()));
		for (int dim = 0; dim < 2; ++dim)
		{
			u[dim] = glm::clamp(u[dim] * scale - pixel[dim], 0.0f, float_one_minus_epsilon);
		}
		assert(u.x >= 0 && u.x <= 1 && u.y >= 0 && u.y <= 1);
		return u;
	};

private:

	float sampleDimension(int dimension)const
	{
		uint32_t seed = std::hash<uint32_t>()(dimension);
		return sobolSample(sobol_index, dimension, SFastOwenScrambler(seed));
	}

	float scale;
	glm::u32vec2 pixel;
	int dimension;
	int64_t sobol_index;
};
