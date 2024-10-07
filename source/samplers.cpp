#include "samplers.h"

CLightSampler::CLightSampler(std::vector<std::shared_ptr<CLight>> input_lights)
{
	std::vector<float> pre_sum;
	const int ipt_lt_num = input_lights.size();
	lights.resize(ipt_lt_num);
	light_cdf.resize(ipt_lt_num);
	pre_sum.resize(ipt_lt_num);

	for (int idx = 0; idx < lights.size(); idx++)
	{
		lights[idx] = input_lights[idx];

		float light_sum = (lights[idx]->l_emit.x + lights[idx]->l_emit.y + lights[idx]->l_emit.z);
		pre_sum[idx] = idx == 0 ? light_sum : (pre_sum[idx - 1] + light_sum);
	}

	float cdf_sum = pre_sum[pre_sum.size() - 1];

	for (int idx = 0; idx < lights.size(); idx++)
	{
		light_cdf[idx] = pre_sum[idx] / cdf_sum;
	}
}

void CLightSampler::Sample(float u, std::shared_ptr<CLight>& light, float& pmf)
{
	float sample_value = u;
	int sample_idx = 0;
	for (sample_idx = 0; sample_idx < lights.size(); sample_idx++)
	{
		if (sample_value < light_cdf[sample_idx])
		{
			break;
		}
	}

	float prob_mass_func = 0.0;
	if (sample_idx == 0)
	{
		prob_mass_func = light_cdf[sample_idx];
	}
	else
	{
		prob_mass_func = light_cdf[sample_idx] - light_cdf[sample_idx - 1];
	}

	light = lights[sample_idx];
	pmf = prob_mass_func;
}
