#pragma once
#include <glm/matrix.hpp>
#include <glm/trigonometric.hpp>
#include <glm/vec2.hpp>
#include <glm/vec3.hpp>
#include <glm/common.hpp>
#include <glm/exponential.hpp>
#include <vector>
#include <string>

class CRGBImage
{
public:
	CRGBImage(glm::u32vec2 ipt_img_sz)
		:image_size(ipt_img_sz)
	{
		output_img.resize(image_size.y * image_size.x);
		memset(output_img.data(), 0, sizeof(glm::vec3) * output_img.size());
	};

	inline glm::u32vec2 getImageSize()const { return image_size; }

	inline void addSample(glm::u32vec2 dst_pos, glm::vec3 L)
	{
		assert(dst_pos.x >= 0 && dst_pos.x <= image_size.x);
		assert(dst_pos.y >= 0 && dst_pos.y <= image_size.y);
		int write_idx = dst_pos.x + dst_pos.y * image_size.x;
		output_img[write_idx] += L;
	}

	inline void finalizeRender(float spp)
	{
		out_tga_data.resize(image_size.y * image_size.x);
		for (int idx_x = 0; idx_x < image_size.x; idx_x++)
		{
			for (int idx_y = 0; idx_y < image_size.y; idx_y++)
			{
				int src_idx = idx_x + idx_y * image_size.x;
				glm::vec3 mapped_data = glm::pow(glm::vec3(output_img[src_idx] / spp), glm::vec3(1.0f / 2.2f));
				glm::vec3 norm_data = glm::clamp(mapped_data, glm::vec3(0.0), glm::vec3(1.0)) * 255.0f;

				int dst_idx = idx_x + (image_size.y - idx_y - 1) * image_size.x;
				out_tga_data[dst_idx] = norm_data;
			}
		}
	}

	inline void* getFinalData() { return out_tga_data.data(); }

private:
	glm::u32vec2 image_size;
	std::vector<glm::vec3> output_img;
	std::vector<glm::u8vec3> out_tga_data;
};

class CPerspectiveCamera
{
public:
	CPerspectiveCamera(glm::mat4x4 tran_mat, float ipt_fov, CRGBImage* ipt_rgb_film)
		:camera_to_world(tran_mat)
		, fov(ipt_fov)
		, rgb_image(ipt_rgb_film) 
	{

		// tran_mat = camera from world
		glm::mat4 world_from_camera = glm::inverse(tran_mat);

		vx = glm::normalize(glm::vec3(world_from_camera[0].x, world_from_camera[0].y, world_from_camera[0].z));
		vy = glm::normalize(glm::vec3(world_from_camera[1].x, world_from_camera[1].y, world_from_camera[1].z));
		glm::vec2 img_wh = glm::vec2(rgb_image->getImageSize().x, rgb_image->getImageSize().y);

		const float fov_scale = 1.0f / tanf(glm::radians(fov) * 0.5);
		film_right_bottom = -vx * 0.5f * img_wh.x - vy * 0.5f * img_wh.y + 0.5f * img_wh.y * fov_scale * glm::normalize(glm::vec3(world_from_camera[2].x, world_from_camera[2].y, world_from_camera[2].z));
		camera_pos = glm::vec3(world_from_camera[3].x, world_from_camera[3].y, world_from_camera[3].z);
	}

	inline CRGBImage* getImage() { return rgb_image; }
	
	inline glm::vec3 getCameraPos()const { return camera_pos; }
	inline glm::vec3 getPixelRayDirection(glm::vec2 pixel_pos) { return glm::normalize(glm::vec3(pixel_pos.x * vx + pixel_pos.y * vy + film_right_bottom)); };

	std::string output_image_scene_path;

private:
	// derived from the transform matrix
	glm::vec3 vx; // normalized
	glm::vec3 vy; // normalized
	glm::vec3 film_right_bottom;
	glm::vec3 camera_pos;
	float fov;

	glm::mat4x4 camera_to_world;
	CRGBImage* rgb_image;
};