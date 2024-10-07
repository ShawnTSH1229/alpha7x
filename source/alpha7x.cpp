#include <embree4/rtcore.h>
#include <stdio.h>
#include <math.h>
#include <limits>
#include <stdio.h>

#if defined(_WIN32)
#  include <conio.h>
#  include <windows.h>
#endif

#include "utils/cxxopts.hpp"
#include "pbrt_parser/parser.h"
#include "scene_builder.h"
#include "integrators.h"

void renderScene(CAlpa7XScene& a7x_scene, std::string& output_image_scene_path)
{
	std::vector<std::shared_ptr<CLight>> lights;
	CPerspectiveCamera* camera = a7x_scene.getCamera();
	camera->output_image_scene_path = output_image_scene_path;

	CSampler* sampler = a7x_scene.getSampler();
	CAccelerator* accel = a7x_scene.createAccelerator(lights);
	std::unique_ptr<CIntegrator> integrator = a7x_scene.createIntegrator(camera, sampler, accel, lights);
	integrator->render();
}


int main(int argc, char* argv[])
{
	cxxopts::Options opts("Alpha7XRender", "Tiny Offline Renderer");
	opts.allow_unrecognised_options();
	opts.add_options()
		("i,input_pbrt_scene", "Input pbrt scene path", cxxopts::value<std::string>())
		("o,output_image", "output image path", cxxopts::value<std::string>())
		("h,help", "Print help message. command example -i X:/path to you project/Alpha7XRenderer/resource/water-caustic/scene-v4.pbrt -o X:/path to you project/Alpha7XRenderer/build/");

	auto opt_result = opts.parse(argc, argv);
	if (opt_result.count("h") || argc < 2 || opt_result.count("i") == 0)
	{
		printf(opts.help().c_str());
		printf("\n");
		exit(-1);
	}

	std::string input_pbrt_scene_path = opt_result["i"].as<std::string>();
	std::string output_image_scene_path = opt_result["o"].as<std::string>();
	
	CAlpa7XScene scene;
	Alpha7XSceneBuilder builder(&scene);
	pbrt::ParseFile(&builder, input_pbrt_scene_path);
	renderScene(scene, output_image_scene_path);

	return 0;
}