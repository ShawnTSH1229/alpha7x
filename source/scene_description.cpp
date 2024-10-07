#include "scene_description.h"
#include "scene_builder.h"
#include "rply/rply.h"


void errorFunction(void* userPtr, enum RTCError error, const char* str)
{
	printf("error %d: %s\n", error, str);
}

CAccelerator::CAccelerator()
{
	rt_device = rtcNewDevice(NULL);
	if (!rt_device)
	{
		printf("error %d: cannot create device\n", rtcGetDeviceError(NULL));
	}
	rtcSetDeviceErrorFunction(rt_device, errorFunction, NULL);

	rt_scene = rtcNewScene(rt_device);
}

CAccelerator::~CAccelerator()
{
	for (auto* p_iter : scene_materials)
	{
		delete p_iter;
	}

	for (auto& geo_iter : scene_geometries)
	{
		rtcReleaseGeometry(geo_iter.geometry);
	}

	rtcReleaseScene(rt_scene);
	rtcReleaseDevice(rt_device);
}

CSurfaceInterraction CAccelerator::intersection(CRay ray)
{
	RTCIntersectArguments args;
	rtcInitIntersectArguments(&args);

	RTCRayHit embree_ray;
	embree_ray.ray.org_x = ray.origin.x;
	embree_ray.ray.org_y = ray.origin.y;
	embree_ray.ray.org_z = ray.origin.z;
	embree_ray.ray.dir_x = ray.direction.x;
	embree_ray.ray.dir_y = ray.direction.y;
	embree_ray.ray.dir_z = ray.direction.z;
	embree_ray.ray.tnear = 0;
	embree_ray.ray.tfar = std::numeric_limits<float>::max();
	embree_ray.ray.time = 0;
	embree_ray.ray.mask = -1;
	embree_ray.hit.u = embree_ray.hit.v = 0;
	embree_ray.hit.geomID = RTC_INVALID_GEOMETRY_ID;
	embree_ray.hit.instID[0] = RTC_INVALID_GEOMETRY_ID;
	embree_ray.hit.primID = RTC_INVALID_GEOMETRY_ID;

	rtcIntersect1(rt_scene, &embree_ray, &args);

	CSurfaceInterraction surface_interaction;

	if (embree_ray.hit.geomID != RTC_INVALID_GEOMETRY_ID)
	{
		const SA7XGeometry& scene_geometry = scene_geometries[embree_ray.hit.geomID];
		glm::vec3 hit_normal = glm::normalize(glm::vec3(embree_ray.hit.Ng_x, embree_ray.hit.Ng_y, embree_ray.hit.Ng_z));
		int mat_idx = scene_geometry.material_idx;

		surface_interaction.is_hit = true;
		surface_interaction.position = ray.origin + embree_ray.ray.tfar * ray.direction;
		surface_interaction.norm = faceForward(ray.direction, hit_normal);
		surface_interaction.material = scene_materials[mat_idx];
		surface_interaction.wo = -ray.direction;
	}

	return surface_interaction;
}

bool CAccelerator::traceVisibilityRay(CRay ray, float max_t)
{
	RTCOccludedArguments sargs;
	rtcInitOccludedArguments(&sargs);

	RTCRay visibility_ray;
	visibility_ray.org_x = ray.origin.x;
	visibility_ray.org_y = ray.origin.y;
	visibility_ray.org_z = ray.origin.z;
	visibility_ray.dir_x = ray.direction.x;
	visibility_ray.dir_y = ray.direction.y;
	visibility_ray.dir_z = ray.direction.z;
	//visibility_ray.tnear = 1e-5;
	visibility_ray.tnear = 0.0f;
	visibility_ray.tfar = max_t - 1e-5;
	visibility_ray.time = 0;
	visibility_ray.mask = -1;
	rtcOccluded1(rt_scene, &visibility_ray, &sargs);
	if (visibility_ray.tfar > 0.0)
	{
		return true;
	}
	return false;
}

void ply_error_call_back(p_ply ply, const char* message)
{
	printf("rply error:%s\n", message);
}

int rply_vertex_callback(p_ply_argument argument) 
{
	float* buffer;
	long index, flags;

	ply_get_argument_user_data(argument, (void**)&buffer, &flags);
	ply_get_argument_element(argument, nullptr, &index);

	int stride = (flags & 0x0F0) >> 4;
	int offset = flags & 0x00F;

	buffer[index * stride + offset] = (float)ply_get_argument_value(argument);

	return 1;
}

int rply_face_callback(p_ply_argument argument) 
{
	std::vector<int>* triangle_indices;
	long flags, length, value_index;

	ply_get_argument_user_data(argument, (void**)&triangle_indices, &flags);
	ply_get_argument_property(argument, nullptr, &length, &value_index);

	if (length != 3 ) { assert(false); }
	if (value_index < 0) { return 1; }
	if (value_index >= 0)
	{
		triangle_indices->push_back((int)ply_get_argument_value(argument));
	}

	return 1;
}

RTCGeometry CAccelerator::readPLY(const std::string& file_name, int ID, glm::mat4x4 transform)
{
	p_ply ply = ply_open(file_name.c_str(), ply_error_call_back, 0, nullptr);
	if (!ply) { assert(false); }
	if (ply_read_header(ply) == 0) { assert(false); }

	p_ply_element element = nullptr;
	size_t vertex_count = 0, face_count = 0;

	while ((element = ply_get_next_element(ply, element)) != nullptr) 
	{
		const char* name;
		long nInstances;

		ply_get_element_info(element, &name, &nInstances);
		if (strcmp(name, "vertex") == 0) { vertex_count = nInstances; }
		else if (strcmp(name, "face") == 0) { face_count = nInstances; }
	}

	if (vertex_count == 0 || face_count == 0) { assert(false); }

	std::vector<glm::vec3> positions;
	std::vector<glm::vec3> normals;
	std::vector<int> triangle_indices;

	positions.resize(vertex_count);
	normals.resize(vertex_count);
	triangle_indices.reserve(face_count * 3);

	int stide = 3;

	long vtx_x_elem_num = ply_set_read_cb(ply, "vertex", "x", rply_vertex_callback, positions.data(), (stide << 4) | 0);
	long vtx_y_elem_num = ply_set_read_cb(ply, "vertex", "y", rply_vertex_callback, positions.data(), (stide << 4) | 1);
	long vtx_z_elem_num = ply_set_read_cb(ply, "vertex", "z", rply_vertex_callback, positions.data(), (stide << 4) | 2);
	if (vtx_x_elem_num == 0 || vtx_y_elem_num == 0 || vtx_z_elem_num == 0) { assert(false); }
	

	long vtx_indices_num = ply_set_read_cb(ply, "face", "vertex_indices", rply_face_callback, &triangle_indices, 0);
	if (vtx_indices_num == 0) { assert(false); }

	if (ply_read(ply) == 0) { assert(false); }
	ply_close(ply);

	RTCGeometry geom = rtcNewGeometry(rt_device, RTC_GEOMETRY_TYPE_TRIANGLE);
	float* geo_vertices = (float*)rtcSetNewGeometryBuffer(geom, RTC_BUFFER_TYPE_VERTEX, 0, RTC_FORMAT_FLOAT3, 3 * sizeof(float), positions.size());
	unsigned* geo_indices = (unsigned*)rtcSetNewGeometryBuffer(geom, RTC_BUFFER_TYPE_INDEX, 0, RTC_FORMAT_UINT3, 3 * sizeof(unsigned), triangle_indices.size() / 3);
	
	for (int vtx_idx = 0; vtx_idx < vertex_count; vtx_idx++)
	{
		glm::vec4 transformed_vtx = glm::vec4(positions[vtx_idx], 1.0) * transform;
		positions[vtx_idx] = glm::vec3(transformed_vtx.x, transformed_vtx.y, transformed_vtx.z);
	}

	glm::mat3x3 normal_mat = glm::transpose(glm::inverse(transform));
	for (int vtx_idx = 0; vtx_idx < vertex_count; vtx_idx++)
	{
		normals[vtx_idx] = normal_mat * normals[vtx_idx];
	}

	memcpy(geo_vertices, positions.data(), positions.size() * sizeof(glm::vec3));
	for (int idx = 0; idx < triangle_indices.size(); idx++)
	{
		geo_indices[idx] = unsigned(triangle_indices[idx]);
	}
	rtcCommitGeometry(geom);
	rtcAttachGeometryByID(rt_scene, geom, ID);
	return geom;
}

RTCGeometry CAccelerator::createRTCGeometry(SShapeSceneEntity* shape_entity, int ID, const std::filesystem::path& file_path, glm::mat4x4 transform)
{
	if (shape_entity->name == "trianglemesh")
	{
		std::vector<int> indices = shape_entity->parameters.GetIntArray("indices");
		std::vector<glm::vec3> positions = shape_entity->parameters.GetPoint3fArray("P");
		std::vector<glm::vec3> normals = shape_entity->parameters.GetPoint3fArray("N");
		std::vector<glm::vec2> uvs = shape_entity->parameters.GetPoint2fArray("uv");

		RTCGeometry geom = rtcNewGeometry(rt_device, RTC_GEOMETRY_TYPE_TRIANGLE);
		float* geo_vertices = (float*)rtcSetNewGeometryBuffer(geom, RTC_BUFFER_TYPE_VERTEX, 0, RTC_FORMAT_FLOAT3, 3 * sizeof(float), positions.size());
		unsigned* geo_indices = (unsigned*)rtcSetNewGeometryBuffer(geom, RTC_BUFFER_TYPE_INDEX, 0, RTC_FORMAT_UINT3, 3 * sizeof(unsigned), indices.size() / 3);
		memcpy(geo_vertices, positions.data(), positions.size() * sizeof(glm::vec3));
		for (int idx = 0; idx < indices.size(); idx++)
		{
			geo_indices[idx] = unsigned(indices[idx]);
		}
		rtcCommitGeometry(geom);
		rtcAttachGeometryByID(rt_scene, geom, ID);
		return geom;
	}
	else if (shape_entity->name == "plymesh")
	{
		std::string file_name = shape_entity->parameters.GetOneString("filename", "");
		std::filesystem::path filepath = file_path / std::filesystem::path(file_name);
		return readPLY(filepath.string(),ID, transform);
	}
	return RTCGeometry();
}

void CAccelerator::finalizeRtSceneCreate()
{
	rtcCommitScene(rt_scene);
}

SLightSample CDiffuseAreaLight::SampleLi(glm::vec3 position, glm::vec3 normal, glm::vec2 u)
{
	SShapeSample shape_sample = triangle.sample(position, u);
	if (shape_sample.pdf == 0 || glm::length(glm::abs(shape_sample.inter.position - position)) == 0)
	{
		assert(false);
		return SLightSample();
	}

	glm::vec3 wi = glm::normalize(shape_sample.inter.position - position);
	glm::vec3 le = L(shape_sample.inter.norm, -wi);

	return SLightSample{ wi,le, shape_sample.pdf, shape_sample.inter.position };
}

glm::vec3 CDiffuseAreaLight::sampleLe(const glm::vec2& u1, const glm::vec2& u2, CRay& ray, glm::vec3& normal_light, float& pdf_pos, float& pdf_dir)
{
	pdf_dir = 0;

	SShapeSample shape_sample = triangle.sample(u1);
	glm::vec3 out_dir = sampleConsineHemisphere(u2);
	pdf_dir = cosineHemispherePDF(out_dir.z);

	CTangentBasis tangent_basis = CTangentBasis::fromZ(shape_sample.inter.norm);

	out_dir = tangent_basis.fromLocal(out_dir);
	glm::vec3 le = L(shape_sample.inter.norm, out_dir);

	ray = CRay(shape_sample.inter.position, out_dir);
	normal_light = shape_sample.inter.norm;
	pdf_pos = shape_sample.pdf;

	return le;
}


