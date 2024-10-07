#pragma once
#include <glm/vec3.hpp>
#include <glm/vec2.hpp>
#include <memory>
#include "pbrt/sampling.h"

inline void coordinateSystem(glm::vec3 v1, glm::vec3& v2, glm::vec3& v3)
{
    float sign = v1.z > 0 ? 1.0 : -1.0;
    float a = (-1.0) / (sign + v1.z);
    float b = v1.x * v1.y * a;
    v2 = glm::vec3(1 + sign * a * (v1.x * v1.x), sign * b, -sign * v1.x);
    v3 = glm::vec3(b, sign + a * (v1.y * v1.y), -v1.y);
}

class CTangentBasis
{
public:
    CTangentBasis() :x(glm::vec3(1, 0, 0)), y(glm::vec3(0, 1, 0)), z(glm::vec3(0, 0, 1)) {};
    CTangentBasis(glm::vec3 x, glm::vec3 y, glm::vec3 z) :x(x), y(y), z(z) {};

    static CTangentBasis fromZ(glm::vec3 ipt_z)
    {
        glm::vec3 opt_x, opt_y;
        coordinateSystem(ipt_z, opt_x, opt_y);
        return CTangentBasis(opt_x, opt_y, ipt_z);
    }

    glm::vec3 toLocal(glm::vec3 v)const
    {
        return glm::vec3(glm::dot(v, x), glm::dot(v, y), glm::dot(v, z));
    }

    glm::vec3 fromLocal(glm::vec3 v)const
    {
        return glm::vec3(v.x * x + v.y * y + v.z * z);
    }

    glm::vec3 x, y, z;
};

enum EBxDFFlags
{
	FG_None = 0,
    FG_Reflection = 1 << 0,
    FG_Transmission = 1 << 1,

    FG_Diffuse = 1 << 2,
    FG_Glossy = 1 << 3,
    FG_Specular = 1 << 4,

    FG_DiffuseReflection = FG_Diffuse | FG_Reflection,
    FG_SpecularReflection = FG_Specular | FG_Reflection,
    FG_SpecularTransmission = FG_Specular | FG_Transmission,
};

enum EBxDFReflTransFlags
{
    RF_None = 0,
    RF_Reflection = 1 << 0,
    RF_Transmission = 1 << 1,
    BXDF_RF_All = RF_Reflection | RF_Transmission
};

inline bool isNonSpecular(EBxDFFlags flag) { return flag & (EBxDFFlags::FG_Diffuse | EBxDFFlags::FG_Glossy); }
inline bool isDiffuse(EBxDFFlags flag) { return flag & EBxDFFlags::FG_Diffuse; }
inline bool isGlossy(EBxDFFlags flag) { return flag & EBxDFFlags::FG_Glossy; }

struct SBSDFSample
{
    SBSDFSample() = default;
    SBSDFSample(glm::vec3 f, glm::vec3 wi, float pdf, EBxDFFlags flag, float eta = 1.0)
        :f(f), wi(wi), eta(eta), pdf(pdf), flag(flag) {};

    inline bool isTransmission() { return false; };

    glm::vec3 f;
    glm::vec3 wi;
    float eta;
    float pdf = 0.0;
    EBxDFFlags flag;
};

enum ETransportMode
{
    TM_Radiance,
    TM_Importance,
};

// include BSDF and BTDF
class CBxDF
{
public:
    CBxDF() = default;
    virtual glm::vec3 f(glm::vec3 wo, glm::vec3 wi, ETransportMode transport_mode) = 0;
    virtual float pdf(glm::vec3 wo, glm::vec3 wi, ETransportMode transport_mode, EBxDFReflTransFlags reflect_flag) = 0;
    virtual  std::shared_ptr<SBSDFSample> sample_f(glm::vec3 wo, float u, glm::vec2 u2, ETransportMode transport_mode, EBxDFReflTransFlags reflect_flag) = 0;
    virtual EBxDFFlags flags()const = 0;
private:
};

// lambert
class CDiffuseBxDF : public CBxDF
{
public:
    CDiffuseBxDF(glm::vec3 reflectance)
        :reflectance(reflectance) {}

    inline glm::vec3 f(glm::vec3 wo, glm::vec3 wi, ETransportMode transport_mode)
    {
        if (!sameHemiSphere(wo, wi))
        {
            return glm::vec3(0, 0, 0);
        }
        return reflectance / glm::pi<float>();
    }

    inline float pdf(glm::vec3 wo, glm::vec3 wi, ETransportMode transport_mode, EBxDFReflTransFlags reflect_flag = EBxDFReflTransFlags::RF_Reflection)
    {
        if (!sameHemiSphere(wo, wi) || !(reflect_flag & EBxDFReflTransFlags::RF_Reflection))
        {
            return 0;
        }

        return cosineHemispherePDF(std::abs(wi.z));
    }

    inline std::shared_ptr<SBSDFSample> sample_f(glm::vec3 wo, float u, glm::vec2 u2, ETransportMode transport_mode, EBxDFReflTransFlags reflect_flag)
    {
        if (!(reflect_flag & EBxDFReflTransFlags::RF_Reflection))
        {
            return nullptr;
        }

        glm::vec3 wi = sampleConsineHemisphere(u2);

        if (wo.z < 0)
        {
            wi.z *= -1.0;
        }

        float pdf = cosineHemispherePDF(std::abs(wi.z));
        return std::make_shared<SBSDFSample>(reflectance / glm::pi<float>(), wi, pdf, EBxDFFlags::FG_DiffuseReflection);
    }

    EBxDFFlags flags()const
    {
        return FG_DiffuseReflection;
    }
private:
    glm::vec3 reflectance;
};

class CDielectricBxDF : public CBxDF
{
public:
    CDielectricBxDF(float eta) :eta(eta){}

    inline glm::vec3 f(glm::vec3 wo, glm::vec3 wi, ETransportMode transport_mode) {return glm::vec3(0, 0, 0);}
    inline float pdf(glm::vec3 wo, glm::vec3 wi, ETransportMode transport_mode, EBxDFReflTransFlags reflect_flag) { return 0.0f; }

    std::shared_ptr<SBSDFSample> sample_f(glm::vec3 wo, float u, glm::vec2 u2, ETransportMode transport_mode, EBxDFReflTransFlags reflect_flag);
    
    EBxDFFlags flags()const
    {
        EBxDFFlags bxdf_flag = (eta == 1) ? FG_Transmission : EBxDFFlags(FG_Reflection | FG_Transmission);
        bxdf_flag = EBxDFFlags(bxdf_flag | FG_Specular);
        return bxdf_flag;
    }

private:
    float eta;
};

class CBSDF
{
public:
    CBSDF() :normal(glm::vec3(0, 0, 0)) {};
    CBSDF(glm::vec3 ipt_normal, std::shared_ptr<CBxDF> ipt_bxdf)
        :normal(ipt_normal)
        , bxdf(ipt_bxdf)
    {
        tangent_basis = CTangentBasis::fromZ(normal);
        assert(ipt_bxdf.get() != nullptr);
    };

    glm::vec3 f(glm::vec3 wo, glm::vec3 wi, ETransportMode transport_mode = ETransportMode::TM_Radiance)const
    {
        glm::vec3 wo_local = tangent_basis.toLocal(wo);
        glm::vec3 wi_local = tangent_basis.toLocal(wi);

        if (wo_local.z == 0)
        {
            return glm::vec3(0, 0, 0);
        }
        return bxdf->f(wo_local, wi_local, transport_mode);
    }

    float pdf(glm::vec3 ipt_wo, glm::vec3 ipt_wi, ETransportMode transport_mode = ETransportMode::TM_Radiance, EBxDFReflTransFlags sample_flags = EBxDFReflTransFlags::BXDF_RF_All)const
    {
        glm::vec3 wo_local = tangent_basis.toLocal(ipt_wo);
        glm::vec3 wi_local = tangent_basis.toLocal(ipt_wi);
        if (wo_local.z == 0)
        {
            return 0;
        }
        return bxdf->pdf(wo_local, wi_local, transport_mode, sample_flags);
    }

    std::shared_ptr<SBSDFSample> sample_f(glm::vec3 wo_world, float u, glm::vec2 u2, ETransportMode transport_mode = ETransportMode::TM_Radiance, EBxDFReflTransFlags flags = EBxDFReflTransFlags::BXDF_RF_All)
    {
        glm::vec3 wo_local = tangent_basis.toLocal(wo_world);
        if (wo_local.z == 0 || !(bxdf->flags() & BXDF_RF_All))
        {
            return nullptr;
        }

        std::shared_ptr<SBSDFSample> bs = bxdf->sample_f(wo_local, u, u2, transport_mode, flags);
        bs->wi = tangent_basis.fromLocal(bs->wi);
        bs->wi = glm::normalize(bs->wi);
        return bs;
    }

    inline EBxDFFlags flags()const
    {
        return bxdf->flags();
    }

    inline bool isValid() const
    {
        return bxdf.get() != nullptr;
    }

private:
    CTangentBasis  tangent_basis;
    glm::vec3 normal;
    std::shared_ptr<CBxDF> bxdf;
};

