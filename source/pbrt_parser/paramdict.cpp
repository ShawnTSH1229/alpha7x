// pbrt is Copyright(c) 1998-2020 Matt Pharr, Wenzel Jakob, and Greg Humphreys.
// The pbrt source code is licensed under the Apache License, Version 2.0.
// SPDX: Apache-2.0

#include "paramdict.h"

namespace pbrt
{
    template <>
    struct ParameterTypeTraits<ParameterType::Boolean> {
        static constexpr char typeName[] = "bool";
        static constexpr int nPerItem = 1;
        using ReturnType = uint8_t;
        static bool Convert(const uint8_t* v) { return *v; }
        static const auto& GetValues(const ParsedParameter& param) { return param.bools; }
    };

    constexpr char ParameterTypeTraits<ParameterType::Boolean>::typeName[];

    template <>
    struct ParameterTypeTraits<ParameterType::Float> {
        static constexpr char typeName[] = "float";
        static constexpr int nPerItem = 1;
        using ReturnType = float;
        static float Convert(const float* v) { return *v; }
        static const auto& GetValues(const ParsedParameter& param) { return param.floats; }
    };

    constexpr char ParameterTypeTraits<ParameterType::Float>::typeName[];

    template <>
    struct ParameterTypeTraits<ParameterType::Integer> {
        static constexpr char typeName[] = "integer";
        static constexpr int nPerItem = 1;
        using ReturnType = int;
        static int Convert(const int* i) { return *i; }
        static const auto& GetValues(const ParsedParameter& param) { return param.ints; }
    };

    constexpr char ParameterTypeTraits<ParameterType::Integer>::typeName[];

    template <>
    struct ParameterTypeTraits<ParameterType::Point2f> {
        static constexpr char typeName[] = "point2";
        static constexpr int nPerItem = 2;
        using ReturnType = glm::vec2;
        static glm::vec2 Convert(const float* v) {
            return glm::vec2(v[0], v[1]);
        }
        static const auto& GetValues(const ParsedParameter& param) { return param.floats; }
    };

    constexpr char ParameterTypeTraits<ParameterType::Point2f>::typeName[];

    template <>
    struct ParameterTypeTraits<ParameterType::Vector2f> {
        static constexpr char typeName[] = "vector2";
        static constexpr int nPerItem = 2;
        using ReturnType = glm::vec2;
        static glm::vec2 Convert(const float* v) {
            return glm::vec2(v[0], v[1]);
        }
        static const auto& GetValues(const ParsedParameter& param) { return param.floats; }
    };

    constexpr char ParameterTypeTraits<ParameterType::Vector2f>::typeName[];

    // Point3f ParameterTypeTraits Definition
    template <>
    struct ParameterTypeTraits<ParameterType::Point3f> {
        // ParameterType::Point3f Type Traits
        using ReturnType = glm::vec3;

        static constexpr char typeName[] = "point3";

        static const auto& GetValues(const ParsedParameter& param) { return param.floats; }

        static constexpr int nPerItem = 3;

        static glm::vec3 Convert(const float* f) {
            return glm::vec3(f[0], f[1], f[2]);
        }
    };

    constexpr char ParameterTypeTraits<ParameterType::Point3f>::typeName[];

    template <>
    struct ParameterTypeTraits<ParameterType::Vector3f> {
        static constexpr char typeName[] = "vector3";
        static constexpr int nPerItem = 3;
        using ReturnType = glm::vec3;
        static glm::vec3 Convert(const float* v) {
            return glm::vec3(v[0], v[1], v[2]);
        }
        static const auto& GetValues(const ParsedParameter& param) { return param.floats; }
    };

    constexpr char ParameterTypeTraits<ParameterType::Vector3f>::typeName[];

    template <>
    struct ParameterTypeTraits<ParameterType::Normal3f> {
        static constexpr char typeName[] = "normal";
        static constexpr int nPerItem = 3;
        using ReturnType = glm::vec3;
        static glm::vec3 Convert(const float* v) {
            return glm::vec3(v[0], v[1], v[2]);
        }
        static const auto& GetValues(const ParsedParameter& param) { return param.floats; }
    };

    constexpr char ParameterTypeTraits<ParameterType::Normal3f>::typeName[];

    template <>
    struct ParameterTypeTraits<ParameterType::String> {
        static constexpr char typeName[] = "string";
        static constexpr int nPerItem = 1;
        using ReturnType = std::string;
        static std::string Convert(const std::string* s) { return *s; }
        static const auto& GetValues(const ParsedParameter& param) { return param.strings; }
    };

    constexpr char ParameterTypeTraits<ParameterType::String>::typeName[];

pbrt::ParameterDictionary::ParameterDictionary(ParsedParameterVector p)
	: params(std::move(p))
{
}

template<ParameterType PT>
inline typename ParameterTypeTraits<PT>::ReturnType pbrt::ParameterDictionary::lookupSingle(const std::string& name, typename ParameterTypeTraits<PT>::ReturnType defaultValue) const
{
    // Search _params_ for parameter _name_
    using traits = ParameterTypeTraits<PT>;
    for (const ParsedParameter* p : params) {
        if (p->name != name || p->type != traits::typeName)
            continue;
        // Extract parameter values from _p_
        const auto& values = traits::GetValues(*p);

        // Issue error if an incorrect number of parameter values were provided
        if (values.empty())
        {
            assert(false);
        }
        if (values.size() != traits::nPerItem)
        {
            assert(false);
        }

        // Return parameter values as _ReturnType_
        p->lookedUp = true;
        return traits::Convert(values.data());
    }

    return defaultValue;
}

float ParameterDictionary::GetOneFloat(const std::string& name, float def) const {
    return lookupSingle<ParameterType::Float>(name, def);
}

int ParameterDictionary::GetOneInt(const std::string& name, int def) const {
    return lookupSingle<ParameterType::Integer>(name, def);
}

bool ParameterDictionary::GetOneBool(const std::string& name, bool def) const {
    return lookupSingle<ParameterType::Boolean>(name, def);
}

glm::vec2 ParameterDictionary::GetOnePoint2f(const std::string& name, glm::vec2 def) const {
    return lookupSingle<ParameterType::Point2f>(name, def);
}

glm::vec2 ParameterDictionary::GetOneVector2f(const std::string& name,
    glm::vec2 def) const {
    return lookupSingle<ParameterType::Vector2f>(name, def);
}

glm::vec3 ParameterDictionary::GetOneVector3f(const std::string& name,
    glm::vec3 def) const {
    return lookupSingle<ParameterType::Vector3f>(name, def);
}

glm::vec3 ParameterDictionary::GetOneNormal3f(const std::string& name,
    glm::vec3 def) const {
    return lookupSingle<ParameterType::Normal3f>(name, def);
}

std::string ParameterDictionary::GetOneString(const std::string& name,
    const std::string& def) const {
    return lookupSingle<ParameterType::String>(name, def);
}

template <typename ReturnType, typename ValuesType, typename C>
static std::vector<ReturnType> returnArray(const ValuesType& values,
    const ParsedParameter& param, int nPerItem,
    C convert) {
    if (values.empty())
        assert(false);
    if (values.size() % nPerItem)
        assert(false);

    param.lookedUp = true;
    size_t n = values.size() / nPerItem;
    std::vector<ReturnType> v(n);
    for (size_t i = 0; i < n; ++i)
        v[i] = convert(&values[nPerItem * i]);
    return v;
}

template <typename ReturnType, typename G, typename C>
std::vector<ReturnType> ParameterDictionary::lookupArray(const std::string& name,
    ParameterType type,
    const char* typeName,
    int nPerItem, G getValues,
    C convert) const {
    for (const ParsedParameter* p : params)
        if (p->name == name && p->type == typeName)
            return returnArray<ReturnType>(getValues(*p), *p, nPerItem, convert);

    return {};
}

template <ParameterType PT>
std::vector<typename ParameterTypeTraits<PT>::ReturnType>
ParameterDictionary::lookupArray(const std::string& name) const {
    using traits = ParameterTypeTraits<PT>;
    return lookupArray<typename traits::ReturnType>(
        name, PT, traits::typeName, traits::nPerItem, traits::GetValues, traits::Convert);
}

std::vector<float> ParameterDictionary::GetFloatArray(const std::string& name) const {
    return lookupArray<ParameterType::Float>(name);
}

std::vector<int> ParameterDictionary::GetIntArray(const std::string& name) const {
    return lookupArray<ParameterType::Integer>(name);
}

std::vector<uint8_t> ParameterDictionary::GetBoolArray(const std::string& name) const {
    return lookupArray<ParameterType::Boolean>(name);
}

std::vector<glm::vec2> ParameterDictionary::GetPoint2fArray(const std::string& name) const {
    return lookupArray<ParameterType::Point2f>(name);
}

std::vector<glm::vec2> ParameterDictionary::GetVector2fArray(
    const std::string& name) const {
    return lookupArray<ParameterType::Vector2f>(name);
}

std::vector<glm::vec3> ParameterDictionary::GetPoint3fArray(const std::string& name) const {
    return lookupArray<ParameterType::Point3f>(name);
}

std::vector<glm::vec3> ParameterDictionary::GetVector3fArray(
    const std::string& name) const {
    return lookupArray<ParameterType::Vector3f>(name);
}

std::vector<glm::vec3> ParameterDictionary::GetNormal3fArray(
    const std::string& name) const {
    return lookupArray<ParameterType::Normal3f>(name);
}

std::vector<std::string> ParameterDictionary::GetStringArray(
    const std::string& name) const {
    return lookupArray<ParameterType::String>(name);
}

glm::vec3 ParameterDictionary::GetRGBColor(const std::string& name) const
{
    for (const ParsedParameter* p : params)
    {
        if (p->name == name && p->type == "rgb")
        {
            return glm::vec3(p->floats[0], p->floats[1], p->floats[2]);
        }
    }


    return glm::vec3();
}

}

