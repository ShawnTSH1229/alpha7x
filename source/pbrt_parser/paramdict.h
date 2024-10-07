// pbrt is Copyright(c) 1998-2020 Matt Pharr, Wenzel Jakob, and Greg Humphreys.
// The pbrt source code is licensed under the Apache License, Version 2.0.
// SPDX: Apache-2.0

#ifndef PBRT_PARAMDICT_H
#define PBRT_PARAMDICT_H

#include <limits>
#include <map>
#include <memory>
#include <string>
#include <vector>
#include <array>

#include <glm/vec2.hpp>
#include <glm/vec3.hpp>

namespace pbrt {
    // ParsedParameter Definition
    class ParsedParameter {
    public:
        // ParsedParameter Public Methods
        ParsedParameter() {}

        void AddFloat(float v);
        void AddInt(int i);
        void AddString(std::string_view str);
        void AddBool(bool v);

        std::string ToString() const;

        // ParsedParameter Public Members
        std::string type, name;
        std::vector<float> floats;
        std::vector<int> ints;
        std::vector<std::string> strings;
        std::vector<uint8_t> bools;
        mutable bool lookedUp = false;
        bool mayBeUnused = false;
    };

    // ParsedParameterVector Definition
    using ParsedParameterVector = std::vector<ParsedParameter*>;

    // ParameterType Definition
    enum class ParameterType {
        Boolean,
        Float,
        Integer,
        Point2f,
        Vector2f,
        Point3f,
        Vector3f,
        Normal3f,
        Spectrum,
        String,
        Texture
    };

    template <ParameterType PT>
    struct ParameterTypeTraits {};

    // ParameterDictionary Definition
    class ParameterDictionary {
    public:
        ParameterDictionary() = default;
        ParameterDictionary(ParsedParameterVector params);

        float GetOneFloat(const std::string& name, float def) const;
        int GetOneInt(const std::string& name, int def) const;
        bool GetOneBool(const std::string& name, bool def) const;
        std::string GetOneString(const std::string& name, const std::string& def) const;

        glm::vec2 GetOnePoint2f(const std::string& name, glm::vec2 def) const;
        glm::vec2 GetOneVector2f(const std::string& name, glm::vec2 def) const;
        glm::vec3 GetOnePoint3f(const std::string& name, glm::vec3 def) const;
        glm::vec3 GetOneVector3f(const std::string& name, glm::vec3 def) const;
        glm::vec3 GetOneNormal3f(const std::string& name, glm::vec3 def) const;

        std::vector<float> GetFloatArray(const std::string& name) const;
        std::vector<int> GetIntArray(const std::string& name) const;
        std::vector<uint8_t> GetBoolArray(const std::string& name) const;

        std::vector< glm::vec2> GetPoint2fArray(const std::string& name) const;
        std::vector< glm::vec2> GetVector2fArray(const std::string& name) const;
        std::vector< glm::vec3> GetPoint3fArray(const std::string& name) const;
        std::vector< glm::vec3> GetVector3fArray(const std::string& name) const;
        std::vector< glm::vec3> GetNormal3fArray(const std::string& name) const;
        std::vector<std::string> GetStringArray(const std::string& name) const;

        glm::vec3 GetRGBColor(const std::string& name)const;

        inline ParsedParameterVector& getParameters() { return params; }
    private:

        // ParameterDictionary Private Methods
        template <ParameterType PT>
        typename ParameterTypeTraits<PT>::ReturnType lookupSingle(
            const std::string& name,
            typename ParameterTypeTraits<PT>::ReturnType defaultValue) const;

        template <ParameterType PT>
        std::vector<typename ParameterTypeTraits<PT>::ReturnType> lookupArray(
            const std::string& name) const;

        template <typename ReturnType, typename G, typename C>
        std::vector<ReturnType> lookupArray(const std::string& name, ParameterType type,
            const char* typeName, int nPerItem, G getValues,
            C convert) const;

        // ParameterDictionary Private Members
        ParsedParameterVector params;
    };
};



#endif