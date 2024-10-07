// pbrt is Copyright(c) 1998-2020 Matt Pharr, Wenzel Jakob, and Greg Humphreys.
// The pbrt source code is licensed under the Apache License, Version 2.0.
// SPDX: Apache-2.0

#ifndef PBRT_PARSER_H
#define PBRT_PARSER_H

//A7x:[BEGIN]
#include <filesystem>
//A7x:[END]

#include <functional>
#include <map>
#include <memory>
#include <string>
#include <string_view>
#include <optional>
#include "paramdict.h"

namespace pbrt {

	// ParserTarget Definition
	class ParserTarget {
	public:
        // ParserTarget Interface
        virtual void Scale(float sx, float sy, float sz) = 0;

        virtual void Shape(const std::string& name, ParsedParameterVector params) = 0;

        virtual ~ParserTarget() {};

        virtual void Option(const std::string& name, const std::string& value) = 0;

        virtual void Identity() = 0;
        virtual void Translate(float dx, float dy, float dz) = 0;
        virtual void Rotate(float angle, float ax, float ay, float az) = 0;
        virtual void LookAt(float ex, float ey, float ez, float lx, float ly, float lz,
            float ux, float uy, float uz) = 0;
        virtual void ConcatTransform(float transform[16]) = 0;
        virtual void Transform(float* transform) = 0;
        virtual void CoordinateSystem(const std::string&) = 0;
        virtual void CoordSysTransform(const std::string&) = 0;
        virtual void ActiveTransformAll() = 0;
        virtual void ActiveTransformEndTime() = 0;
        virtual void ActiveTransformStartTime() = 0;
        virtual void TransformTimes(float start, float end) = 0;

        virtual void ColorSpace(const std::string& n) = 0;
        virtual void PixelFilter(const std::string& name, ParsedParameterVector params) = 0;
        virtual void Film(const std::string& type, ParsedParameterVector params) = 0;
        virtual void Accelerator(const std::string& name, ParsedParameterVector params) = 0;
        virtual void Integrator(const std::string& name, ParsedParameterVector params) = 0;
        virtual void Camera(const std::string&, ParsedParameterVector params) = 0;
        virtual void MakeNamedMedium(const std::string& name, ParsedParameterVector params) = 0;
        virtual void MediumInterface(const std::string& insideName,  const std::string& outsideName) = 0;
        virtual void Sampler(const std::string& name, ParsedParameterVector params) = 0;

        virtual void WorldBegin() = 0;
        virtual void AttributeBegin() = 0;
        virtual void AttributeEnd() = 0;
        virtual void Attribute(const std::string& target, ParsedParameterVector params) = 0;
        virtual void Texture(const std::string& name, const std::string& type,
            const std::string& texname, ParsedParameterVector params) = 0;
        virtual void Material(const std::string& name, ParsedParameterVector params) = 0;
        virtual void MakeNamedMaterial(const std::string& name, ParsedParameterVector params) = 0;
        virtual void NamedMaterial(const std::string& name) = 0;
        virtual void LightSource(const std::string& name, ParsedParameterVector params) = 0;
        virtual void AreaLightSource(const std::string& name, ParsedParameterVector params) = 0;
        virtual void ReverseOrientation() = 0;
        virtual void ObjectBegin(const std::string& name) = 0;
        virtual void ObjectEnd() = 0;
        virtual void ObjectInstance(const std::string& name) = 0;

        virtual void EndOfFiles() = 0;

		//A7x:[BEGIN]
		virtual void SetSearchPath(const std::filesystem::path searchpath) = 0;
		//A7x:[END]

	protected:
	};

	// Scene Parsing Declarations
	void ParseFile(ParserTarget* target, const std::string& filename);

	// Token Definition
	struct Token {
		Token() = default;
		Token(std::string_view token) : token(token){}
		std::string ToString() const;
		std::string_view token;
	};

	// Tokenizer Definition
	class Tokenizer {
	public:
		Tokenizer(void* ptr, size_t len,
			std::function<void(const char*)> errorCallback);

		~Tokenizer();

		static std::unique_ptr<Tokenizer> CreateFromFile(
			const std::string& filename,
			std::function<void(const char*)> errorCallback);
		
		std::optional<Token> Next();
	private:
		// Tokenizer Private Methods
		void CheckUTF(const void* ptr, int len) const;

		int getChar() {
			if (pos == end)
				return EOF;
			int ch = *pos++;
			return ch;
		}
		void ungetChar() {
			--pos;
		}

		// Tokenizer Private Members
		// This function is called if there is an error during lexing.
		std::function<void(const char*)> errorCallback;

		char* srcDataPtr;

		// Pointers to the current position in the file and one past the end of
		// the file.
		const char* pos, * end;

		// If there are escaped characters in the string, we can't just return
		// a std::string_view into the mapped file. In that case, we handle the
		// escaped characters and return a std::string_view to sEscaped.  (And
		// thence, std::string_views from previous calls to Next() must be invalid
		// after a subsequent call, since we may reuse sEscaped.)
		std::string sEscaped;
	};
}
#endif