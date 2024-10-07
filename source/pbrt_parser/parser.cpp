// pbrt is Copyright(c) 1998-2020 Matt Pharr, Wenzel Jakob, and Greg Humphreys.
// The pbrt source code is licensed under the Apache License, Version 2.0.
// SPDX: Apache-2.0

#include "parser.h"

#include <fstream>
#include <assert.h>
#include <format>
#include <double-conversion/double-conversion.h>

//A7X:[BEGIN]
#include <filesystem>
//A7X:[END]

namespace pbrt
{
    ///////////////////////////////////////////////////////////////////////////
// ParsedParameter

    void ParsedParameter::AddFloat(float v) {
        assert(ints.empty() && strings.empty() && bools.empty());
        floats.push_back(v);
    }

    void ParsedParameter::AddInt(int i) {
        assert(floats.empty() && strings.empty() && bools.empty());
        ints.push_back(i);
    }

    void ParsedParameter::AddString(std::string_view str) {
        assert(floats.empty() && ints.empty() && bools.empty());
        strings.push_back({ str.begin(), str.end() });
    }

    void ParsedParameter::AddBool(bool v) {
        assert(floats.empty() && ints.empty() && strings.empty());
        bools.push_back(v);
    }

    std::string ParsedParameter::ToString() const {
        std::string str;
        str += std::string("\"") + type + " " + name + std::string("\" [ ");
        if (!floats.empty())
            for (float d : floats)
                str += std::format("%f ", d);
        else if (!ints.empty())
            for (int i : ints)
                str += std::format("%d ", i);
        else if (!strings.empty())
            for (const auto& s : strings)
                str += '\"' + s + "\" ";
        else if (!bools.empty())
            for (bool b : bools)
                str += b ? "true " : "false ";
        str += "] ";

        return str;
    }

    static std::string toString(std::string_view s) {
        return std::string(s.data(), s.size());
    }

    // Tokenizer Implementation
    static char decodeEscaped(int ch) {
        switch (ch) {
        case EOF:
            assert(false);
            break;
        case 'b':
            return '\b';
        case 'f':
            return '\f';
        case 'n':
            return '\n';
        case 'r':
            return '\r';
        case 't':
            return '\t';
        case '\\':
            return '\\';
        case '\'':
            return '\'';
        case '\"':
            return '\"';
        default:
            assert(false);
        }
        return 0;  // NOTREACHED
    }

    static double_conversion::StringToDoubleConverter floatParser(
        double_conversion::StringToDoubleConverter::ALLOW_HEX, 0. /* empty string value */,
        0. /* junk string value */, nullptr /* infinity symbol */, nullptr /* NaN symbol */);

   

    void Tokenizer::CheckUTF(const void* ptr, int len) const {
        const unsigned char* c = (const unsigned char*)ptr;
        // https://en.wikipedia.org/wiki/Byte_order_mark
        if (len >= 2 && ((c[0] == 0xfe && c[1] == 0xff) || (c[0] == 0xff && c[1] == 0xfe)))
            errorCallback("File is encoded with UTF-16, which is not currently "
                "supported by pbrt (https://github.com/mmp/pbrt-v4/issues/136).");
    }


    Tokenizer::Tokenizer(void* ptr, size_t len, std::function<void(const char*)> errorCallback)
    {
        srcDataPtr = (char*)ptr;
        pos = (const char*)ptr;
        end = pos + len;
        CheckUTF(ptr, len);
    }

    Tokenizer::~Tokenizer()
    {
        free(srcDataPtr);
    }

    std::unique_ptr<Tokenizer> Tokenizer::CreateFromFile(const std::string& filename, std::function<void(const char*)> errorCallback)
    {
        std::ifstream file(filename, std::ios::ate);
        if (!file.is_open())
        {
            return std::unique_ptr<Tokenizer>(nullptr);
        }

        std::streamsize size = file.tellg();
        file.seekg(0, std::ios::beg);
        char* data = (char*)malloc(size);
        file.read(data, size);
        return std::make_unique<Tokenizer>(((void*)data), size_t(size), errorCallback);
    }

    std::optional<Token> Tokenizer::Next()
    {
        while (true) {
            const char* tokenStart = pos;

            int ch = getChar();
            if (ch == EOF)
                return std::nullopt;
            else if (ch == ' ' || ch == '\n' || ch == '\t' || ch == '\r') {
                // nothing
            }
            else if (ch == '"') {
                // scan to closing quote
                bool haveEscaped = false;
                while ((ch = getChar()) != '"') {
                    if (ch == EOF) {
                        errorCallback("premature EOF");
                        return {};
                    }
                    else if (ch == '\n') {
                        errorCallback("unterminated string");
                        return {};
                    }
                    else if (ch == '\\') {
                        haveEscaped = true;
                        // Grab the next character
                        if ((ch = getChar()) == EOF) {
                            errorCallback("premature EOF");
                            return {};
                        }
                    }
                }

                if (!haveEscaped)
                    return Token({ tokenStart, size_t(pos - tokenStart) });
                else {
                    sEscaped.clear();
                    for (const char* p = tokenStart; p < pos; ++p) {
                        if (*p != '\\')
                            sEscaped.push_back(*p);
                        else {
                            ++p;
                            sEscaped.push_back(decodeEscaped(*p));
                        }
                    }
                    return Token({ sEscaped.data(), sEscaped.size() });
                }
            } else if (ch == '[' || ch == ']') {
                return Token({ tokenStart, size_t(1) });
            } else if (ch == '#') {
                // comment: scan to EOL (or EOF)
                while ((ch = getChar()) != EOF) {
                    if (ch == '\n' || ch == '\r') {
                        ungetChar();
                        break;
                    }
                }

                return Token({ tokenStart, size_t(pos - tokenStart) });
            }
            else {
                // Regular statement or numeric token; scan until we hit a
                // space, opening quote, or bracket.
                while ((ch = getChar()) != EOF) {
                    if (ch == ' ' || ch == '\n' || ch == '\t' || ch == '\r' || ch == '"' ||
                        ch == '[' || ch == ']') {
                        ungetChar();
                        break;
                    }
                }
                return Token({ tokenStart, size_t(pos - tokenStart) });
            }
        };

        return std::nullopt;
    }

    static int parseInt(const Token& t) {
        bool negate = t.token[0] == '-';

        int index = 0;
        if (t.token[0] == '+' || t.token[0] == '-')
            ++index;

        int64_t value = 0;
        while (index < t.token.size()) {
            if (!(t.token[index] >= '0' && t.token[index] <= '9'))
                assert(false);
            value = 10 * value + (t.token[index] - '0');
            ++index;

            if (value > std::numeric_limits<int>::max())
                assert(false);
            else if (value < std::numeric_limits<int>::lowest())
                assert(false);
        }

        return negate ? -value : value;
    }

    static double parseFloat(const Token& t) {
        // Fast path for a single digit
        if (t.token.size() == 1) {
            if (!(t.token[0] >= '0' && t.token[0] <= '9'))
                assert(false);
            return t.token[0] - '0';
        }

        // Copy to a buffer so we can NUL-terminate it, as strto[idf]() expect.
        char buf[64];
        char* bufp = buf;
        std::unique_ptr<char[]> allocBuf;
        if (t.token.size() + 1 >= sizeof(buf)) {
            // This should be very unusual, but is necessary in case we get a
            // goofball number with lots of leading zeros, for example.
            allocBuf = std::make_unique<char[]>(t.token.size() + 1);
            bufp = allocBuf.get();
        }

        std::copy(t.token.begin(), t.token.end(), bufp);
        bufp[t.token.size()] = '\0';

        // Can we just use strtol?
        auto isInteger = [](std::string_view str) {
            for (char ch : str)
                if (!(ch >= '0' && ch <= '9'))
                    return false;
            return true;
            };

        int length = 0;
        double val;
        if (isInteger(t.token)) {
            char* endptr;
            val = double(strtol(bufp, &endptr, 10));
            length = endptr - bufp;
        }
        else if (sizeof(float) == sizeof(float))
            val = floatParser.StringToFloat(bufp, t.token.size(), &length);
        else
            val = floatParser.StringToDouble(bufp, t.token.size(), &length);

        if (length == 0)
            assert(false);

        return val;
    }

    inline bool isQuotedString(std::string_view str) {
        return str.size() >= 2 && str[0] == '"' && str.back() == '"';
    }

    static std::string_view dequoteString(const Token& t) {
        if (!isQuotedString(t.token))
            assert(false);

        std::string_view str = t.token;
        str.remove_prefix(1);
        str.remove_suffix(1);
        return str;
    }

    constexpr int TokenOptional = 0;
    constexpr int TokenRequired = 1;

    template <typename Next, typename Unget>
    static ParsedParameterVector parseParameters(
        Next nextToken, Unget ungetToken, bool formatting,
        const std::function<void(const Token& token, const char*)>& errorCallback) {
        ParsedParameterVector parameterVector;

        while (true) {
            std::optional<Token> t = nextToken(TokenOptional);
            if (!t.has_value())
                return parameterVector;

            if (!isQuotedString(t->token)) {
                ungetToken(*t);
                return parameterVector;
            }

            ParsedParameter* param = new ParsedParameter();

            std::string_view decl = dequoteString(*t);

            auto skipSpace = [&decl](std::string_view::const_iterator iter) {
                while (iter != decl.end() && (*iter == ' ' || *iter == '\t'))
                    ++iter;
                return iter;
                };
            // Skip to the next whitespace character (or the end of the string).
            auto skipToSpace = [&decl](std::string_view::const_iterator iter) {
                while (iter != decl.end() && *iter != ' ' && *iter != '\t')
                    ++iter;
                return iter;
                };

            auto typeBegin = skipSpace(decl.begin());
            if (typeBegin == decl.end())
                assert(false);

            // Find end of type declaration
            auto typeEnd = skipToSpace(typeBegin);
            param->type.assign(typeBegin, typeEnd);

            if (formatting) {  // close enough: upgrade...
                if (param->type == "point")
                    param->type = "point3";
                if (param->type == "vector")
                    param->type = "vector3";
                if (param->type == "color")
                    param->type = "rgb";
            }

            auto nameBegin = skipSpace(typeEnd);
            if (nameBegin == decl.end())
                assert(false);

            auto nameEnd = skipToSpace(nameBegin);
            param->name.assign(nameBegin, nameEnd);

            enum ValType { Unknown, String, Bool, Float, Int } valType = Unknown;

            if (param->type == "integer")
                valType = Int;

            auto addVal = [&](const Token& t) {
                if (isQuotedString(t.token)) {
                    switch (valType) {
                    case Unknown:
                        valType = String;
                        break;
                    case String:
                        break;
                    case Float:
                        errorCallback(t, "expected floating-point value");
                    case Int:
                        errorCallback(t, "expected integer value");
                    case Bool:
                        errorCallback(t, "expected Boolean value");
                    }

                    param->AddString(dequoteString(t));
                }
                else if (t.token[0] == 't' && t.token == "true") {
                    switch (valType) {
                    case Unknown:
                        valType = Bool;
                        break;
                    case String:
                        errorCallback(t, "expected string value");
                    case Float:
                        errorCallback(t, "expected floating-point value");
                    case Int:
                        errorCallback(t, "expected integer value");
                    case Bool:
                        break;
                    }

                    param->AddBool(true);
                }
                else if (t.token[0] == 'f' && t.token == "false") {
                    switch (valType) {
                    case Unknown:
                        valType = Bool;
                        break;
                    case String:
                        errorCallback(t, "expected string value");
                    case Float:
                        errorCallback(t, "expected floating-point value");
                    case Int:
                        errorCallback(t, "expected integer value");
                    case Bool:
                        break;
                    }

                    param->AddBool(false);
                }
                else {
                    switch (valType) {
                    case Unknown:
                        valType = Float;
                        break;
                    case String:
                        errorCallback(t, "expected string value");
                    case Float:
                        break;
                    case Int:
                        break;
                    case Bool:
                        errorCallback(t, "expected Boolean value");
                    }

                    if (valType == Int)
                        param->AddInt(parseInt(t));
                    else
                        param->AddFloat(parseFloat(t));
                }
                };

            Token val = *nextToken(TokenRequired);

            if (val.token == "[") {
                while (true) {
                    val = *nextToken(TokenRequired);
                    if (val.token == "]")
                        break;
                    addVal(val);
                }
            }
            else {
                addVal(val);
            }

            if (formatting && param->type == "bool") {
                for (const auto& b : param->strings) {
                    if (b == "true")
                        param->bools.push_back(true);
                    else if (b == "false")
                        param->bools.push_back(false);
                    else
                        assert(false);
                }
                param->strings.clear();
            }

            parameterVector.push_back(param);
        }

        return parameterVector;
    }

    void parse(ParserTarget* target, std::unique_ptr<pbrt::Tokenizer> t)
    {
        std::unique_ptr<Tokenizer> fileTokenizer = std::move(t);

        std::optional<Token> ungetToken;

        // nextToken is a little helper function that handles the file stack,
        // returning the next token from the current file until reaching EOF,
        // at which point it switches to the next file (if any).
        std::function<std::optional<Token>(int)> nextToken;
        nextToken = [&](int flags) -> std::optional<Token> {
            if (ungetToken.has_value())
                return std::exchange(ungetToken, {});

            std::optional<Token> tok = fileTokenizer->Next();

            if (!tok) {
                return std::optional<Token>();
            }
            else if (tok->token[0] == '#') {
                return nextToken(flags);
            }
            else
                // Regular token; success.
                return tok;
            };

        auto unget = [&](Token t) {
            assert(!ungetToken.has_value());
            ungetToken = t;
            };

        bool formatting = false;
        // Helper function for pbrt API entrypoints that take a single string
        // parameter and a ParameterVector (e.g. pbrtShape()).
        auto basicParamListEntrypoint =
            [&](void (ParserTarget::* apiFunc)(const std::string&, ParsedParameterVector)) {
            Token t = *nextToken(TokenRequired);
            std::string_view dequoted = dequoteString(t);
            std::string n = toString(dequoted);
            ParsedParameterVector parameterVector = parseParameters(
                nextToken, unget, formatting, [&](const Token& t, const char* msg) {
                    std::string token = toString(t.token);
                    printf(token.c_str());
                });
            (target->*apiFunc)(n, std::move(parameterVector));
            };

        auto syntaxError = [&](const Token& t) {
            printf("Unknown directive: %s", toString(t.token).c_str());
            };

        std::optional<Token> tok;

        while (true) {
            tok = nextToken(TokenOptional);
            if (!tok.has_value())
                break;

            switch (tok->token[0]) {
            case 'A':
                if (tok->token == "AttributeBegin")
                    target->AttributeBegin();
                else if (tok->token == "AttributeEnd")
                    target->AttributeEnd();
                else if (tok->token == "Attribute")
                    basicParamListEntrypoint(&ParserTarget::Attribute);
                else if (tok->token == "ActiveTransform") {
                    Token a = *nextToken(TokenRequired);
                    if (a.token == "All")
                        target->ActiveTransformAll();
                    else if (a.token == "EndTime")
                        target->ActiveTransformEndTime();
                    else if (a.token == "StartTime")
                        target->ActiveTransformStartTime();
                    else
                        syntaxError(*tok);
                }
                else if (tok->token == "AreaLightSource")
                    basicParamListEntrypoint(&ParserTarget::AreaLightSource);
                else if (tok->token == "Accelerator")
                    basicParamListEntrypoint(&ParserTarget::Accelerator);
                else
                    syntaxError(*tok);
                break;
            case 'C':
                if (tok->token == "ConcatTransform") {
                    if (nextToken(TokenRequired)->token != "[")
                        syntaxError(*tok);
                    float m[16];
                    for (int i = 0; i < 16; ++i)
                        m[i] = parseFloat(*nextToken(TokenRequired));
                    if (nextToken(TokenRequired)->token != "]")
                        syntaxError(*tok);
                    target->ConcatTransform(m);
                }
                else if (tok->token == "CoordinateSystem") {
                    std::string_view n = dequoteString(*nextToken(TokenRequired));
                    target->CoordinateSystem(toString(n));
                }
                else if (tok->token == "CoordSysTransform") {
                    std::string_view n = dequoteString(*nextToken(TokenRequired));
                    target->CoordSysTransform(toString(n));
                }
                else if (tok->token == "ColorSpace") {
                    std::string_view n = dequoteString(*nextToken(TokenRequired));
                    target->ColorSpace(toString(n));
                }
                else if (tok->token == "Camera")
                    basicParamListEntrypoint(&ParserTarget::Camera);
                else
                    syntaxError(*tok);
                break;

            case 'F':
                if (tok->token == "Film")
                    basicParamListEntrypoint(&ParserTarget::Film);
                else
                    syntaxError(*tok);
                break;

            case 'I':
                if (tok->token == "Integrator")
                    basicParamListEntrypoint(&ParserTarget::Integrator);
                else if (tok->token == "Identity")
                    target->Identity();
                else
                    syntaxError(*tok);
                break;

            case 'L':
                if (tok->token == "LightSource")
                    basicParamListEntrypoint(&ParserTarget::LightSource);
                else if (tok->token == "LookAt") {
                    float v[9];
                    for (int i = 0; i < 9; ++i)
                        v[i] = parseFloat(*nextToken(TokenRequired));
                    target->LookAt(v[0], v[1], v[2], v[3], v[4], v[5], v[6], v[7], v[8]);
                }
                else
                    syntaxError(*tok);
                break;

            case 'M':
                if (tok->token == "MakeNamedMaterial")
                    basicParamListEntrypoint(&ParserTarget::MakeNamedMaterial);
                else if (tok->token == "MakeNamedMedium")
                    basicParamListEntrypoint(&ParserTarget::MakeNamedMedium);
                else if (tok->token == "Material")
                    basicParamListEntrypoint(&ParserTarget::Material);
                else if (tok->token == "MediumInterface") {
                    std::string_view n = dequoteString(*nextToken(TokenRequired));
                    std::string names[2];
                    names[0] = toString(n);

                    // Check for optional second parameter
                    std::optional<Token> second = nextToken(TokenOptional);
                    if (second.has_value()) {
                        if (isQuotedString(second->token))
                            names[1] = toString(dequoteString(*second));
                        else {
                            unget(*second);
                            names[1] = names[0];
                        }
                    }
                    else
                        names[1] = names[0];

                    target->MediumInterface(names[0], names[1]);
                }
                else
                    syntaxError(*tok);
                break;

            case 'N':
                if (tok->token == "NamedMaterial") {
                    std::string_view n = dequoteString(*nextToken(TokenRequired));
                    target->NamedMaterial(toString(n));
                }
                else
                    syntaxError(*tok);
                break;

            case 'O':
                if (tok->token == "ObjectBegin") {
                    std::string_view n = dequoteString(*nextToken(TokenRequired));
                    target->ObjectBegin(toString(n));
                }
                else if (tok->token == "ObjectEnd")
                    target->ObjectEnd();
                else if (tok->token == "ObjectInstance") {
                    std::string_view n = dequoteString(*nextToken(TokenRequired));
                    target->ObjectInstance(toString(n));
                }
                else if (tok->token == "Option") {
                    std::string name = toString(dequoteString(*nextToken(TokenRequired)));
                    std::string value = toString(nextToken(TokenRequired)->token);
                    target->Option(name, value);
                }
                else
                    syntaxError(*tok);
                break;

            case 'P':
                if (tok->token == "PixelFilter")
                    basicParamListEntrypoint(&ParserTarget::PixelFilter);
                else
                    syntaxError(*tok);
                break;

            case 'R':
                if (tok->token == "ReverseOrientation")
                    target->ReverseOrientation();
                else if (tok->token == "Rotate") {
                    float v[4];
                    for (int i = 0; i < 4; ++i)
                        v[i] = parseFloat(*nextToken(TokenRequired));
                    target->Rotate(v[0], v[1], v[2], v[3]);
                }
                else
                    syntaxError(*tok);
                break;

            case 'S':
                if (tok->token == "Shape")
                    basicParamListEntrypoint(&ParserTarget::Shape);
                else if (tok->token == "Sampler")
                    basicParamListEntrypoint(&ParserTarget::Sampler);
                else if (tok->token == "Scale") {
                    float v[3];
                    for (int i = 0; i < 3; ++i)
                        v[i] = parseFloat(*nextToken(TokenRequired));
                    target->Scale(v[0], v[1], v[2]);
                }
                else
                    syntaxError(*tok);
                break;

            case 'T':
                if (tok->token == "TransformBegin") {
                    target->AttributeBegin();
                }
                else if (tok->token == "TransformEnd") {
                    target->AttributeEnd();
                }
                else if (tok->token == "Transform") {
                    if (nextToken(TokenRequired)->token != "[")
                        syntaxError(*tok);
                    float m[16];
                    for (int i = 0; i < 16; ++i)
                        m[i] = parseFloat(*nextToken(TokenRequired));
                    if (nextToken(TokenRequired)->token != "]")
                        syntaxError(*tok);
                    target->Transform(m);
                }
                else if (tok->token == "Translate") {
                    float v[3];
                    for (int i = 0; i < 3; ++i)
                        v[i] = parseFloat(*nextToken(TokenRequired));
                    target->Translate(v[0], v[1], v[2]);
                }
                else if (tok->token == "TransformTimes") {
                    float v[2];
                    for (int i = 0; i < 2; ++i)
                        v[i] = parseFloat(*nextToken(TokenRequired));
                    target->TransformTimes(v[0], v[1]);
                }
                else if (tok->token == "Texture") {
                    std::string_view n = dequoteString(*nextToken(TokenRequired));
                    std::string name = toString(n);
                    n = dequoteString(*nextToken(TokenRequired));
                    std::string type = toString(n);

                    Token t = *nextToken(TokenRequired);
                    std::string_view dequoted = dequoteString(t);
                    std::string texName = toString(dequoted);
                    ParsedParameterVector params = parseParameters(
                        nextToken, unget, formatting, [&](const Token& t, const char* msg) {
                            std::string token = toString(t.token);
                            std::string str = std::format("%s: %s", token, msg);
                            printf(str.c_str());
                        });

                    target->Texture(name, type, texName, std::move(params));
                }
                else
                    syntaxError(*tok);
                break;

            case 'W':
                if (tok->token == "WorldBegin")
                    target->WorldBegin();
                else if (tok->token == "WorldEnd" && formatting)
                    ;  // just swallow it
                else
                    syntaxError(*tok);
                break;

            default:
                break;
                syntaxError(*tok);
            };
                    
        };
    }

    void ParseFile(ParserTarget* target, const std::string& filename)
    {
        auto tokError = [](const char* msg) {
            printf("%s\n", msg);
            };

         //A7X:[BEGIN]
        std:: filesystem::path path(filename);
        if (!std::filesystem::is_directory(path))
        {
            path = path.parent_path();
        }
        //A7X:[END]

        target->SetSearchPath(path.string());
        std::unique_ptr<Tokenizer> t = Tokenizer::CreateFromFile(filename, tokError);
        parse(target, std::move(t));
        target->EndOfFiles();
    }

    std::string Token::ToString() const
    {
        return std::string(token.data(), token.size());
    }
}

