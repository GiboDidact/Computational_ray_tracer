// pbrt is Copyright(c) 1998-2020 Matt Pharr, Wenzel Jakob, and Greg Humphreys.
// The pbrt source code is licensed under the Apache License, Version 2.0.
// SPDX: Apache-2.0
//EDITED

#pragma once
#include "pbrt.h"

// A special present from windgi.h on Windows...
#ifdef RGB
#undef RGB
#endif  // RGB

namespace pbrt {

    // RGB Definition
    class RGB
    {
    public:
        // RGB Public Methods
        RGB() = default;

        RGB(float r, float g, float b) : r(r), g(g), b(b) {}

        RGB(glm::vec3 vec) : r(vec.x), g(vec.y), b(vec.z) {}


        RGB& operator+=(RGB s)
        {
            r += s.r;
            g += s.g;
            b += s.b;
            return *this;
        }

        RGB operator+(RGB s) const
        {
            RGB ret = *this;
            return ret += s;
        }

        RGB& operator-=(RGB s)
        {
            r -= s.r;
            g -= s.g;
            b -= s.b;
            return *this;
        }

        RGB operator-(RGB s) const
        {
            RGB ret = *this;
            return ret -= s;
        }

        friend RGB operator-(float a, RGB s) { return { a - s.r, a - s.g, a - s.b }; }


        RGB& operator*=(RGB s)
        {
            r *= s.r;
            g *= s.g;
            b *= s.b;
            return *this;
        }

        RGB operator*(RGB s) const
        {
            RGB ret = *this;
            return ret *= s;
        }

        RGB operator*(float a) const
        {
            //DCHECK(!IsNaN(a));
            return { a * r, a * g, a * b };
        }

        RGB& operator*=(float a)
        {
            //DCHECK(!IsNaN(a));
            r *= a;
            g *= a;
            b *= a;
            return *this;
        }

        friend RGB operator*(float a, RGB s) { return s * a; }


        RGB& operator/=(RGB s)
        {
            r /= s.r;
            g /= s.g;
            b /= s.b;
            return *this;
        }

        RGB operator/(RGB s) const
        {
            RGB ret = *this;
            return ret /= s;
        }

        RGB& operator/=(float a)
        {
            //DCHECK(!IsNaN(a));
            //DCHECK_NE(a, 0);
            r /= a;
            g /= a;
            b /= a;
            return *this;
        }

        RGB operator/(float a) const
        {
            RGB ret = *this;
            return ret /= a;
        }


        RGB operator-() const { return { -r, -g, -b }; }


        float Average() const { return (r + g + b) / 3; }


        bool operator==(RGB s) const { return r == s.r && g == s.g && b == s.b; }

        bool operator!=(RGB s) const { return r != s.r || g != s.g || b != s.b; }

        float operator[](int c) const
        {
            //DCHECK(c >= 0 && c < 3);
            if (c == 0)
                return r;
            else if (c == 1)
                return g;
            return b;
        }

        float& operator[](int c)
        {
            //DCHECK(c >= 0 && c < 3);
            if (c == 0)
                return r;
            else if (c == 1)
                return g;
            return b;
        }

        glm::vec3 getglm() { return glm::vec3(r, g, b); }

        std::string ToString() const;

        // RGB Public Members
        float r = 0, g = 0, b = 0;
    };


    inline RGB max(RGB a, RGB b)
    {
        return RGB(std::max(a.r, b.r), std::max(a.g, b.g), std::max(a.b, b.b));
    }

    inline RGB Lerp(float t, RGB s1, RGB s2)
    {
        return (1 - t) * s1 + t * s2;
    }

    // RGB Inline Functions
    template <typename U, typename V>
    inline RGB Clamp(RGB rgb, U min, V max)
    {
        return RGB(pbrt::Clamp(rgb.r, min, max), pbrt::Clamp(rgb.g, min, max),
            pbrt::Clamp(rgb.b, min, max));
    }

    inline RGB ClampZero(RGB rgb)
    {
        return RGB(std::max(0.0f, rgb.r), std::max(0.0f, rgb.g),
            std::max(0.0f, rgb.b));
    }

    // XYZ Definition
    class XYZ
    {
    public:
        // XYZ Public Methods
        XYZ() = default;

        XYZ(float X, float Y, float Z) : X(X), Y(Y), Z(Z) {}

        XYZ(glm::vec3 vec) : X(vec.x), Y(vec.y), Z(vec.z) {}

        float Average() const { return (X + Y + Z) / 3; }

        glm::vec3 getglm() { return glm::vec3(X, Y, Z); }

        glm::vec2 xy() const { return glm::vec2(X / (X + Y + Z), Y / (X + Y + Z)); }


        static XYZ FromxyY(glm::vec2 xy, float Y = 1)
        {
            if (xy.y == 0)
                return XYZ(0, 0, 0);
            return XYZ(xy.x * Y / xy.y, Y, (1 - xy.x - xy.y) * Y / xy.y);
        }


        XYZ& operator+=(const XYZ& s)
        {
            X += s.X;
            Y += s.Y;
            Z += s.Z;
            return *this;
        }

        XYZ operator+(const XYZ& s) const
        {
            XYZ ret = *this;
            return ret += s;
        }


        XYZ& operator-=(const XYZ& s)
        {
            X -= s.X;
            Y -= s.Y;
            Z -= s.Z;
            return *this;
        }

        XYZ operator-(const XYZ& s) const
        {
            XYZ ret = *this;
            return ret -= s;
        }

        friend XYZ operator-(float a, const XYZ& s) { return { a - s.X, a - s.Y, a - s.Z }; }


        XYZ& operator*=(const XYZ& s)
        {
            X *= s.X;
            Y *= s.Y;
            Z *= s.Z;
            return *this;
        }

        XYZ operator*(const XYZ& s) const
        {
            XYZ ret = *this;
            return ret *= s;
        }

        XYZ operator*(float a) const
        {
            //DCHECK(!IsNaN(a));
            return { a * X, a * Y, a * Z };
        }

        XYZ& operator*=(float a)
        {
            //DCHECK(!IsNaN(a));
            X *= a;
            Y *= a;
            Z *= a;
            return *this;
        }


        XYZ& operator/=(const XYZ& s)
        {
            X /= s.X;
            Y /= s.Y;
            Z /= s.Z;
            return *this;
        }

        XYZ operator/(const XYZ& s) const
        {
            XYZ ret = *this;
            return ret /= s;
        }

        XYZ& operator/=(float a)
        {
            //DCHECK(!IsNaN(a));
            //DCHECK_NE(a, 0);
            X /= a;
            Y /= a;
            Z /= a;
            return *this;
        }

        XYZ operator/(float a) const
        {
            XYZ ret = *this;
            return ret /= a;
        }


        XYZ operator-() const { return { -X, -Y, -Z }; }


        bool operator==(const XYZ& s) const { return X == s.X && Y == s.Y && Z == s.Z; }

        bool operator!=(const XYZ& s) const { return X != s.X || Y != s.Y || Z != s.Z; }

        float operator[](int c) const
        {
            //DCHECK(c >= 0 && c < 3);
            if (c == 0)
                return X;
            else if (c == 1)
                return Y;
            return Z;
        }

        float& operator[](int c)
        {
            //DCHECK(c >= 0 && c < 3);
            if (c == 0)
                return X;
            else if (c == 1)
                return Y;
            return Z;
        }

        std::string ToString() const;

        // XYZ Public Members
        float X = 0, Y = 0, Z = 0;
    };


    inline XYZ operator*(float a, const XYZ& s)
    {
        return s * a;
    }

    template <typename U, typename V>
    inline XYZ Clamp(const XYZ& xyz, U min, V max)
    {
        return XYZ(pbrt::Clamp(xyz.X, min, max), pbrt::Clamp(xyz.Y, min, max),
            pbrt::Clamp(xyz.Z, min, max));
    }

    inline XYZ ClampZero(const XYZ& xyz)
    {
        return XYZ(std::max(0.0f, xyz.X), std::max(0.0f, xyz.Y),
            std::max(0.0f, xyz.Z));
    }

    inline XYZ Lerp(float t, const XYZ& s1, const XYZ& s2)
    {
        return (1 - t) * s1 + t * s2;
    }


    // RGBSigmoidPolynomial Definition
    class RGBSigmoidPolynomial
    {
    public:
        // RGBSigmoidPolynomial Public Methods
        RGBSigmoidPolynomial() = default;

        RGBSigmoidPolynomial(float c0, float c1, float c2) : c0(c0), c1(c1), c2(c2) {}
       
        std::string ToString() const;

        float operator()(float lambda) const
        {
            return s(EvaluatePolynomial(lambda, c2, c1, c0));
        }

        float polynomial(float lambda) const
        {
            return EvaluatePolynomial(lambda, c2, c1, c0);
        }

        float MaxValue() const
        {
            float result = std::max((*this)(360), (*this)(830));
            float lambda = -c1 / (2 * c0);
            if (lambda >= 360 && lambda <= 830)
                result = std::max(result, (*this)(lambda));
            return result;
        }

        // RGBSigmoidPolynomial Private Methods

        static float s(float x)
        {
            if (std::isinf(x))
                return x > 0 ? 1 : 0;
            return .5f + x / (2 * std::sqrt(1 + (x * x)));
        };
    private:
        // RGBSigmoidPolynomial Private Members
        float c0, c1, c2;
    };

    
    // RGBToSpectrumTable Definition
    class RGBToSpectrumTable
    {
    public:
        // RGBToSpectrumTable Public Constants
        static constexpr int res = 64;
   
        // RGBToSpectrumTable Public Methods
        RGBToSpectrumTable(const float* zNodes, float* _coeffs)
            : zNodes(zNodes) , coeffs(_coeffs) {}

        RGBSigmoidPolynomial operator()(RGB rgb) const;

        static void Init();


        static RGBToSpectrumTable* sRGB;
       // static const RGBToSpectrumTable* DCI_P3;
       // static const RGBToSpectrumTable* Rec2020;
      //  static const RGBToSpectrumTable* ACES2065_1;

        std::string ToString() const;

    private:
        // RGBToSpectrumTable Private Members
        const float* zNodes;
        float* coeffs;    
    };
    

    /*
    // ColorEncoding Definitions
    class LinearColorEncoding;
    class sRGBColorEncoding;
    class GammaColorEncoding;

    class ColorEncoding
        : public TaggedPointer<LinearColorEncoding, sRGBColorEncoding, GammaColorEncoding> {
    public:
        using TaggedPointer::TaggedPointer;
        // ColorEncoding Interface
         inline void ToLinear(pstd::span<const uint8_t> vin,
            pstd::span<float> vout) const;
         inline void FromLinear(pstd::span<const float> vin,
            pstd::span<uint8_t> vout) const;

         inline float TofloatLinear(float v) const;

        std::string ToString() const;

        static const ColorEncoding Get(const std::string& name, Allocator alloc);

        static ColorEncoding Linear;
        static ColorEncoding sRGB;

        static void Init(Allocator alloc);
    };

    class LinearColorEncoding {
    public:

            void ToLinear(pstd::span<const uint8_t> vin, pstd::span<float> vout) const {
            //DCHECK_EQ(vin.size(), vout.size());
            for (size_t i = 0; i < vin.size(); ++i)
                vout[i] = vin[i] / 255.f;
        }


            float TofloatLinear(float v) const { return v; }


            void FromLinear(pstd::span<const float> vin, pstd::span<uint8_t> vout) const {
            //DCHECK_EQ(vin.size(), vout.size());
            for (size_t i = 0; i < vin.size(); ++i)
                vout[i] = uint8_t(Clamp(vin[i] * 255.f + 0.5f, 0, 255));
        }

        std::string ToString() const { return "[ LinearColorEncoding ]"; }
    };

    class sRGBColorEncoding {
    public:
        // sRGBColorEncoding Public Methods

            void ToLinear(pstd::span<const uint8_t> vin, pstd::span<float> vout) const;

            float TofloatLinear(float v) const;

            void FromLinear(pstd::span<const float> vin, pstd::span<uint8_t> vout) const;

        std::string ToString() const { return "[ sRGBColorEncoding ]"; }
    };

    class GammaColorEncoding {
    public:

            GammaColorEncoding(float gamma);


            void ToLinear(pstd::span<const uint8_t> vin, pstd::span<float> vout) const;

            float TofloatLinear(float v) const;

            void FromLinear(pstd::span<const float> vin, pstd::span<uint8_t> vout) const;

        std::string ToString() const;

    private:
        float gamma;
        pstd::array<float, 256> applyLUT;
        pstd::array<float, 1024> inverseLUT;
    };

    inline void ColorEncoding::ToLinear(pstd::span<const uint8_t> vin,
        pstd::span<float> vout) const {
        auto tolin = [&](auto ptr) { return ptr->ToLinear(vin, vout); };
        Dispatch(tolin);
    }

    inline float ColorEncoding::TofloatLinear(float v) const {
        auto tfl = [&](auto ptr) { return ptr->TofloatLinear(v); };
        return Dispatch(tfl);
    }

    inline void ColorEncoding::FromLinear(pstd::span<const float> vin,
        pstd::span<uint8_t> vout) const {
        auto fl = [&](auto ptr) { return ptr->FromLinear(vin, vout); };
        Dispatch(fl);
    }


        inline float LinearToSRGB(float value) {
        if (value <= 0.0031308f)
            return 12.92f * value;
        // Minimax polynomial approximation from enoki's color.h.
        float sqrtValue = SafeSqrt(value);
        float p = EvaluatePolynomial(sqrtValue, -0.0016829072605308378f, 0.03453868659826638f,
            0.7642611304733891f, 2.0041169284241644f,
            0.7551545191665577f, -0.016202083165206348f);
        float q = EvaluatePolynomial(sqrtValue, 4.178892964897981e-7f,
            -0.00004375359692957097f, 0.03467195408529984f,
            0.6085338522168684f, 1.8970238036421054f, 1.f);
        return p / q * value;
    }


        inline uint8_t LinearToSRGB8(float value, float dither = 0) {
        if (value <= 0)
            return 0;
        if (value >= 1)
            return 255;
        return Clamp(pstd::round(255.f * LinearToSRGB(value) + dither), 0, 255);
    }


        inline float SRGBToLinear(float value) {
        if (value <= 0.04045f)
            return value * (1 / 12.92f);
        // Minimax polynomial approximation from enoki's color.h.
        float p = EvaluatePolynomial(value, -0.0163933279112946f, -0.7386328024653209f,
            -11.199318357635072f, -47.46726633009393f,
            -36.04572663838034f);
        float q = EvaluatePolynomial(value, -0.004261480793199332f, -19.140923959601675f,
            -59.096406619244426f, -18.225745396846637f, 1.f);
        return p / q * value;
    }

    extern PBRT_CONST float SRGBToLinearLUT[256];


        inline float SRGB8ToLinear(uint8_t value) {
        return SRGBToLinearLUT[value];
    }

    // White Balance Definitions
    // clang-format off
    // These are the Bradford transformation matrices.
    const SquareMatrix<3> LMSFromXYZ(0.8951, 0.2664, -0.1614,
        -0.7502, 1.7135, 0.0367,
        0.0389, -0.0685, 1.0296);
    const SquareMatrix<3> XYZFromLMS(0.986993, -0.147054, 0.159963,
        0.432305, 0.51836, 0.0492912,
        -0.00852866, 0.0400428, 0.968487);
    // clang-format on

    inline SquareMatrix<3> WhiteBalance(glm::vec2 srcWhite, glm::vec2 targetWhite) {
        // Find LMS coefficients for source and target white
        XYZ srcXYZ = XYZ::FromxyY(srcWhite), dstXYZ = XYZ::FromxyY(targetWhite);
        auto srcLMS = LMSFromXYZ * srcXYZ, dstLMS = LMSFromXYZ * dstXYZ;

        // Return white balancing matrix for source and target white
        SquareMatrix<3> LMScorrect = SquareMatrix<3>::Diag(
            dstLMS[0] / srcLMS[0], dstLMS[1] / srcLMS[1], dstLMS[2] / srcLMS[2]);
        return XYZFromLMS * LMScorrect * LMSFromXYZ;
    }
    */

    // White Balance Definitions
    // clang-format off
    // These are the Bradford transformation matrices.
    const glm::mat3 LMSFromXYZ(glm::vec3(0.8951, -0.7502, 0.0389), 
                                glm::vec3(0.2664, 1.7135, -0.0685), 
                                glm::vec3(-0.1614, 0.0367, 1.0296));
   
    const glm::mat3 XYZFromLMS(glm::vec3(0.986993, 0.432305, -0.00852866),
                                 glm::vec3(-0.147054, 0.51836, 0.0400428),
                                 glm::vec3(0.159963, 0.0492912, 0.968487));
  
    // clang-format on

    inline glm::mat3 WhiteBalance(glm::vec2 srcWhite, glm::vec2 targetWhite)
    {
        // Find LMS coefficients for source and target white
        XYZ srcXYZ = XYZ::FromxyY(srcWhite), dstXYZ = XYZ::FromxyY(targetWhite);
        auto srcLMS = LMSFromXYZ * srcXYZ.getglm(), dstLMS = LMSFromXYZ * dstXYZ.getglm();

        // Return white balancing matrix for source and target white
        glm::mat3 LMScorrect = glm::mat3(glm::vec3(dstLMS[0] / srcLMS[0], 0, 0), 
                                         glm::vec3(0, dstLMS[1] / srcLMS[1], 0), 
                                         glm::vec3(0, 0, dstLMS[2] / srcLMS[2]));

        return XYZFromLMS * LMScorrect * LMSFromXYZ;
    }
}  // namespace pbrt

