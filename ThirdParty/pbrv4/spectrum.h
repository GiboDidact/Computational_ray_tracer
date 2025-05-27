// pbrt is Copyright(c) 1998-2020 Matt Pharr, Wenzel Jakob, and Greg Humphreys.
// The pbrt source code is licensed under the Apache License, Version 2.0.
// SPDX: Apache-2.0
//EDITED

// PhysLight code contributed by Anders Langlands and Luca Fascione
// Copyright (c) 2020, Weta Digital, Ltd.
// SPDX-License-Identifier: Apache-2.0

#pragma once
#include "pbrt.h"
#include "color.h"
#include "../../RayTracer/Sampling.h"


namespace pbrt {
    constexpr float Lambda_min = 360, Lambda_max = 830;

    static constexpr int NSpectrumSamples = 8;

    static constexpr float CIE_Y_integral = 106.856895;

    // Spectrum Function Declarations
    static float Blackbody(float lambda, float T)
    {
        if (T <= 0)
            return 0;
        const float c = 299792458.f;
        const float h = 6.62606957e-34f;
        const float kb = 1.3806488e-23f;
        // Return emitted radiance for blackbody at wavelength _lambda_
        float l = lambda * 1e-9f;
        float Le = (2 * h * c * c) / (std::pow(l, 5) * (std::exp((h * c) / (l * kb * T)) - 1));
        //CHECK(!IsNaN(Le));
        return Le;
    }

    class Spectrum
    {
    public:
        // Spectrum Interface
        virtual std::string ToString() const = 0;

        virtual float Query(float lambda) const = 0;

        virtual float MaxValue() const = 0;

        virtual SampledSpectrum Sample(const SampledWavelengths& lambda) const = 0;
    };
    
    // SampledSpectrum Definition
    class SampledSpectrum
    {
    public:
        // SampledSpectrum Public Methods

        SampledSpectrum operator+(const SampledSpectrum& s) const
        {
            SampledSpectrum ret = *this;
            return ret += s;
        }

        SampledSpectrum& operator-=(const SampledSpectrum& s)
        {
            for (int i = 0; i < NSpectrumSamples; ++i)
                values[i] -= s.values[i];
            return *this;
        }

        SampledSpectrum operator-(const SampledSpectrum& s) const
        {
            SampledSpectrum ret = *this;
            return ret -= s;
        }

        friend SampledSpectrum operator-(float a, const SampledSpectrum& s)
        {
            //DCHECK(!IsNaN(a));
            SampledSpectrum ret;
            for (int i = 0; i < NSpectrumSamples; ++i)
                ret.values[i] = a - s.values[i];
            return ret;
        }


        SampledSpectrum& operator*=(const SampledSpectrum& s)
        {
            for (int i = 0; i < NSpectrumSamples; ++i)
                values[i] *= s.values[i];
            return *this;
        }

        SampledSpectrum operator*(const SampledSpectrum& s) const
        {
            SampledSpectrum ret = *this;
            return ret *= s;
        }

        SampledSpectrum operator*(float a) const
        {
            //DCHECK(!IsNaN(a));
            SampledSpectrum ret = *this;
            for (int i = 0; i < NSpectrumSamples; ++i)
                ret.values[i] *= a;
            return ret;
        }

        SampledSpectrum& operator*=(float a)
        {
            //DCHECK(!IsNaN(a));
            for (int i = 0; i < NSpectrumSamples; ++i)
                values[i] *= a;
            return *this;
        }

        friend SampledSpectrum operator*(float a, const SampledSpectrum& s) { return s * a; }


        SampledSpectrum& operator/=(const SampledSpectrum& s)
        {
            for (int i = 0; i < NSpectrumSamples; ++i)
            {
                //DCHECK_NE(0, s.values[i]);
                values[i] /= s.values[i];
            }
            return *this;
        }

        SampledSpectrum operator/(const SampledSpectrum& s) const
        {
            SampledSpectrum ret = *this;
            return ret /= s;
        }

        SampledSpectrum& operator/=(float a)
        {
            //DCHECK_NE(a, 0);
            //DCHECK(!IsNaN(a));
            for (int i = 0; i < NSpectrumSamples; ++i)
                values[i] /= a;
            return *this;
        }

        SampledSpectrum operator/(float a) const
        {
            SampledSpectrum ret = *this;
            return ret /= a;
        }


        SampledSpectrum operator-() const
        {
            SampledSpectrum ret;
            for (int i = 0; i < NSpectrumSamples; ++i)
                ret.values[i] = -values[i];
            return ret;
        }

        bool operator==(const SampledSpectrum& s) const { return values == s.values; }

        bool operator!=(const SampledSpectrum& s) const { return values != s.values; }

        std::string ToString() const;


        bool HasNaNs() const
        {
            for (int i = 0; i < NSpectrumSamples; ++i)
                if (std::isnan(values[i]))
                    return true;
            return false;
        }


        XYZ ToXYZ(const SampledWavelengths& lambda) const;

        RGB ToRGB(const SampledWavelengths& lambda, const RGBColorSpace& cs) const;

        float y(const SampledWavelengths& lambda) const;

        SampledSpectrum() = default;

        explicit SampledSpectrum(float c) { values.fill(c); }

        SampledSpectrum(std::span<const float> v)
        {
            //DCHECK_EQ(NSpectrumSamples, v.size());
            for (int i = 0; i < NSpectrumSamples; ++i)
                values[i] = v[i];
        }


        float operator[](int i) const
        {
            //DCHECK(i >= 0 && i < NSpectrumSamples);
            return values[i];
        }

        float& operator[](int i)
        {
            //DCHECK(i >= 0 && i < NSpectrumSamples);
            return values[i];
        }


        explicit operator bool() const
        {
            for (int i = 0; i < NSpectrumSamples; ++i)
                if (values[i] != 0)
                    return true;
            return false;
        }

        SampledSpectrum& operator+=(const SampledSpectrum& s)
        {
            for (int i = 0; i < NSpectrumSamples; ++i)
                values[i] += s.values[i];
            return *this;
        }


        float MinComponentValue() const
        {
            float m = values[0];
            for (int i = 1; i < NSpectrumSamples; ++i)
                m = std::min(m, values[i]);
            return m;
        }

        float MaxComponentValue() const
        {
            float m = values[0];
            for (int i = 1; i < NSpectrumSamples; ++i)
                m = std::max(m, values[i]);
            return m;
        }

        float Average() const
        {
            float sum = values[0];
            for (int i = 1; i < NSpectrumSamples; ++i)
                sum += values[i];
            return sum / NSpectrumSamples;
        }

    private:
        //friend struct SOA<SampledSpectrum>;
        std::array<float, NSpectrumSamples> values; //these are the reflectance scalar values between 0 and 1?
    };


    // SampledWavelengths Definitions
    class SampledWavelengths
    {
    public:
        // SampledWavelengths Public Methods

        bool operator==(const SampledWavelengths& swl) const
        {
            return lambda == swl.lambda && pdf == swl.pdf;
        }

        bool operator!=(const SampledWavelengths& swl) const
        {
            return lambda != swl.lambda || pdf != swl.pdf;
        }

        std::string ToString() const;

        static SampledWavelengths SampleUniform(float u, float lambda_min = Lambda_min,
            float lambda_max = Lambda_max)
        {
            SampledWavelengths swl;
            // Sample first wavelength using _u_
            swl.lambda[0] = Lerp(u, lambda_min, lambda_max);

            // Initialize _lambda_ for remaining wavelengths
            float delta = (lambda_max - lambda_min) / NSpectrumSamples;
            for (int i = 1; i < NSpectrumSamples; ++i)
            {
                swl.lambda[i] = swl.lambda[i - 1] + delta;
                if (swl.lambda[i] > lambda_max)
                    swl.lambda[i] = lambda_min + (swl.lambda[i] - lambda_max);
            }

            // Compute PDF for sampled wavelengths
            for (int i = 0; i < NSpectrumSamples; ++i)
                swl.pdf[i] = 1 / (lambda_max - lambda_min);

            return swl;
        }

        float operator[](int i) const { return lambda[i]; }

        float& operator[](int i) { return lambda[i]; }

        SampledSpectrum PDF() const
        {
            return SampledSpectrum(pdf);
        }

        void TerminateSecondary()
        {
            if (SecondaryTerminated())
                return;
            // Update wavelength probabilities for termination
            for (int i = 1; i < NSpectrumSamples; ++i)
                pdf[i] = 0;
            pdf[0] /= NSpectrumSamples;
        }


        bool SecondaryTerminated() const
        {
            for (int i = 1; i < NSpectrumSamples; ++i)
                if (pdf[i] != 0)
                    return false;
            return true;
        }


        static SampledWavelengths SampleVisible(float u)
        {
            SampledWavelengths swl;
            for (int i = 0; i < NSpectrumSamples; ++i)
            {
                // Compute _up_ for $i$th wavelength sample
                float up = u + float(i) / NSpectrumSamples;
                if (up > 1)
                    up -= 1;

                swl.lambda[i] = SampleVisibleWavelengths(up);
                swl.pdf[i] = VisibleWavelengthsPDF(swl.lambda[i]);
            }
            return swl;
        }

        void tamperlambda(int index, float val) { lambda[index] = val; }
    private:
        // SampledWavelengths Private Members
        //friend struct SOA<SampledWavelengths>;
        std::array<float, NSpectrumSamples> lambda, pdf;
    };


    class DenselySampledSpectrum;
    namespace Spectra {
        DenselySampledSpectrum D(float T);
    }   //namespace Spectra

    float SpectrumToPhotometric(Spectrum* s);

    XYZ SpectrumToXYZ(Spectrum* s);


// Spectrum Definitions
    class ConstantSpectrum : public Spectrum
    {
    public:

        ConstantSpectrum(float c) : c(c) {}

        float Query(float lambda) const override { return c; }
        // ConstantSpectrum Public Methods

        SampledSpectrum Sample(const SampledWavelengths&) const override;

        float MaxValue() const override { return c; }

        std::string ToString() const override;

    private:
        float c;
    };

    class DenselySampledSpectrum : public Spectrum
    {
    public:
        // DenselySampledSpectrum Public Methods
        DenselySampledSpectrum(int lambda_min = Lambda_min, int lambda_max = Lambda_max) :
            lambda_min(lambda_min), lambda_max(lambda_max), values(lambda_max - lambda_min + 1) {}
        //DenselySampledSpectrum(Spectrum* s) : DenselySampledSpectrum(s, Lambda_min, Lambda_max) {}
        DenselySampledSpectrum(const DenselySampledSpectrum& s) :
            lambda_min(s.lambda_min), lambda_max(s.lambda_max), values(s.values.begin(), s.values.end()) {}

        SampledSpectrum Sample(const SampledWavelengths& lambda) const override
        {
            SampledSpectrum s;
            for (int i = 0; i < NSpectrumSamples; ++i)
            {
                int offset = std::lround(lambda[i]) - lambda_min;
                if (offset < 0 || offset >= values.size())
                    s[i] = 0;
                else
                    s[i] = values[offset];
            }
            return s;
        }


        void Scale(float s)
        {
            for (float& v : values)
                v *= s;
        }


        float MaxValue() const override { return *std::max_element(values.begin(), values.end()); }

        std::string ToString() const override;

        DenselySampledSpectrum(Spectrum* spec, int lambda_min = Lambda_min, int lambda_max = Lambda_max) :
            lambda_min(lambda_min), lambda_max(lambda_max), values(lambda_max - lambda_min + 1)
        {
            //CHECK_GE(lambda_max, lambda_min);
            if (spec)
                for (int lambda = lambda_min; lambda <= lambda_max; ++lambda)
                    values[lambda - lambda_min] = spec->Query(lambda);
        }

        template <typename F>
        static DenselySampledSpectrum SampleFunction(F func, int lambda_min = Lambda_min, int lambda_max = Lambda_max)
        {
            DenselySampledSpectrum s(lambda_min, lambda_max);
            for (int lambda = lambda_min; lambda <= lambda_max; ++lambda)
                s.values[lambda - lambda_min] = func(lambda);
            return s;
        }

        float Query(float lambda) const override
        {
            //DCHECK_GT(lambda, 0);
            int offset = std::lround(lambda) - lambda_min;
            if (offset < 0 || offset >= values.size())
                return 0;
            return values[offset];
        }


        bool operator==(const DenselySampledSpectrum& d) const
        {
            if (lambda_min != d.lambda_min || lambda_max != d.lambda_max ||
                values.size() != d.values.size())
                return false;
            for (size_t i = 0; i < values.size(); ++i)
                if (values[i] != d.values[i])
                    return false;
            return true;
        }

    private:
        friend struct std::hash<pbrt::DenselySampledSpectrum>;
        // DenselySampledSpectrum Private Members
        int lambda_min, lambda_max;
        std::vector<float> values;
    };

    class PiecewiseLinearSpectrum : public Spectrum
    {
    public:
        // PiecewiseLinearSpectrum Public Methods
        PiecewiseLinearSpectrum() = default;

        void Scale(float s)
        {
            for (float& v : values)
                v *= s;
        }

        float MaxValue() const override;


        SampledSpectrum Sample(const SampledWavelengths& lambda) const override
        {
            SampledSpectrum s;
            for (int i = 0; i < NSpectrumSamples; ++i)
                s[i] = (*this).Query(lambda[i]);
            return s;
        }

        float Query(float lambda) const override;

        std::string ToString() const override;

        PiecewiseLinearSpectrum(std::span<const float> lambdas,
            std::span<const float> values);

        static std::optional<Spectrum*> Read(const std::string& filename);

        static PiecewiseLinearSpectrum* FromInterleaved(std::span<const float> samples,
            bool normalize);

    private:
        // PiecewiseLinearSpectrum Private Members
        std::vector<float> lambdas, values;
    };

    class BlackbodySpectrum : public Spectrum
    {
    public:
        // BlackbodySpectrum Public Methods

        BlackbodySpectrum(float T) : T(T)
        {
            // Compute blackbody normalization constant for given temperature
            float lambdaMax = 2.8977721e-3f / T;
            normalizationFactor = 1 / Blackbody(lambdaMax * 1e9f, T);
        }

        float Query(float lambda) const override
        {
            return Blackbody(lambda, T) * normalizationFactor;
        }

        SampledSpectrum Sample(const SampledWavelengths& lambda) const override
        {
            SampledSpectrum s;
            for (int i = 0; i < NSpectrumSamples; ++i)
                s[i] = Blackbody(lambda[i], T) * normalizationFactor;
            return s;
        }

        float MaxValue() const override { return 1.f; }

        std::string ToString() const override;

    private:
        // BlackbodySpectrum Private Members
        float T;
        float normalizationFactor;
    };


    //pass in color space and rgb and it gives you the reflectance(0,1) spectral distribution
    class RGBAlbedoSpectrum : public Spectrum
    {
    public:
        // RGBAlbedoSpectrum Public Methods

        float Query(float lambda) const override { return rsp(lambda); }

        float MaxValue() const override { return rsp.MaxValue(); }

        RGBAlbedoSpectrum(const RGBColorSpace& cs, RGB rgb);

        SampledSpectrum Sample(const SampledWavelengths& lambda) const override
        {
            SampledSpectrum s;
            for (int i = 0; i < NSpectrumSamples; ++i)
                s[i] = rsp(lambda[i]);
            return s;
        }

        std::string ToString() const override;

    private:
        // RGBAlbedoSpectrum Private Members
        RGBSigmoidPolynomial rsp;
    };

    class RGBUnboundedSpectrum : public Spectrum
    {
    public:
        // RGBUnboundedSpectrum Public Methods

        float Query(float lambda) const override { return scale * rsp(lambda); }

        float MaxValue() const override { return scale * rsp.MaxValue(); }


        RGBUnboundedSpectrum(const RGBColorSpace& cs, RGB rgb);


        RGBUnboundedSpectrum() : rsp(0, 0, 0), scale(0) {}


        SampledSpectrum Sample(const SampledWavelengths& lambda) const override
        {
            SampledSpectrum s;
            for (int i = 0; i < NSpectrumSamples; ++i)
                s[i] = scale * rsp(lambda[i]);
            return s;
        }

        std::string ToString() const override;

    private:
        // RGBUnboundedSpectrum Private Members
        float scale = 1;
        RGBSigmoidPolynomial rsp;
    };

    class RGBIlluminantSpectrum : public Spectrum
    {
    public:
        // RGBIlluminantSpectrum Public Methods
        RGBIlluminantSpectrum() = default;

        RGBIlluminantSpectrum(const RGBColorSpace& cs, RGB rgb);


        float Query(float lambda) const override
        {
            if (!illuminant)
                return 0;
            return scale * rsp(lambda) * (*illuminant).Query(lambda);
        }


        float MaxValue() const override
        {
            if (!illuminant)
                return 0;
            return scale * rsp.MaxValue() * illuminant->MaxValue();
        }


        const DenselySampledSpectrum* Illuminant() const { return illuminant; }


        SampledSpectrum Sample(const SampledWavelengths& lambda) const override
        {
            if (!illuminant)
                return SampledSpectrum(0);
            SampledSpectrum s;
            for (int i = 0; i < NSpectrumSamples; ++i)
                s[i] = scale * rsp(lambda[i]);
            return s * illuminant->Sample(lambda);
        }

        std::string ToString() const override;

    private:
        // RGBIlluminantSpectrum Private Members
        float scale;
        RGBSigmoidPolynomial rsp;
        const DenselySampledSpectrum* illuminant;
    };



    // SampledSpectrum Inline Functions
    inline SampledSpectrum SafeDiv(SampledSpectrum a, SampledSpectrum b)
    {
        SampledSpectrum r;
        for (int i = 0; i < NSpectrumSamples; ++i)
            r[i] = (b[i] != 0) ? a[i] / b[i] : 0.;
        return r;
    }

    template <typename U, typename V>
    inline SampledSpectrum Clamp(const SampledSpectrum& s, U low, V high)
    {
        SampledSpectrum ret;
        for (int i = 0; i < NSpectrumSamples; ++i)
            ret[i] = pbrt::Clamp(s[i], low, high);
        //DCHECK(!ret.HasNaNs());
        return ret;
    }

    inline SampledSpectrum ClampZero(const SampledSpectrum& s)
    {
        SampledSpectrum ret;
        for (int i = 0; i < NSpectrumSamples; ++i)
            ret[i] = std::max(0.0f, s[i]);
        //DCHECK(!ret.HasNaNs());
        return ret;
    }

    inline SampledSpectrum Sqrt(const SampledSpectrum& s)
    {
        SampledSpectrum ret;
        for (int i = 0; i < NSpectrumSamples; ++i)
            ret[i] = std::sqrt(s[i]);
        //DCHECK(!ret.HasNaNs());
        return ret;
    }

    inline SampledSpectrum SafeSqrt(const SampledSpectrum& s)
    {
        SampledSpectrum ret;
        for (int i = 0; i < NSpectrumSamples; ++i)
            ret[i] = SafeSqrt(s[i]);
        //DCHECK(!ret.HasNaNs());
        return ret;
    }

    inline SampledSpectrum Pow(const SampledSpectrum& s, float e)
    {
        SampledSpectrum ret;
        for (int i = 0; i < NSpectrumSamples; ++i)
            ret[i] = std::pow(s[i], e);
        return ret;
    }

    inline SampledSpectrum Exp(const SampledSpectrum& s)
    {
        SampledSpectrum ret;
        for (int i = 0; i < NSpectrumSamples; ++i)
            ret[i] = std::exp(s[i]);
        //DCHECK(!ret.HasNaNs());
        return ret;
    }

    inline SampledSpectrum FastExp(const SampledSpectrum& s)
    {
        SampledSpectrum ret;
        for (int i = 0; i < NSpectrumSamples; ++i)
            ret[i] = std::exp(s[i]);
        //DCHECK(!ret.HasNaNs());
        return ret;
    }

    inline SampledSpectrum Bilerp(std::array<float, 2> p, std::span<const SampledSpectrum> v)
    {
        return ((1 - p[0]) * (1 - p[1]) * v[0] + p[0] * (1 - p[1]) * v[1] +
            (1 - p[0]) * p[1] * v[2] + p[0] * p[1] * v[3]);
    }

    inline SampledSpectrum Lerp(float t, const SampledSpectrum& s1, const SampledSpectrum& s2)
    {
        return (1 - t) * s1 + t * s2;
    }

    // Spectral Data Declarations
    namespace Spectra {

        void Init();

        inline const DenselySampledSpectrum& X()
        {
            extern DenselySampledSpectrum* x;
            return *x;
        }

        inline const DenselySampledSpectrum& Y()
        {
            extern DenselySampledSpectrum* y;
            return *y;
        }

        inline const DenselySampledSpectrum& Z()
        {
            extern DenselySampledSpectrum* z;
            return *z;
        }

    }  // namespace Spectra

    // Spectral Function Declarations
    Spectrum* GetNamedSpectrum(std::string name);

    std::string FindMatchingNamedSpectrum(Spectrum s);

    namespace Spectra {
         inline const DenselySampledSpectrum& X();
         inline const DenselySampledSpectrum& Y();
         inline const DenselySampledSpectrum& Z();
    }  // namespace Spectra

    // Spectrum Inline Functions
    inline float InnerProduct(Spectrum* f, Spectrum* g)
    {
        float integral = 0;
        for (float lambda = Lambda_min; lambda <= Lambda_max; ++lambda)
            integral += f->Query(lambda) * g->Query(lambda);
        return integral;
    }

    inline float InnerProduct(const Spectrum* f, const Spectrum* g)
    {
        float integral = 0;
        for (float lambda = Lambda_min; lambda <= Lambda_max; ++lambda)
            integral += f->Query(lambda) * g->Query(lambda);
        return integral;
    }

    inline float InnerProduct(Spectrum* f, const Spectrum* g)
    {
        float integral = 0;
        for (float lambda = Lambda_min; lambda <= Lambda_max; ++lambda)
            integral += f->Query(lambda) * g->Query(lambda);
        return integral;
    }

    // Spectrum Inline Method Definitions
    /*
    inline float Spectrum::operator()(float lambda) const
    {
        auto op = [&](auto ptr) { return (*ptr)(lambda); };
        return Dispatch(op);
    }

    inline SampledSpectrum Spectrum::Sample(const SampledWavelengths& lambda) const
    {
        auto samp = [&](auto ptr) { return ptr->Sample(lambda); };
        return Dispatch(samp);
    }

    inline float Spectrum::MaxValue() const
    {
        auto max = [&](auto ptr) { return ptr->MaxValue(); };
        return Dispatch(max);
    }
    */

}  // namespace pbrt

namespace std {

    template <>
    struct hash<pbrt::DenselySampledSpectrum>
    {
        size_t operator()(const pbrt::DenselySampledSpectrum& s) const
        {
            return pbrt::HashBuffer(s.values.data(), s.values.size());
        }
    };

}  // namespace std