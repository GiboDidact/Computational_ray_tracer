// pbrt is Copyright(c) 1998-2020 Matt Pharr, Wenzel Jakob, and Greg Humphreys.
// The pbrt source code is licensed under the Apache License, Version 2.0.
// SPDX: Apache-2.0
//EDITED

#pragma once
#ifndef PBRT_FILTERS_H
#define PBRT_FILTERS_H

#include "../../pch.h"
#include "pbrt.h"
#include "../../RayTracer/Sampling.h"

namespace pbrt {

    // FilterSample Definition
    struct FilterSample 
    {
        glm::vec2 p;
        float weight; //f[pi] / pdf, f is the function of the filter, pi is domain value, pdf(pi)
    };

    class Filter 
    {
    public:
        // Filter Interface
        // static Filter Create(const std::string& name, const ParameterDictionary& parameters, const FileLoc* loc);

        virtual inline glm::vec2 Radius() const = 0;

        virtual inline float Evaluate(glm::vec2 p) const = 0;

        virtual inline float Integral() const = 0;

        virtual inline FilterSample Sample(glm::vec2 u) const = 0;

        virtual std::string ToString() const = 0;
    };

    /*
    //helps with sampling filters with difficult functions, I can just use my reject or inverse method
    class FilterSampler 
    {
    public:
        // FilterSampler Public Methods
        FilterSampler(Filter filter);
        std::string ToString() const;

        FilterSample Sample(glm::vec2 u) const 
        {
            float pdf;
            glm::ivec2 pi;
            glm::vec2 p = distrib.Sample(u, &pdf, &pi);
            return FilterSample{ p, f[pi] / pdf };
        }

    private:
        // FilterSampler Private Members
        Bounds2f domain;
        Array2D<float> f;
        PiecewiseConstant2D distrib;
    };
    */

    // BoxFilter Definition
    class BoxFilter : public Filter
    {
    public:
        // BoxFilter Public Methods
        BoxFilter(glm::vec2 _radius = glm::vec2(0.5, 0.5)) : radius(_radius) {}

        //static BoxFilter* Create(const ParameterDictionary& parameters, const FileLoc* loc);

        glm::vec2 Radius() const override { return radius; }

        std::string ToString() const override;

        float Evaluate(glm::vec2 p) const override
        {
            return (std::abs(p.x) <= radius.x && std::abs(p.y) <= radius.y) ? 1 : 0;
        }

        FilterSample Sample(glm::vec2 u) const override
        {
            glm::vec2 p(Lerp(u[0], -radius.x, radius.x), Lerp(u[1], -radius.y, radius.y));
            return { p, float(1) };
        }

        float Integral() const override { return 2 * radius.x * 2 * radius.y; }

    private:
        glm::vec2 radius;
    };

    // GaussianFilter Definition
    class GaussianFilter : public Filter
    {
    public:
        // GaussianFilter Public Methods
        GaussianFilter(glm::vec2 radius, float sigma = 0.5f)
            : radius(radius),
            sigma(sigma),
            expX(Gaussian(radius.x, 0, sigma)),
            expY(Gaussian(radius.y, 0, sigma)),
            inv_sampler_x([=](float x) -> float { return this->EvaluateX(x); }, -radius.x, radius.x, 10000),
            inv_sampler_y([=](float y) -> float { return this->EvaluateY(y); }, -radius.y, radius.y, 10000) {}

        //static GaussianFilter* Create(const ParameterDictionary& parameters, const FileLoc* loc);

        glm::vec2 Radius() const override { return radius; }

        std::string ToString() const override;

        float Evaluate(glm::vec2 p) const override
        {
            return (std::max<float>(0, Gaussian(p.x, 0, sigma) - expX) *
                std::max<float>(0, Gaussian(p.y, 0, sigma) - expY));
        }

        
        float Integral() const override
        {
            return ((GaussianIntegral(-radius.x, radius.x, 0, sigma) - 2 * radius.x * expX) *
                (GaussianIntegral(-radius.y, radius.y, 0, sigma) - 2 * radius.y * expY));
        }

        
        //the filtering integral sampling uses f(x,y)/p(x,y)
        FilterSample Sample(glm::vec2 u) const override 
        { 
            FilterSample fs;
            fs.p.x = inv_sampler_x.Sample(u.x, true);
            fs.p.y = inv_sampler_y.Sample(u.y, true);
            fs.weight = Evaluate(fs.p) / (inv_sampler_x.PDF(fs.p.x) * inv_sampler_y.PDF(fs.p.y));//or use 1
            return fs;
            
            /*
            FilterSample fs;
            fs.p = Continuous_Rejection_Sample2D([=](glm::vec2 pos) -> float { return this->Evaluate(pos); }, -radius.x, radius.x, -radius.y, radius.y, 0);
            fs.weight = Evaluate(fs.p) / (inv_sampler_x.PDF(fs.p.x) * inv_sampler_y.PDF(fs.p.y));//or use 1
            return fs;
            */
        }

    private:
        // GaussianFilter Private Members
        glm::vec2 radius;
        float sigma, expX, expY;
        //FilterSampler sampler;
        Continuous_Inversion_Sampler inv_sampler_x;
        Continuous_Inversion_Sampler inv_sampler_y;

        float EvaluateX(float x) const
        {
            return std::max<float>(0, Gaussian(x, 0, sigma) - expX);
        }
        float EvaluateY(float y) const
        {
            return std::max<float>(0, Gaussian(y, 0, sigma) - expY);
        }
    };

    /*
    // MitchellFilter Definition
    class MitchellFilter : public Filter
    {
    public:
        // MitchellFilter Public Methods
        MitchellFilter(glm::vec2 radius, float b = 1.f / 3.f, float c = 1.f / 3.f)
            : radius(radius), b(b), c(c), sampler(this) {}

        //static MitchellFilter* Create(const ParameterDictionary& parameters, const FileLoc* loc);

        
        glm::vec2 Radius() const override { return radius; }

        std::string ToString() const override;

        
        float Evaluate(glm::vec2 p) const  override
        {
            return Mitchell1D(2 * p.x / radius.x) * Mitchell1D(2 * p.y / radius.y);
        }

        
        FilterSample Sample(glm::vec2 u) const override  { return sampler.Sample(u); }

        
        float Integral() const override { return radius.x * radius.y / 4; }

    private:
        // MitchellFilter Private Methods
        float Mitchell1D(float x) const 
        {
            x = std::abs(x);
            if (x <= 1)
                return ((12 - 9 * b - 6 * c) * x * x * x + (-18 + 12 * b + 6 * c) * x * x +
                    (6 - 2 * b)) *
                (1.f / 6.f);
            else if (x <= 2)
                return ((-b - 6 * c) * x * x * x + (6 * b + 30 * c) * x * x +
                    (-12 * b - 48 * c) * x + (8 * b + 24 * c)) *
                (1.f / 6.f);
            else
                return 0;
        }

        // MitchellFilter Private Members
        glm::vec2 radius;
        float b, c;
        //FilterSampler sampler;
    };
    */

    // LanczosSincFilter Definition
    class LanczosSincFilter : public Filter
    {
    public:
        // LanczosSincFilter Public Methods
        LanczosSincFilter(glm::vec2 radius, float tau = 3.f)
            : radius(radius), tau(tau), 
            inv_sampler_x([=](float x) -> float { return this->EvaluateX(x); }, -radius.x, radius.x, 2000),
            inv_sampler_y([=](float y) -> float { return this->EvaluateY(y); }, -radius.y, radius.y, 2000) {}

        //static LanczosSincFilter* Create(const ParameterDictionary& parameters, const FileLoc* loc);

        
        glm::vec2 Radius() const override { return radius; }

        std::string ToString() const override;

        float Evaluate(glm::vec2 p) const override
        {
            return WindowedSinc(p.x, radius.x, tau) * WindowedSinc(p.y, radius.y, tau);
        }

        FilterSample Sample(glm::vec2 u) const override 
        { 
            FilterSample fs;
            fs.p.x = inv_sampler_x.Sample();
            fs.p.y = inv_sampler_y.Sample();
            fs.weight = Evaluate(fs.p) / (inv_sampler_x.PDF(fs.p.x) * inv_sampler_y.PDF(fs.p.y));//or use 1
            return fs;
        }

        float Integral() const override;

    private:
        // LanczosSincFilter Private Members
        glm::vec2 radius;
        float tau;
        //FilterSampler sampler;
        Continuous_Inversion_Sampler inv_sampler_x;
        Continuous_Inversion_Sampler inv_sampler_y;

        float EvaluateX(float x) const
        {
            return WindowedSinc(x, radius.x, tau);
        }
        float EvaluateY(float y) const
        {
            return  WindowedSinc(y, radius.y, tau);
        }
    };

    // TriangleFilter Definition
    class TriangleFilter : public Filter
    {
    public:
        // TriangleFilter Public Methods
        TriangleFilter(glm::vec2 radius) : radius(radius) {}

       // static TriangleFilter* Create(const ParameterDictionary& parameters, const FileLoc* loc);
   
        glm::vec2 Radius() const override { return radius; }

        std::string ToString() const override;
        
        float Evaluate(glm::vec2 p) const override
        {
            return std::max<float>(0, radius.x - std::abs(p.x)) *
                std::max<float>(0, radius.y - std::abs(p.y));
        }
        
        FilterSample Sample(glm::vec2 u) const override
        {
            //internal u
            return { glm::vec2(SampleTent(u[0], radius.x), SampleTent(u[1], radius.y)),
                    float(1) };
        }

        float Integral() const override { return std::pow(radius.x,2) * std::pow(radius.y,2); }

    private:
        glm::vec2 radius;
    };

    /*
    inline float Filter::Evaluate(glm::vec2 p) const 
    {
        auto eval = [&](auto ptr) { return ptr->Evaluate(p); };
        return Dispatch(eval);
    }

    inline FilterSample Filter::Sample(glm::vec2 u) const 
    {
        auto sample = [&](auto ptr) { return ptr->Sample(u); };
        return Dispatch(sample);
    }

    inline glm::vec2 Filter::Radius() const
    {
        auto radius = [&](auto ptr) { return ptr->Radius(); };
        return Dispatch(radius);
    }

    inline float Filter::Integral() const 
    {
        auto integral = [&](auto ptr) { return ptr->Integral(); };
        return Dispatch(integral);
    }
    */

}  // namespace pbrt

#endif  // PBRT_FILTERS_H
