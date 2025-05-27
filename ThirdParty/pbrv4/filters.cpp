// pbrt is Copyright(c) 1998-2020 Matt Pharr, Wenzel Jakob, and Greg Humphreys.
// The pbrt source code is licensed under the Apache License, Version 2.0.
// SPDX: Apache-2.0
//EDITED

#include "../../pch.h"
#include "filters.h"
#include "rng.h"

namespace pbrt {

    // Box Filter Method Definitions
    std::string BoxFilter::ToString() const 
    {
        return "[ BoxFilter radius: " + std::to_string(radius.x) + ", " + std::to_string(radius.y);
    }

    // Gaussian Filter Method Definitions
    std::string GaussianFilter::ToString() const 
    {
        return "[ GaussianFilter radius: " + std::to_string(radius.x) + ", " + std::to_string(radius.y) +
            " sigma: " + std::to_string(sigma) + " expX: " + std::to_string(expX) + " expY: " + std::to_string(expY);
    }

    // Sinc Filter Method Definitions
    float LanczosSincFilter::Integral() const 
    {
        float sum = 0;
        int sqrtSamples = 64;
        int nSamples = sqrtSamples * sqrtSamples;
        float area = 2 * radius.x * 2 * radius.y;
        RNG rng;
        for (int y = 0; y < sqrtSamples; ++y) {
            for (int x = 0; x < sqrtSamples; ++x) {
                glm::vec2 u((x + rng.Uniform<float>()) / sqrtSamples,
                    (y + rng.Uniform<float>()) / sqrtSamples);
                glm::vec2 p(Lerp(u.x, -radius.x, radius.x), Lerp(u.y, -radius.y, radius.y));
                sum += Evaluate(p);
            }
        }
        return sum / nSamples * area;
    }

    std::string LanczosSincFilter::ToString() const 
    {
        return "[ LanczosSincFilter radius: " + std::to_string(radius.x) + ", " + std::to_string(radius.y)
            + " tau: " + std::to_string(tau);
    }

    // Triangle Filter Method Definitions
    std::string TriangleFilter::ToString() const 
    {
        return "[ TriangleFilter radius: " + std::to_string(radius.x) + ", " + std::to_string(radius.y);
    }

    /*
        
    std::string Filter::ToString() const 
    {
        if (!ptr())
            return "(nullptr)";

        auto ts = [&](auto ptr) { return ptr->ToString(); };
        return DispatchCPU(ts);
    }

    BoxFilter* BoxFilter::Create(const ParameterDictionary& parameters, const FileLoc* loc,
        Allocator alloc) {
        Float xw = parameters.GetOneFloat("xradius", 0.5f);
        Float yw = parameters.GetOneFloat("yradius", 0.5f);
        return alloc.new_object<BoxFilter>(Vector2f(xw, yw));
    }
    
    GaussianFilter* GaussianFilter::Create(const ParameterDictionary& parameters,
        const FileLoc* loc, Allocator alloc) {
        // Find common filter parameters
        Float xw = parameters.GetOneFloat("xradius", 1.5f);
        Float yw = parameters.GetOneFloat("yradius", 1.5f);
        Float sigma = parameters.GetOneFloat("sigma", 0.5f);  // equivalent to old alpha = 2
        return alloc.new_object<GaussianFilter>(Vector2f(xw, yw), sigma, alloc);
    }
    
    // Mitchell Filter Method Definitions
    std::string MitchellFilter::ToString() const {
        return StringPrintf("[ MitchellFilter radius: %s b: %f c: %f sampler: %s ]", radius,
            b, c, sampler);
    }

    MitchellFilter* MitchellFilter::Create(const ParameterDictionary& parameters,
        const FileLoc* loc, Allocator alloc) {
        // Find common filter parameters
        Float xw = parameters.GetOneFloat("xradius", 2.f);
        Float yw = parameters.GetOneFloat("yradius", 2.f);
        Float B = parameters.GetOneFloat("B", 1.f / 3.f);
        Float C = parameters.GetOneFloat("C", 1.f / 3.f);
        return alloc.new_object<MitchellFilter>(Vector2f(xw, yw), B, C, alloc);
    }
    
    LanczosSincFilter* LanczosSincFilter::Create(const ParameterDictionary& parameters,
        const FileLoc* loc, Allocator alloc) {
        Float xw = parameters.GetOneFloat("xradius", 4.);
        Float yw = parameters.GetOneFloat("yradius", 4.);
        Float tau = parameters.GetOneFloat("tau", 3.f);
        return alloc.new_object<LanczosSincFilter>(Vector2f(xw, yw), tau, alloc);
    }

    TriangleFilter* TriangleFilter::Create(const ParameterDictionary& parameters,
        const FileLoc* loc, Allocator alloc) {
        // Find common filter parameters
        Float xw = parameters.GetOneFloat("xradius", 2.f);
        Float yw = parameters.GetOneFloat("yradius", 2.f);
        return alloc.new_object<TriangleFilter>(Vector2f(xw, yw));
    }

    Filter Filter::Create(const std::string& name, const ParameterDictionary& parameters,
        const FileLoc* loc, Allocator alloc) {
        Filter filter = nullptr;
        if (name == "box")
            filter = BoxFilter::Create(parameters, loc, alloc);
        else if (name == "gaussian")
            filter = GaussianFilter::Create(parameters, loc, alloc);
        else if (name == "mitchell")
            filter = MitchellFilter::Create(parameters, loc, alloc);
        else if (name == "sinc")
            filter = LanczosSincFilter::Create(parameters, loc, alloc);
        else if (name == "triangle")
            filter = TriangleFilter::Create(parameters, loc, alloc);
        else
            ErrorExit(loc, "%s: filter type unknown.", name);

        if (!filter)
            ErrorExit(loc, "%s: unable to create filter.", name);

        parameters.ReportUnused();
        return filter;
    }

    // FilterSampler Method Definitions
    FilterSampler::FilterSampler(Filter filter, Allocator alloc)
        : domain(Point2f(-filter.Radius()), Point2f(filter.Radius())),
        f(int(32 * filter.Radius().x), int(32 * filter.Radius().y), alloc),
        distrib(alloc) {
        // Tabularize unnormalized filter function in _f_
        for (int y = 0; y < f.YSize(); ++y)
            for (int x = 0; x < f.XSize(); ++x) {
                Point2f p =
                    domain.Lerp(Point2f((x + 0.5f) / f.XSize(), (y + 0.5f) / f.YSize()));
                f(x, y) = filter.Evaluate(p);
            }

        // Compute sampling distribution for filter
        distrib = PiecewiseConstant2D(f, domain, alloc);
    }

    std::string FilterSampler::ToString() const {
        return StringPrintf("[ FilterSampler domain: %s f: %s distrib: %s ]", domain, f,
            distrib);
    }
    */
}  // namespace pbrt