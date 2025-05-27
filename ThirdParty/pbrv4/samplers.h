// pbrt is Copyright(c) 1998-2020 Matt Pharr, Wenzel Jakob, and Greg Humphreys.
// The pbrt source code is licensed under the Apache License, Version 2.0.
// SPDX: Apache-2.0
//EDITED

#pragma once
#include "../../pch.h"
#include "rng.h"
#include "hash.h"

/*
If you look at the unbiased natural uniform monte sampler it seems like you can't improve variance, the variance
is simply 1/Nvar(g(x)), however if you use fourier analysis it turns out you can reduce variance on the method of generating
the uniform! Because we generate uniform on computer, different methods actually have different variances

. So these classes generate uniform variables of d dimension, but the way they do it is different as to reduce variance

. they also have an option to make the sampling deterministic, that way we can reproduce the same image, each pixel/sample gets the same variables

(directly based off pbr4 sampler class and code)
*/

namespace pbrt {

    class Sampler
    {
    public:
        virtual int SamplesPerPixel() const = 0;
        virtual void StartPixelSample(glm::ivec2 p, int sampleIndex, int dimension = 0) = 0;

        virtual float Get1D() = 0;
        virtual glm::vec2 Get2D() = 0;
        virtual glm::vec2 GetPixel2D() = 0;

        //virtual std::string ToString() const = 0;
    };

    class IndependentSampler : public Sampler
    {
    public:
        IndependentSampler(int _samplesPerPixel, int _seed = 0) : samplesPerPixel(_samplesPerPixel), seed(_seed) {}

        static constexpr const char* Name() { return "IndependentSampler"; }

        int SamplesPerPixel() const override { return samplesPerPixel; }

        void StartPixelSample(glm::ivec2 p, int sampleIndex, int dimension) override
        {
            rng.SetSequence(Hash(p, seed));
            rng.Advance(sampleIndex * 65536ull + dimension);
        }

        float Get1D() override { return rng.Uniform<float>(); }
        glm::vec2 Get2D() override { return { rng.Uniform<float>(), rng.Uniform<float>() }; }
        glm::vec2 GetPixel2D() override { return Get2D(); }

        //std::string ToString() const override;

    private:
        int samplesPerPixel, seed;
        RNG rng;
    };


    //the x/ypixel samples is the grid size created, its the x by y grid that is creates the spaced out samples from
    class StratifiedSampler : public Sampler {
    public:
        // StratifiedSampler Public Methods
        StratifiedSampler(int xPixelSamples, int yPixelSamples, bool jitter, int seed = 0)
            : xPixelSamples(xPixelSamples),
            yPixelSamples(yPixelSamples),
            seed(seed),
            jitter(jitter) {}

        static constexpr const char* Name() { return "StratifiedSampler"; }


        int SamplesPerPixel() const override { return xPixelSamples * yPixelSamples; }

        void StartPixelSample(glm::ivec2 p, int index, int dim) override
        {
            //if you have more indices than the grid hash will loop forever
            if (jitter == false && index >= SamplesPerPixel())
            {
                std::cout << "more sampels than pixels for strat\n";
                return;
            }
            pixel = p;
            sampleIndex = index;
            dimension = dim;
            rng.SetSequence(Hash(p, seed));
            rng.Advance(sampleIndex * 65536ull + dimension);
        }

        float Get1D() override
        {
            // Compute _stratum_ index for current pixel and dimension
            uint64_t hash = Hash(pixel, dimension, seed);
            int stratum = Helper::PermutationElement(sampleIndex, SamplesPerPixel(), hash);

            ++dimension;
            float delta = jitter ? rng.Uniform<float>() : 0.5f;
            return (stratum + delta) / SamplesPerPixel();
        }


        glm::vec2 Get2D() override
        {
            if (sampleIndex >= SamplesPerPixel())
            {
                return glm::vec2(0, 0);
            }
            // Compute _stratum_ index for current pixel and dimension
            uint64_t hash = Hash(pixel, dimension, seed);
            int stratum = Helper::PermutationElement(sampleIndex, SamplesPerPixel(), hash);

            dimension += 2;
            int x = stratum % xPixelSamples, y = stratum / xPixelSamples;
            float dx = jitter ? rng.Uniform<float>() : 0.5f;
            float dy = jitter ? rng.Uniform<float>() : 0.5f;

            return { (x + dx) / xPixelSamples, (y + dy) / yPixelSamples };
        }

        glm::vec2 GetPixel2D() override { return Get2D(); }

        //std::string ToString() const override;

    private:
        // StratifiedSampler Private Members
        int xPixelSamples, yPixelSamples, seed;
        bool jitter;
        RNG rng;
        glm::ivec2 pixel;
        int sampleIndex = 0, dimension = 0;
    };







    //SOBOL
    static constexpr float FloatOneMinusEpsilon = 0x1.fffffep-1;


    struct NoRandomizer {
        uint32_t operator()(uint32_t v) const { return v; }
    };

    // BinaryPermuteScrambler Definition
    struct BinaryPermuteScrambler {

        BinaryPermuteScrambler(uint32_t perm) : permutation(perm) {}

        uint32_t operator()(uint32_t v) const { return permutation ^ v; }
        uint32_t permutation;
    };

    // FastOwenScrambler Definition
    struct FastOwenScrambler {
        FastOwenScrambler(uint32_t seed) : seed(seed) {}
        // FastOwenScrambler Public Methods

        uint32_t operator()(uint32_t v) const {
            v = Helper::ReverseBits32(v);
            v ^= v * 0x3d20adea;
            v += seed;
            v *= (seed >> 16) | 1;
            v ^= v * 0x05526c56;
            v ^= v * 0x53a22864;
            return Helper::ReverseBits32(v);
        }

        uint32_t seed;
    };

    // OwenScrambler Definition
    struct OwenScrambler {
        OwenScrambler(uint32_t seed) : seed(seed) {}
        // OwenScrambler Public Methods
        uint32_t operator()(uint32_t v) const {
            if (seed & 1)
                v ^= 1u << 31;
            for (int b = 1; b < 32; ++b) {
                // Apply Owen scrambling to binary digit _b_ in _v_
                uint32_t mask = (~0u) << (32 - b);
                if ((uint32_t)Helper::MixBits((v & mask) ^ seed) & (1u << b))
                    v ^= 1u << (31 - b);
            }
            return v;
        }

        uint32_t seed;
    };

    enum class RandomizeStrategy { None, PermuteDigits, FastOwen, Owen };

    template <typename R>
    inline float SobolSample(int64_t a, int dimension, R randomizer)
    {
        //DCHECK_LT(dimension, NSobolDimensions);
       // DCHECK(a >= 0 && a < (1ull << SobolMatrixSize));
        // Compute initial Sobol\+$'$ sample _v_ using generator matrices
        uint32_t v = 0;
        for (int i = dimension * Helper::SobolMatrixSize; a != 0; a >>= 1, i++)
            if (a & 1)
                v ^= Helper::SobolMatrices32[i];

        // Randomize Sobol\+$'$ sample and return floating-point value
        v = randomizer(v);
        return std::min(v * 0x1p-32f, FloatOneMinusEpsilon);
    }

    inline uint64_t SobolIntervalToIndex(uint32_t m, uint64_t frame, glm::ivec2 p)
    {
        if (m == 0)
            return frame;

        const uint32_t m2 = m << 1;
        uint64_t index = uint64_t(frame) << m2;

        uint64_t delta = 0;
        for (int c = 0; frame; frame >>= 1, ++c)
            if (frame & 1)  // Add flipped column m + c + 1.
                delta ^= Helper::VdCSobolMatrices[m - 1][c];

        // flipped b
        uint64_t b = (((uint64_t)((uint32_t)p.x) << m) | ((uint32_t)p.y)) ^ delta;

        for (int c = 0; b; b >>= 1, ++c)
            if (b & 1)  // Add column 2 * m - c.
                index ^= Helper::VdCSobolMatricesInv[m - 1][c];

        return index;
    }

    class SobolSampler : public Sampler
    {
    public:
        // SobolSampler Public Methods
        SobolSampler(int samplesPerPixel, glm::ivec2 fullResolution, RandomizeStrategy randomize, int seed = 0)
            : samplesPerPixel(samplesPerPixel), seed(seed), randomize(randomize)
        {
            // if (!IsPowerOf2(samplesPerPixel))
               //  Warning("Non power-of-two sample count %d will perform suboptimally with the "
                  //   "SobolSampler.",
                   //  samplesPerPixel);

            scale = Helper::RoundUpPow2(std::max(fullResolution.x, fullResolution.y));
        }


        static constexpr const char* Name() { return "SobolSampler"; }


        int SamplesPerPixel() const override { return samplesPerPixel; }


        void StartPixelSample(glm::ivec2 p, int sampleIndex, int dim) override
        {
            pixel = p;
            dimension = std::max(2, dim);
            sobolIndex = SobolIntervalToIndex(std::log2(scale), sampleIndex, pixel); //Log2()
        }


        float Get1D() override
        {
            if (dimension >= Helper::NSobolDimensions)
                dimension = 2;
            return SampleDimension(dimension++);
        }


        glm::vec2 Get2D() override
        {
            if (dimension + 1 >= Helper::NSobolDimensions)
                dimension = 2;
            glm::vec2 u(SampleDimension(dimension), SampleDimension(dimension + 1));
            dimension += 2;
            return u;
        }

        glm::vec2 GetPixel2D() override
        {
            glm::vec2 u(SobolSample(sobolIndex, 0, NoRandomizer()),
                SobolSample(sobolIndex, 1, NoRandomizer()));
            // Remap Sobol\+$'$ dimensions used for pixel samples
            for (int dim = 0; dim < 2; ++dim)
            {
                // DCHECK_RARE(1e-7, u[dim] * scale - pixel[dim] < 0);
                 //DCHECK_RARE(1e-7, u[dim] * scale - pixel[dim] > 1);
                u[dim] = glm::clamp(u[dim] * scale - pixel[dim], 0.0f, OneMinusEpsilon);
            }

            return u;
        }

        //std::string ToString() const override;

    private:
        // SobolSampler Private Methods
        float SampleDimension(int dimension) const
        {
            // Return un-randomized Sobol\+$'$ sample if appropriate
            if (randomize == RandomizeStrategy::None)
                return SobolSample(sobolIndex, dimension, NoRandomizer());

            // Return randomized Sobol\+$'$ sample using _randomize_
            uint32_t hash = Hash(dimension, seed);
            if (randomize == RandomizeStrategy::PermuteDigits)
                return SobolSample(sobolIndex, dimension, BinaryPermuteScrambler(hash));
            else if (randomize == RandomizeStrategy::FastOwen)
                return SobolSample(sobolIndex, dimension, FastOwenScrambler(hash));
            else
                return SobolSample(sobolIndex, dimension, OwenScrambler(hash));
        }

        // SobolSampler Private Members
        int samplesPerPixel, scale, seed;
        RandomizeStrategy randomize;
        glm::vec2 pixel;
        int dimension;
        int64_t sobolIndex;
    };

}  // namespace pbrt