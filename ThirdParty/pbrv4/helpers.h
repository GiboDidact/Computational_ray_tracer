// pbrt is Copyright(c) 1998-2020 Matt Pharr, Wenzel Jakob, and Greg Humphreys.
// The pbrt source code is licensed under the Apache License, Version 2.0.
// SPDX: Apache-2.0
//EDITED
#pragma once
#include <cmath>

namespace pbrt {

    // CompensatedFloat Definition
    struct CompensatedFloat {
    public:
        // CompensatedFloat Public Methods
        CompensatedFloat(float v, float err = 0) : v(v), err(err) {}
        
        explicit operator float() const { return v + err; }
        
        explicit operator double() const { return double(v) + double(err); }
        std::string ToString() const;

        float v, err;
    };

    //fused multiply add
    inline  float FMA(float a, float b, float c)
    {
        return std::fma(a, b, c);
    }
    inline  double FMA(double a, double b, double c)
    {
        return std::fma(a, b, c);
    }
    inline  long double FMA(long double a, long double b, long double c)
    {
        return std::fma(a, b, c);
    }
    
    template <typename T>
    inline  typename std::enable_if_t<std::is_integral_v<T>, T> FMA(T a, T b,
        T c) {
        return a * b + c;
    }
    //FMA for Intervals

     inline uint32_t FloatToBits(float f)
     {
        return std::bit_cast<uint32_t>(f);
     }

     static constexpr float MachineEpsilon = std::numeric_limits<float>::epsilon() * 0.5;

    inline constexpr float gamma(int n) {
         return (n * MachineEpsilon) / (1 - n * MachineEpsilon);
     }

    template <typename Ta, typename Tb, typename Tc, typename Td>
    inline auto DifferenceOfProducts(Ta a, Tb b, Tc c, Td d) {
        auto cd = c * d;
        auto differenceOfProducts = FMA(a, b, -cd);
        auto error = FMA(-c, d, cd);
        return differenceOfProducts + error;
    }

    inline int MaxComponentIndex(glm::vec3 t) {
        return (t.x > t.y) ? ((t.x > t.z) ? 0 : 2) : ((t.y > t.z) ? 1 : 2);
    }

    inline float MaxComponentValue(glm::vec3 t) {
        using std::max;
        return max({ t.x, t.y, t.z });
    }


    inline  float ErfInv(float a)
    {
        // https://stackoverflow.com/a/49743348
        float p;
        float t = std::log(std::max(FMA(a, -a, 1.0f), std::numeric_limits<float>::min()));
        //CHECK(!IsNaN(t) && !IsInf(t));
        if (std::abs(t) > 6.125f) {          // maximum ulp error = 2.35793
            p = 3.03697567e-10f;             //  0x1.4deb44p-32
            p = FMA(p, t, 2.93243101e-8f);   //  0x1.f7c9aep-26
            p = FMA(p, t, 1.22150334e-6f);   //  0x1.47e512p-20
            p = FMA(p, t, 2.84108955e-5f);   //  0x1.dca7dep-16
            p = FMA(p, t, 3.93552968e-4f);   //  0x1.9cab92p-12
            p = FMA(p, t, 3.02698812e-3f);   //  0x1.8cc0dep-9
            p = FMA(p, t, 4.83185798e-3f);   //  0x1.3ca920p-8
            p = FMA(p, t, -2.64646143e-1f);  // -0x1.0eff66p-2
            p = FMA(p, t, 8.40016484e-1f);   //  0x1.ae16a4p-1
        }
        else {                             // maximum ulp error = 2.35456
            p = 5.43877832e-9f;              //  0x1.75c000p-28
            p = FMA(p, t, 1.43286059e-7f);   //  0x1.33b458p-23
            p = FMA(p, t, 1.22775396e-6f);   //  0x1.49929cp-20
            p = FMA(p, t, 1.12962631e-7f);   //  0x1.e52bbap-24
            p = FMA(p, t, -5.61531961e-5f);  // -0x1.d70c12p-15
            p = FMA(p, t, -1.47697705e-4f);  // -0x1.35be9ap-13
            p = FMA(p, t, 2.31468701e-3f);   //  0x1.2f6402p-9
            p = FMA(p, t, 1.15392562e-2f);   //  0x1.7a1e4cp-7
            p = FMA(p, t, -2.32015476e-1f);  // -0x1.db2aeep-3
            p = FMA(p, t, 8.86226892e-1f);   //  0x1.c5bf88p-1
        }
        return a * p;
    }

    template <typename T, typename U, typename V>
    inline  constexpr T Clamp(T val, U low, V high)
    {
        if (val < low)
            return T(low);
        else if (val > high)
            return T(high);
        else
            return val;
    }
    
    template <typename Float, typename C>
    inline  constexpr Float EvaluatePolynomial(Float t, C c) {
        return c;
    }
    
    
    template <typename Float, typename C, typename... Args>
    inline  constexpr Float EvaluatePolynomial(Float t, C c, Args... cRemaining) {
        return FMA(t, EvaluatePolynomial(t, cRemaining...), c);
    }
    
    /*
    // https://stackoverflow.com/a/10792321
    inline float FastExp(float x)
    {
        // Compute $x'$ such that $\roman{e}^x = 2^{x'}$
        float xp = x * 1.442695041f;

        // Find integer and fractional components of $x'$
        float fxp = std::floor(xp), f = xp - fxp;
        int i = (int)fxp;

        // Evaluate polynomial approximation of $2^f$
        float twoToF = EvaluatePolynomial(f, 1.f, 0.695556856f, 0.226173572f, 0.0781455737f);

        // Scale $2^f$ by $2^i$ and return final result
        int exponent = std::exp(twoToF) + i;
        if (exponent < -126)
            return 0;
        if (exponent > 127)
            return std::numeric_limits<float>::max();
        uint32_t bits = FloatToBits(twoToF);
        bits &= 0b10000000011111111111111111111111u;
        bits |= (exponent + 127) << 23;
        return std::bit_cast<float>(bits);
    }
    */
    inline float Lerp(float x, float a, float b)
    {
        return (1 - x) * a + x * b;
    }
    
    template <typename Predicate>
    inline  size_t FindInterval(size_t sz, const Predicate& pred)
    {
        using ssize_t = std::make_signed_t<size_t>;
        ssize_t size = (ssize_t)sz - 2, first = 1;
        while (size > 0) {
            // Evaluate predicate at midpoint and update _first_ and _size_
            size_t half = (size_t)size >> 1, middle = first + half;
            bool predResult = pred(middle);
            first = predResult ? middle + 1 : first;
            size = predResult ? size - (half + 1) : half;
        }
        return (size_t)Clamp((ssize_t)first - 1, 0, sz - 2);
    }
    
    inline  float SafeSqrt(float x)
    {
        //DCHECK_GE(x, -1e-3f);  // not too negative
        return std::sqrt(std::max(0.f, x));
    }

    inline CompensatedFloat TwoProd(float a, float b)
    {
        float ab = a * b;
        return { ab, FMA(a, b, -ab) };
    }

    inline CompensatedFloat TwoSum(float a, float b)
    {
        float s = a + b, delta = s - a;
        return { s, (a - (s - delta)) + (b - delta) };
    }

    namespace internal {
        // InnerProduct Helper Functions
        template <typename float>
        inline CompensatedFloat InnerProduct(float a, float b) {
            return TwoProd(a, b);
        }

        // Accurate dot products with FMA: Graillat et al.,
        // https://www-pequan.lip6.fr/~graillat/papers/posterRNC7.pdf
        //
        // Accurate summation, dot product and polynomial evaluation in complex
        // floating point arithmetic, Graillat and Menissier-Morain.
        template <typename float, typename... T>
        inline CompensatedFloat InnerProduct(float a, float b, T... terms) {
            CompensatedFloat ab = TwoProd(a, b);
            CompensatedFloat tp = InnerProduct(terms...);
            CompensatedFloat sum = TwoSum(ab.v, tp.v);
            return { sum.v, ab.err + (tp.err + sum.err) };
        }

    }  // namespace internal

    template <typename... T>
    inline std::enable_if_t<std::conjunction_v<std::is_arithmetic<T>...>, float>
        InnerProduct(T... terms) {
        CompensatedFloat ip = internal::InnerProduct(terms...);
        return Float(ip);
    }

    inline float Gaussian(float x, float mu = 0, float sigma = 1) 
    {
        return 1.0f / std::sqrt(2 * Pi * sigma * sigma) *
            std::exp(-std::powf(x - mu,2) / (2 * sigma * sigma));
    }

    inline float GaussianIntegral(float x0, float x1, float mu = 0, float sigma = 1) 
    {
        //DCHECK_GT(sigma, 0);
        float sigmaRoot2 = sigma * float(1.414213562373095);
        return 0.5f * (std::erf((mu - x0) / sigmaRoot2) - std::erf((mu - x1) / sigmaRoot2));
    }

    // http://www.plunk.org/~hatch/rightway.html
    inline float SinXOverX(float x)
    {
        if (1 - x * x == 1)
            return 1;
        return std::sin(x) / x;
    }

    inline float Sinc(float x)
    {
        return SinXOverX(Pi * x);
    }

    inline float WindowedSinc(float x, float radius, float tau) 
    {
        if (std::abs(x) > radius)
            return 0;
        return Sinc(x) * Sinc(x / tau);
    }

    //template <int N>
    //std::optional<glm::mat3> LinearLeastSquares(const float A[][N], const float B[][N], int rows);

    template <int N>
    std::optional<glm::mat3> LinearLeastSquares(const float A[][N], const float B[][N], int rows)
    {
        glm::mat3 AtA = glm::mat3(0.0f);
        glm::mat3 AtB = glm::mat3(0.0f);
        
        for (int i = 0; i < N; ++i)
            for (int j = 0; j < N; ++j)
                for (int r = 0; r < rows; ++r) {
                    AtA[i][j] += A[r][i] * A[r][j];
                    AtB[i][j] += A[r][i] * B[r][j];
                }
       
        glm::mat3 AtAi = glm::inverse(AtA);
        //if (!AtAi)
          //  return {};
        return { glm::transpose(AtAi * AtB) };
    }

}  // namespace pbrt