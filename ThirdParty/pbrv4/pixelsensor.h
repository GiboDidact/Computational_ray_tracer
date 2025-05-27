// pbrt is Copyright(c) 1998-2020 Matt Pharr, Wenzel Jakob, and Greg Humphreys.
// The pbrt source code is licensed under the Apache License, Version 2.0.
// SPDX: Apache-2.0
//EDITED

#pragma once
#ifndef PBRT_FILM_H
#define PBRT_FILM_H

// PhysLight code contributed by Anders Langlands and Luca Fascione
// Copyright (c) 2020, Weta Digital, Ltd.
// SPDX-License-Identifier: Apache-2.0

#include "pbrt.h"
#include "colorspace.h"
#include "spectrum.h"


namespace pbrt {
    /*
    exposure - just a scalar you multiply the rgb by

    rgb response - constructs the spectrumrgb and rgbtoxyz transformations with linear algebra

    white balancing - one algorithm to convert to lms and scale, other one implicitly scales it. just the ratio of green channel wtih illimant or something  
    */
    // PixelSensor Definition
    class PixelSensor {
    public:
        // PixelSensor Public Methods
       // static PixelSensor* Create(const ParameterDictionary& parameters,
           // const RGBColorSpace* colorSpace, Float exposureTime,
           // const FileLoc* loc, Allocator alloc);

        //static PixelSensor* CreateDefault();

        PixelSensor(Spectrum* r, Spectrum* g, Spectrum* b, const RGBColorSpace* outputColorSpace, Spectrum* sensorIllum, float imagingRatio)
            : r_bar(r), g_bar(g), b_bar(b), imagingRatio(imagingRatio) 
        {
            // Compute XYZ from camera RGB matrix
            // Compute _rgbCamera_ values for training swatches
            float rgbCamera[nSwatchReflectances][3];
            for (int i = 0; i < nSwatchReflectances; ++i) {
                RGB rgb = ProjectReflectance<RGB>(swatchReflectances[i], sensorIllum, &r_bar,
                    &g_bar, &b_bar);
                for (int c = 0; c < 3; ++c)
                    rgbCamera[i][c] = rgb[c];
            }

            // Compute _xyzOutput_ values for training swatches
            float xyzOutput[24][3];
            float sensorWhiteG = InnerProduct(sensorIllum, &g_bar);
            float sensorWhiteY = InnerProduct(sensorIllum, &Spectra::Y());
            for (size_t i = 0; i < nSwatchReflectances; ++i) 
            {
                Spectrum* s = swatchReflectances[i];
                XYZ xyz = ProjectReflectance<XYZ>(s, (Spectrum*)&outputColorSpace->illuminant, 
                             (Spectrum*)&Spectra::X(), (Spectrum*)&Spectra::Y(), (Spectrum*)&Spectra::Z()) * (sensorWhiteY / sensorWhiteG);
                for (int c = 0; c < 3; ++c)
                    xyzOutput[i][c] = xyz[c];
            }

            // Initialize _XYZFromSensorRGB_ using linear least squares
            std::optional<glm::mat3> m = LinearLeastSquares(rgbCamera, xyzOutput, nSwatchReflectances);
            if (!m)
                std::cout << "Sensor XYZ from RGB matrix could not be solved." << std::endl;
            XYZFromSensorRGB = *m;
        }

        PixelSensor(const RGBColorSpace* outputColorSpace, Spectrum* sensorIllum, float imagingRatio)
            : r_bar(Spectra::X()), g_bar(Spectra::Y()), b_bar(Spectra::Z()), imagingRatio(imagingRatio) 
        {
            // Compute white balancing matrix for XYZ _PixelSensor_
            if (sensorIllum) {
                glm::vec2 sourceWhite = SpectrumToXYZ(sensorIllum).xy();
                glm::vec2 targetWhite = outputColorSpace->w;
                XYZFromSensorRGB = WhiteBalance(sourceWhite, targetWhite);
            }
        }

        RGB ToSensorRGB(SampledSpectrum L, const SampledWavelengths& lambda) const 
        {
            L = SafeDiv(L, lambda.PDF());
            return imagingRatio * RGB((r_bar.Sample(lambda) * L).Average(),
                (g_bar.Sample(lambda) * L).Average(),
                (b_bar.Sample(lambda) * L).Average());
        }

        // PixelSensor Public Members
        glm::mat3 XYZFromSensorRGB;

    private:
        // PixelSensor Private Methods
        template <typename Triplet>
        static Triplet ProjectReflectance(Spectrum* r, Spectrum* illum, Spectrum* b1, Spectrum* b2, Spectrum* b3);

        // PixelSensor Private Members
        DenselySampledSpectrum r_bar, g_bar, b_bar;
        float imagingRatio;
        static constexpr int nSwatchReflectances = 24;
        static Spectrum* swatchReflectances[nSwatchReflectances];
    };

    // PixelSensor Inline Methods
    template <typename Triplet>
    inline Triplet PixelSensor::ProjectReflectance(Spectrum* refl, Spectrum* illum, Spectrum* b1, Spectrum* b2, Spectrum* b3)
    {
        Triplet result;
        float g_integral = 0;
        for (float lambda = Lambda_min; lambda <= Lambda_max; ++lambda)
        {
            g_integral += b2->Query(lambda) * illum->Query(lambda);
            result[0] += b1->Query(lambda) * refl->Query(lambda) * illum->Query(lambda);
            result[1] += b2->Query(lambda) * refl->Query(lambda) * illum->Query(lambda);
            result[2] += b3->Query(lambda) * refl->Query(lambda) * illum->Query(lambda);
        }
        return result / g_integral;
    }


}  // namespace pbrt

#endif  // PBRT_FILM_H