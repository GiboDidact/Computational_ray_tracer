// pbrt is Copyright(c) 1998-2020 Matt Pharr, Wenzel Jakob, and Greg Humphreys.
// The pbrt source code is licensed under the Apache License, Version 2.0.
// SPDX: Apache-2.0
//EDITED
#pragma once

#include "../../pch.h"
#include "pbrt.h"
#include "spectrum.h"
#include "color.h"

namespace pbrt {
    

    //an rgb color space really just needs the transform from XYZ or to RGB, and to create it you give the 4 chromaticity values for r,g,b,w
    //technically you could give the 3 spectral functions for each color and integrate directly like how we did for XYZ
    // RGBColorSpace Definition
    class RGBColorSpace
    {
    public:
        // RGBColorSpace Public Methods
        RGBColorSpace(glm::vec2 r, glm::vec2 g, glm::vec2 b, Spectrum* illuminant, const RGBToSpectrumTable* rgbToSpectrumTable);
        
        RGBSigmoidPolynomial ToRGBCoeffs(RGB rgb) const;

        static void Init();

        // RGBColorSpace Public Members
        glm::vec2 r, g, b, w;
        DenselySampledSpectrum illuminant;
        glm::mat3 XYZFromRGB, RGBFromXYZ;
        static const RGBColorSpace* sRGB, * DCI_P3, * Rec2020, * ACES2065_1;


        bool operator==(const RGBColorSpace& cs) const
        {
            return (r == cs.r && g == cs.g && b == cs.b && w == cs.w &&
                rgbToSpectrumTable == cs.rgbToSpectrumTable);
        }

        bool operator!=(const RGBColorSpace& cs) const
        {
            return (r != cs.r || g != cs.g || b != cs.b || w != cs.w ||
                rgbToSpectrumTable != cs.rgbToSpectrumTable);
        }

        std::string ToString() const;


        RGB LuminanceVector() const
        {
            return RGB(XYZFromRGB[1][0], XYZFromRGB[1][1], XYZFromRGB[1][2]);
        }

        RGB ToRGB(XYZ xyz) const
        {
            glm::vec3 _rgb = RGBFromXYZ * xyz.getglm();
            return RGB(_rgb.x, _rgb.y, _rgb.z);
        }

        XYZ ToXYZ(RGB rgb) const
        {
            glm::vec3 _xyz = XYZFromRGB * rgb.getglm();
            return XYZ(_xyz.x, _xyz.y, _xyz.z);
        }

        static const RGBColorSpace* GetNamed(std::string name);
        static const RGBColorSpace* Lookup(glm::vec2 r, glm::vec2 g, glm::vec2 b, glm::vec2 w);

    private:
        // RGBColorSpace Private Members
        const RGBToSpectrumTable* rgbToSpectrumTable;
    };

    glm::mat3 ConvertRGBColorSpace(const RGBColorSpace& from, const RGBColorSpace& to);

}  // namespace pbrt

