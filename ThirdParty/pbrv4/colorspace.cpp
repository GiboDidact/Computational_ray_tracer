// pbrt is Copyright(c) 1998-2020 Matt Pharr, Wenzel Jakob, and Greg Humphreys.
// The pbrt source code is licensed under the Apache License, Version 2.0.
// SPDX: Apache-2.0
//EDITED

#include "../../pch.h"
#include "colorspace.h"

namespace pbrt {

    // RGBColorSpace Method Definitions
    //pass in the chromaticity rgb values which represent the colornesss, then the white illuminant with all 3 on to calculuate xyY whitepoint
    RGBColorSpace::RGBColorSpace(glm::vec2 r, glm::vec2 g, glm::vec2 b, Spectrum* illuminant, const RGBToSpectrumTable* rgbToSpec)
                            : r(r), g(g), b(b), illuminant(illuminant), rgbToSpectrumTable(rgbToSpec)
    {
        // Compute whitepoint primaries and XYZ coordinates
        XYZ W = SpectrumToXYZ(illuminant);
        w = W.xy();
        XYZ R = XYZ::FromxyY(r), G = XYZ::FromxyY(g), B = XYZ::FromxyY(b);
        
        // Initialize XYZ color space conversion matrices
        glm::mat3 rgb(R.X, R.Y, R.Z, G.X, G.Y, G.Z, B.X, B.Y, B.Z);
        glm::vec3 _C = glm::inverse(rgb) * W.getglm(); //debug inverse might not exist
        XYZ C(_C.x, _C.y, _C.z);
        
        XYZFromRGB = rgb * glm::diagonal3x3(glm::vec3(C[0], C[1], C[2]));
        RGBFromXYZ = glm::inverse(XYZFromRGB);  //debug inverse might not exist
    }

    glm::mat3 ConvertRGBColorSpace(const RGBColorSpace& from, const RGBColorSpace& to)
    {
        if (from == to)
            return {};
        return to.RGBFromXYZ * from.XYZFromRGB;
    }

    RGBSigmoidPolynomial RGBColorSpace::ToRGBCoeffs(RGB rgb) const
    {
        if (!(rgb.r >= 0 && rgb.g >= 0 && rgb.b >= 0))
            std::cout << "RGBColorSpace::ToRGBCoeffs < 0\n";
        //DCHECK(rgb.r >= 0 && rgb.g >= 0 && rgb.b >= 0);
        return (*rgbToSpectrumTable)(ClampZero(rgb));
    }

    const RGBColorSpace* RGBColorSpace::GetNamed(std::string n)
    {
        std::string name;
        std::transform(n.begin(), n.end(), std::back_inserter(name), ::tolower);
        if (name == "aces2065-1")
            return ACES2065_1;
        else if (name == "rec2020")
            return Rec2020;
        else if (name == "dci-p3")
            return DCI_P3;
        else if (name == "srgb")
            return sRGB;
        else
            return nullptr;
    }

    const RGBColorSpace* RGBColorSpace::Lookup(glm::vec2 r, glm::vec2 g, glm::vec2 b, glm::vec2 w)
    {
        auto closeEnough = [](const glm::vec2& a, const glm::vec2& b) {
            return ((a.x == b.x || std::abs((a.x - b.x) / b.x) < 1e-3) &&
                (a.y == b.y || std::abs((a.y - b.y) / b.y) < 1e-3));
        };
        for (const RGBColorSpace* cs : { ACES2065_1, DCI_P3, Rec2020, sRGB }) {
            if (closeEnough(r, cs->r) && closeEnough(g, cs->g) && closeEnough(b, cs->b) &&
                closeEnough(w, cs->w))
                return cs;
        }
        return nullptr;
    }

    const RGBColorSpace* RGBColorSpace::sRGB;
    const RGBColorSpace* RGBColorSpace::DCI_P3;
    const RGBColorSpace* RGBColorSpace::Rec2020;
    const RGBColorSpace* RGBColorSpace::ACES2065_1;
    //LMS-eye cone color space https://www.strollswithmydog.com/cone-fundamental-lms-color-space/
    //CIE 1931 monochromatic color matching functions, 435.8nm, 546.1nm, 700nm

    void RGBColorSpace::Init()
    {
         // Rec. ITU-R BT.709.3

         sRGB = new RGBColorSpace(glm::vec2(.64, .33), glm::vec2(.3, .6), glm::vec2(.15, .06),
                                  GetNamedSpectrum("stdillum-D65"), RGBToSpectrumTable::sRGB);
         
         /*
         // P3-D65 (display)
         DCI_P3 = new RGBColorSpace(glm::vec2(.68, .32), glm::vec2(.265, .690), glm::vec2(.15, .06),
                                    GetNamedSpectrum("stdillum-D65"), RGBToSpectrumTable::DCI_P3);
         // ITU-R Rec BT.2020
         Rec2020 = new RGBColorSpace(glm::vec2(.708, .292), glm::vec2(.170, .797), glm::vec2(.131, .046),
                                     GetNamedSpectrum("stdillum-D65"), RGBToSpectrumTable::Rec2020);

         ACES2065_1 = new RGBColorSpace(glm::vec2(.7347, .2653), glm::vec2(0., 1.), glm::vec2(.0001, -.077),
                                        GetNamedSpectrum("illum-acesD60"), RGBToSpectrumTable::ACES2065_1);
                                        */
    }

    std::string RGBColorSpace::ToString() const
    {
        return "[ RGBColorSpace r:" + std::to_string(r.x) + ", " + std::to_string(r.y) +
            " g: " + std::to_string(g.x) + ", " + std::to_string(g.y) +
            " b: " + std::to_string(b.x) + ", " + std::to_string(b.y) +
            " w: " + std::to_string(w.x) + ", " + std::to_string(w.y);
        // " illuminant: " + illuminant + " RGBToXYZ: " +
            // XYZFromRGB + " XYZToRGB: "+ RGBFromXYZ;
    // return StringPrintf("[ RGBColorSpace r: %s g: %s b: %s w: %s illuminant: "
       //  "%s RGBToXYZ: %s XYZToRGB: %s ]",
       //  r, g, b, w, illuminant, XYZFromRGB, RGBFromXYZ);
    }

}  // namespace pbrt