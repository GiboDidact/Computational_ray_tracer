// pbrt is Copyright(c) 1998-2020 Matt Pharr, Wenzel Jakob, and Greg Humphreys.
// The pbrt source code is licensed under the Apache License, Version 2.0.
// SPDX: Apache-2.0
//EDITED

#include "../../pch.h"
#include "RNG.h"

#include <cinttypes>

namespace pbrt {

std::string RNG::ToString() const {
    return "RNG::ToString() rng.cpp";
    //  return StringPrintf("[ RNG state: %" PRIu64 " inc: %" PRIu64 " ]", state, inc);
}

}  // namespace pbrt
