#pragma once

#include <array>
#include <cstddef>

namespace fd {

    // Generic stencil container
    template<std::size_t N>
    struct Stencil {
        std::array<int, N> offset;      // grid offsets
        std::array<double, N> coeff;    // weights
    };

    // -----------------------------
    // First derivative, 2nd order
    // Central difference
    // -----------------------------
    constexpr Stencil<3> d1_central_2 = {
        { -1, 0, 1 },
        { -0.5, 0.0, 0.5 }
    };

    // -----------------------------
    // First derivative, 4th order
    // Central difference
    // -----------------------------
    constexpr Stencil<5> d1_central_4 = {
        { -2, -1, 0, 1, 2 },
        {  1.0/12.0, -2.0/3.0, 0.0, 2.0/3.0, -1.0/12.0 }
    };

    // -----------------------------
    // Second derivative, 2nd order
    // Central difference
    // -----------------------------
    constexpr Stencil<3> d2_central_2 = {
        { -1, 0, 1 },
        { 1.0, -2.0, 1.0 }
    };

    // -----------------------------
    // Second derivative, 4th order
    // Central difference
    // -----------------------------
    constexpr fd::Stencil<5> d2_central_4 = {  
        { -2, -1, 0, 1, 2 },
        { -1.0/12.0, 4.0/3.0, -5.0/2.0, 4.0/3.0, -1.0/12.0 }
    };

}
