#pragma once
#include "Base.hpp"

/**
 * Dispersion: 2t [cos(x) + cos(y) + cos(z)]
 * Band: [-6t, 6t]
 */

namespace DOS {
    struct SimpleCubic : public Base {
        SimpleCubic(size_t N);

    protected:
        void compute() final;
    };
}