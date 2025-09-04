#pragma once
#include "Base.hpp"

/**
 * Dispersion: 2t [cos(x) + cos(y) + cos(z)]
 * Band: [-6t, 6t]
 * https://journals.aps.org/prb/abstract/10.1103/PhysRevB.56.13960
 */ 

namespace DOS {
    struct SimpleCubic : public Base {
        SimpleCubic(size_t N, double E_F, double debye);

    protected:
        void compute() final;
    };
}