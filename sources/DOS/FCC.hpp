#pragma once
#include "Base.hpp"

/**
 * Dispersion: 4t [ cos(x/2) cos(y/2) + cos(x/2) cos(z/2) + cos(y/2) cos(z/2) ]
 * Band: From -12t to 4t
 * 
 * DOS formula from  On the Exact Evaluation of the Face-Centred Cubic Lattice Green Function, GS Joyce, RT Delves
 * https://link.springer.com/article/10.1007/s10955-011-0257-0
 */

namespace DOS {
    struct FCC : public Base {
        FCC(size_t N);

    protected:
        void compute() final;
    };
}