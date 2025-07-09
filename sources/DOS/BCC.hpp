#pragma once
#include "Base.hpp"

/**
 * DOS formula from Exact evaluation of the body centred cubic lattice Green function, G S Joyce
 * https://iopscience.iop.org/article/10.1088/0022-3719/4/12/008/meta
 * https://www.jstor.org/stable/52609 -> The P(z) in the paper is not our Green's function, but rather G(z) = P(1/z)/z
 * Band from [-1, 1]
 */

namespace DOS {
    struct BCC : public Base {
        BCC(size_t N);

    protected:
        void compute() final;
    };
}