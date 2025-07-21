#pragma once
#include "Base.hpp"

/**
 * DOS formula from Exact evaluation of the body centred cubic lattice Green function, G S Joyce
 * Dynamical mean field study of the Dirac liquid, S.A. Jafari
 * https://link.springer.com/article/10.1140/epjb/e2009-00128-1
 * Band from [-3, 3]
 */

namespace DOS {
    struct HoneyComb : public Base {
        HoneyComb(size_t N, double E_F, double debye);

    protected:
        void compute() final;
    };
}