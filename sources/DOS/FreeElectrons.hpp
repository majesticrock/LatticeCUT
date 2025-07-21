#pragma once
#include "Base.hpp"


/**
 * Computes the DOS for free electrons in N dimensions.
 * Formula from doi:10.1119/1.10859
 */

namespace DOS {
    struct FreeElectrons : public Base {
        const int dimension;

        FreeElectrons(size_t N, double band_width, int _dimension, double E_F, double debye);

        double normalization(double mass=1, double hbar=1) const;

    protected:
        void compute() final;
    };
}