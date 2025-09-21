#pragma once
#include "Base.hpp"

namespace DOS {
    struct DoublePeak : public Base {
        DoublePeak(size_t N, double E_F, double debye, double peak_weight);

    protected:
        const double rel_peak_weight;
        void compute() final;
    };
}