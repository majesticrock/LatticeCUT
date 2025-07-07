#pragma once
#include "Base.hpp"

namespace DOS {
    struct SimpleCubic : public Base {
        SimpleCubic(size_t N, double band_width);

    protected:
        void compute() final;
    };
}