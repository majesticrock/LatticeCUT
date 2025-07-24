#pragma once
#include <vector>

#define UNIFORM_DISCRETIZATION

namespace DOS {
    struct EnergyRanges {
        constexpr static double inner_fraction = 0.5;
        constexpr static double outer_fraction = 1. - inner_fraction;
        
        double E_min{};
        double E_max{};
        double total_range{};

#ifndef UNIFORM_DISCRETIZATION
        double inner_min{};
        double inner_max{};
        double total_outer_range{};
#endif

        int N{};
#ifndef UNIFORM_DISCRETIZATION
        int inner_discretization{};
        int lower_discretization{};
        int upper_discretization{};

        double inner_dE{};
        double lower_dE{};
        double upper_dE{};
#else
        double dE{};
#endif

        EnergyRanges() = default;
        EnergyRanges(double _E_min, double _E_max, double _inner_min, double _inner_max, int _N);

        int energy_to_index(double epsilon) const;
        double index_to_energy(int index) const;
        double get_dE(int index) const;
        
        std::vector<double> get_abscissa() const;
    };
}