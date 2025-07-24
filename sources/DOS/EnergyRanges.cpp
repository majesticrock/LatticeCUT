#include "EnergyRanges.hpp"
#include <cmath>
#include <cassert>
#include <iostream>
#include <iomanip>

DOS::EnergyRanges::EnergyRanges(double _E_min, double _E_max, double _inner_min, double _inner_max, int _N)
	: E_min{_E_min}, E_max{_E_max}, total_range{E_max - E_min},
#ifndef UNIFORM_DISCRETIZATION
	inner_min{std::max(E_min, _inner_min)}, inner_max{std::min(E_max, _inner_max)}, 
	total_outer_range{ E_max - inner_max + inner_min - E_min }, 
#endif
	N{_N},
#ifndef UNIFORM_DISCRETIZATION
	inner_discretization{static_cast<int>( inner_fraction  * N )}, 
	lower_discretization{static_cast<int>( outer_fraction  * N      * (inner_min - E_min)     / total_outer_range)}, 
	upper_discretization{static_cast<int>( (outer_fraction * N - 1) * (E_max     - inner_max) / total_outer_range)},
	inner_dE{(inner_max - inner_min) / inner_discretization},
	lower_dE{(inner_min - E_min)     / lower_discretization},
	upper_dE{(E_max     - inner_max) / upper_discretization}
#else
	dE{(E_max - E_min) / (N - 1)}
#endif
{
	std::cout << std::setprecision(17) << std::endl;
	std::cout << this->E_min << "   " << this->index_to_energy(0) << std::endl;
	std::cout << this->E_max << "   " << this->index_to_energy(N - 1) << std::endl;
}

int DOS::EnergyRanges::energy_to_index(double epsilon) const
{
#ifdef UNIFORM_DISCRETIZATION
	return static_cast<int>((epsilon - E_min) / dE);
#else
    if (epsilon < inner_min) {
		return static_cast<int>(std::lround((epsilon - E_min) / lower_dE));
	}
	if (epsilon < inner_max) {
		return static_cast<int>(std::lround((epsilon - inner_min) / inner_dE) + lower_discretization);
	}
	return static_cast<int>(std::lround((epsilon - inner_max) / upper_dE) + inner_discretization + lower_discretization);
#endif
}

double DOS::EnergyRanges::index_to_energy(int index) const
{
#ifdef UNIFORM_DISCRETIZATION
	return E_min + dE * index;
#else
    if (index < lower_discretization) {
		return E_min + lower_dE * index;
	}
	index -= lower_discretization;
	if (index < inner_discretization) {
		return inner_min + inner_dE * index;
	}
	index -= inner_discretization;
	return inner_max + upper_dE * index;
#endif
}

double DOS::EnergyRanges::get_dE(int index) const
{
#ifdef UNIFORM_DISCRETIZATION
	return dE;
#else
    if (index < lower_discretization) {
		return lower_dE;
	}
	if (index == lower_discretization){
		return 0.5 * (lower_dE + inner_dE);
	}
	index -= lower_discretization;
	if (index < inner_discretization) {
		return inner_dE;
	}
	if (index == inner_discretization) {
		return 0.5 * (inner_dE + upper_dE);
	}
	return upper_dE;
#endif
}

std::vector<double> DOS::EnergyRanges::get_abscissa() const
{
    std::vector<double> absc(N);
    for(int i = 0; i < N; ++i) {
        absc[i] = index_to_energy(i);
    }
    return absc;
}