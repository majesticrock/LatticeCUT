#include "FreeElectrons.hpp"
#include <cassert>

namespace DOS {
    constexpr double factor_range = 10;

    FreeElectrons::FreeElectrons(size_t N, double band_width, int _dimension, double E_F, double debye) 
        : Base(N, 0, band_width, E_F, factor_range * debye), dimension(_dimension)
    { }

    double FreeElectrons::normalization(double mass/* =1 */, double hbar/* =1 */) const
    {
        return std::pow(2 * mass / (hbar * hbar), 0.5 * dimension) / (std::tgamma(0.5 * dimension) * std::pow(2 * _pi, dimension));
    }

    void FreeElectrons::compute() 
    {
        const double An = normalization(); 
        for (size_t i = 0; i < _dos.size(); i++)
        {
            _dos[i] = _energies.get_dE(i) * An * std::pow(_energies.index_to_energy(i), 0.5 * dimension - 1.0);
        }
    }
}