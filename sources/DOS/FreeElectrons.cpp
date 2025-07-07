#include "FreeElectrons.hpp"
#include <cassert>

namespace DOS {
    FreeElectrons::FreeElectrons(size_t N, double bandwidth, int _dimension) 
        : Base(N, 0, bandwidth), dimension(_dimension)
    { }

    double FreeElectrons::normalization(double mass/* =1 */, double hbar/* =1 */) const
    {
        return std::pow(2 * mass / (hbar * hbar), 0.5 * dimension) / (std::tgamma(0.5 * dimension) * std::pow(2 * _pi, dimension));
    }

    void FreeElectrons::compute() 
    {
        const double An = normalization(); 
        const double step = (_max_energy - _min_energy) / (_dos.size() - 1);
        for (size_t i = 0; i < _dos.size(); i++)
        {
            _dos[i] = An * std::pow(i * step + _min_energy, 0.5 * dimension - 1.0);
        }
    }
}