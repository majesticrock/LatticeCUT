#include "HoneyComb.hpp"
#include <boost/math/special_functions/ellint_1.hpp>

namespace DOS {
    constexpr double factor_range = 10;

    HoneyComb::HoneyComb(size_t N, double E_F, double debye) 
        : Base(N, -1, 1, E_F, factor_range * debye)
    { }

    _internal_precision capital_Z0(double energy) {
        energy = std::abs(energy);
        if (energy < 1) {
            return (1.+energy)*(1.+energy) - 0.25*(energy*energy-1.)*(energy*energy-1.);
        }
        else {
            return 4. * energy;
        }
    }

    _internal_precision capital_Z1(double energy) {
        energy = std::abs(energy);
        if (energy < 1) {
            return 4. * energy;
        }
        else {
            return (1.+energy)*(1.+energy) - 0.25*(energy*energy-1.)*(energy*energy-1.);
        }
    }

    void HoneyComb::compute() 
    {
        const double ratio = 6. / _energies.total_range;
        for(int k = 0; k < _dos.size(); ++k) {
            const _internal_precision z = ratio * _energies.index_to_energy(k);
            _dos[k] = ratio * _energies.get_dE(k) * (LONG_1_PI * LONG_1_PI) * (std::abs(z) / std::sqrt(capital_Z0(z))) 
                        * boost::math::ellint_1(std::sqrt(capital_Z1(z) / capital_Z0(z)));
        }
    }
}