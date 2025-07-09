#include "HoneyComb.hpp"
#include <boost/math/special_functions/ellint_1.hpp>

namespace DOS {
    HoneyComb::HoneyComb(size_t N) 
        : Base(N, -3, 3)
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
        const _internal_precision dE = static_cast<_internal_precision>((_max_energy - _min_energy) / _dos.size());

        for(int k = 0; k < _dos.size(); ++k) {
            const _internal_precision z = _min_energy + k * dE;
            _dos[k] = (LONG_1_PI*LONG_1_PI) * (std::abs(z) / std::sqrt(capital_Z0(z))) 
                        * boost::math::ellint_1(std::sqrt(capital_Z1(z) / capital_Z0(z)));
        }

        double sum{};
		for (const auto& d : _dos)
			sum += d;
		sum *= dE;
		std::cout << "DOS Norm = " << sum << std::endl;
    }
}