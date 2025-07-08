#include "FCC.hpp"
#include <complex>
#include <iostream>
#include <boost/math/special_functions/hypergeometric_pFq.hpp>

namespace DOS {
    using dos_complex = std::complex<_internal_precision>;

    FCC::FCC(size_t N, double band_width)
        : Base(N, -12. * band_width / 16., 4. * band_width / 16.)
    { }

    dos_complex k_minus(dos_complex const& one_over_z) {
        return one_half 
            - 2.L * one_over_z * std::pow(1.L + one_over_z, -three_halves) 
            - one_half * (1.L - one_over_z) * std::pow(1.L + one_over_z, -three_halves) * std::sqrt(1.L - 3.L * one_over_z);
    }

    // DOS formula from https://link.springer.com/article/10.1007/s10955-011-0257-0
    void FCC::compute()
    {
        const _internal_precision EPS = 1e-6;
        const _internal_precision dE = static_cast<_internal_precision>((_max_energy - _min_energy) / _dos.size());
        for(int k = 0; k < _dos.size(); ++k) {
            const dos_complex z = static_cast<_internal_precision>(16. / (_max_energy - _min_energy)) 
                * dos_complex{_min_energy + k * dE, EPS};
            const dos_complex one_over_z = 1.L / z;

            const dos_complex ellint = one_half * LONG_PI * boost::math::hypergeometric_pFq({one_half, one_half}, {_internal_precision{1}}, k_minus(z));
            const dos_complex G = FOUR_OVER_PI_SQR * one_over_z * std::pow(1.L + one_over_z, -three_halves)
                * (2.L - (1.L - std::sqrt(1.L - 3.L * one_over_z))) * ellint * ellint;

            _dos[k] = (-1. / _pi) * static_cast<double>(std::imag(G));
        }
        
        double sum{};
		for (const auto& d : _dos)
			sum += d;
		sum *= dE;
		std::cout << "DOS Norm = " << sum << std::endl;
    }
}