#include "FCC.hpp"
#include <complex>
#include <iostream>
#include <mrock/utility/Numerics/hypergeometric_2F1.hpp>

namespace DOS {
    using dos_complex = std::complex<_internal_precision>;
    constexpr double factor_range = 10;

    FCC::FCC(size_t N, double E_F, double debye)
        : Base(N, -0.999, 3, E_F, factor_range * debye)
    { }

    dos_complex k_minus(dos_complex const& one_over_z) {
        return one_half - 2.L * one_over_z * std::pow(1.L + one_over_z, -three_halves) 
            - one_half * (1.L - one_over_z) * std::pow(1.L + one_over_z, -three_halves) * std::sqrt(1.L - 3.L * one_over_z);
    }

    //DOS formula from https://link.springer.com/article/10.1007/s10955-011-0257-0
    void FCC::compute()
    {
        const _internal_precision EPS = 1e-10;
        const _internal_precision ratio = 4. / _energies.total_range;
        for(int k = 0; k < _dos.size(); ++k) {
            const dos_complex z = dos_complex{ratio * _energies.index_to_energy(k), EPS};
            const dos_complex one_over_z = 1.L / z;

            const dos_complex ellint = one_half * LONG_PI 
                * mrock::utility::Numerics::hypergeometric_2F1(one_half, one_half, _internal_precision{1}, k_minus(one_over_z));
            const dos_complex G = FOUR_OVER_PI_SQR * one_over_z * std::pow(1.L + one_over_z, -three_halves)
                * (2.L - std::sqrt(1.L - 3.L * one_over_z)) * ellint * ellint;

            _dos[k] =  ratio * _energies.get_dE(k) * (-1. / _pi) * static_cast<double>(std::imag(G));
        }
    }
}