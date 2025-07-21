#include "BCC.hpp"
#include <complex>
#include <boost/math/special_functions/ellint_1.hpp>
#include <mrock/utility/Numerics/hypergeometric_2F1.hpp>
#include <cassert>

namespace DOS {
    using dos_complex = std::complex<_internal_precision>;
    constexpr double factor_range = 10;

    BCC::BCC(size_t N, double E_F, double debye)
        : Base(N, -1., 1., E_F, factor_range * debye)
    { }

    dos_complex elliptical_integral(dos_complex k) {
        return one_half * LONG_PI  * mrock::utility::Numerics::hypergeometric_2F1(one_half, one_half, _internal_precision{1}, k);
    }

    // https://link.springer.com/article/10.1023/A:1015753428116
    void BCC::compute()
    {
        const _internal_precision ratio = 2. / _energies.total_range;
        for(int k = 0; k < _dos.size(); ++k) {
            const dos_complex energy = ratio * _energies.index_to_energy(k);
            const dos_complex ell  = elliptical_integral(one_half + std::sqrt(1.L - 1.L / (energy * energy)));
            _dos[k] = ratio * _energies.get_dE(k) * FOUR_OVER_PI_SQR * LONG_1_PI * (ell * ell / energy).imag();
        }
    }
}