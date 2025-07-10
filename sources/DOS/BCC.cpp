#include "BCC.hpp"
#include <complex>
#include <boost/math/special_functions/ellint_1.hpp>
#include <mrock/utility/Numerics/hypergeometric_2F1.hpp>
#include <cassert>

namespace DOS {
    using dos_complex = std::complex<_internal_precision>;

    BCC::BCC(size_t N)
        : Base(N, -1., 1.)
    { }

    dos_complex elliptical_integral(dos_complex k) {
        return one_half * LONG_PI  * mrock::utility::Numerics::hypergeometric_2F1(one_half, one_half, _internal_precision{1}, k);
    }

    // https://link.springer.com/article/10.1023/A:1015753428116
    void BCC::compute()
    {
        const _internal_precision dE = static_cast<_internal_precision>((_max_energy - _min_energy) / _dos.size());
        for(int k = 0; k < _dos.size(); ++k) {
            const dos_complex energy = {_min_energy + dE * k, 1e-8};
            const dos_complex ell  = elliptical_integral(one_half + std::sqrt(1.L - 1.L / (energy * energy)));
            _dos[k] = FOUR_OVER_PI_SQR * LONG_1_PI * (ell * ell / energy).imag();
        }
    }
}