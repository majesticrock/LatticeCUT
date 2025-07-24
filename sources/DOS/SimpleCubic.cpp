#include "SimpleCubic.hpp"
#include <boost/math/special_functions/ellint_1.hpp>
#include <boost/math/special_functions/ellint_2.hpp>
#include <boost/math/special_functions/pow.hpp>
#include <boost/math/quadrature/tanh_sinh.hpp>
#include <boost/math/quadrature/gauss_kronrod.hpp>
#include <iostream>

namespace DOS {
	constexpr double factor_range = 25;

    SimpleCubic::SimpleCubic(size_t N, double E_F, double debye)
        : Base(N, -1, 1, E_F, factor_range * debye)
    { }

	_CONST_LONG_FLOATING CUT_OFF = 1e-6L;
    _CONST_LONG_FLOATING R_AT_2 = (LONG_PI - 4.L) / 8.L;
	_CONST_LONG_FLOATING R_AT_2_1 = (5.L * LONG_PI / 64.L - 0.25L); 

    inline _internal_precision R(abscissa_t x) {
		x *= 0.5;
		// Special, analytically known cases:
		if (abs(x) < CUT_OFF) {
			return LOG_4;
		}
		return static_cast<_internal_precision>(boost::math::ellint_1(sqrt_1_minus_x_squared(x)) + 0.5 * log(x * x));
	}

	inline _internal_precision derivative_R(const abscissa_t& x) {
		// Special, analytically known cases:
		if (abs(x) < CUT_OFF) {
			return 0.0;
		}
		if (abs(x + 2) < CUT_OFF) {
			return R_AT_2 + R_AT_2_1 * static_cast<_internal_precision>(x + 2);
		}
		if (abs(x - 2) < CUT_OFF) {
			return -R_AT_2 + R_AT_2_1 * static_cast<_internal_precision>(x - 2);
		}

		const abscissa_t ALPHA = sqrt_1_minus_x_squared(0.5 * x);
		return static_cast<_internal_precision>((x * x * (1 - boost::math::ellint_1(ALPHA)) + 4 * (boost::math::ellint_2(ALPHA) - 1)) / (x * (x + 2) * (x - 2)));
	}

    inline _internal_precision I_1(const abscissa_t& epsilon) {
		const abscissa_t lower_bound{ std::max(abscissa_t{ -1 }, abscissa_t{ -2 } - epsilon) };
		const abscissa_t upper_bound{ std::min(abscissa_t{ 1 }, abscissa_t{ 2 } - epsilon) };
		_internal_precision ret = static_cast<_internal_precision>(asin(upper_bound) * R(upper_bound + epsilon) - asin(lower_bound) * R(lower_bound + epsilon));

		auto integrand = [&epsilon](const abscissa_t& phi) {
			return asin(phi) * derivative_R(phi + epsilon);
			};

		ret -= static_cast<_internal_precision>(boost::math::quadrature::gauss_kronrod<abscissa_t, 30>::integrate(integrand, lower_bound, upper_bound, 10, 1e-12));
		return ret;
	}

	inline _internal_precision I_2(const _internal_precision& epsilon) {
		// For some magical reason this integral is constant for epsilon in [-1, 1]
		// I_2 = -pi * ln(4)
		if (epsilon >= -1 && epsilon <= 1) {
			return (-LONG_PI * LOG_4);
		}
		const _internal_precision lower_bound{ std::max(_internal_precision{ -1 }, _internal_precision{ -2 } - epsilon) };
		const _internal_precision upper_bound{ std::min(_internal_precision{ 1 }, _internal_precision{ 2 } - epsilon) };

		auto integrand = [&epsilon](const _internal_precision& phi) {
			return log(0.25 * (epsilon + phi) * (epsilon + phi)) / sqrt_1_minus_x_squared(phi);
			};
		boost::math::quadrature::tanh_sinh<_internal_precision> integrator;
		return 0.5 * static_cast<_internal_precision>(integrator.integrate(integrand, lower_bound, upper_bound));
	}

    void SimpleCubic::compute()
    {
		const double ratio = 6. / _energies.total_range;
		_dos.front() = 0.0;
		_dos.back() = 0.0;
        for (int k = 1; k < _dos.size() - 1; ++k) {
			const _internal_precision epsilon = ratio * _energies.index_to_energy(k);

			try {
				_dos[k] = ratio * _energies.get_dE(k) * static_cast<double>(boost::math::pow<3>(LONG_1_PI) * (I_1(epsilon) - I_2(epsilon)));
			} catch (const boost::wrapexcept<std::overflow_error>& ex) {
				std::cerr << "Died while computing elliptic integrals in sc-DOS at epsilon=" << epsilon << std::endl;
				throw ex;
			}
		}
    }
}