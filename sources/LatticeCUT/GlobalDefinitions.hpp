#pragma once
#define _USE_MATH_DEFINES

#include <Eigen/Dense>
#include <cmath>
#include <complex>
#include <stdint.h>
#include <limits>

namespace LatticeCUT {
    using l_float = double;
    using l_complex = std::complex<l_float>;

    constexpr l_complex imaginary_unit{0., 1.};
    constexpr l_float pi   = l_float(M_PI); // pi

    constexpr l_float fermi_function(l_float energy, l_float beta) noexcept {
        if (beta < 0) { // negative values are treated as infinity
            return (energy != 0 ? (energy < 0 ? l_float{1} : l_float{}) : l_float{0.5} );
        }
        return 1. / (1. + std::exp(beta * energy));
    }

    template<class... Args>
    l_float norm(Args... args) {
        return std::sqrt((... + (args * args)));
    }

    constexpr l_float SQRT_PRECISION = 2.384185791015625e-07;
	constexpr l_float PRECISION = 5.684341886080802e-14;

    /* This function abuses the structure of our desired precision:
	*  The mantissa is empty, i.e., we can solely rely on the exponent.
	*  If the exponent of <number> is >= 0b01111010011, |number| >= precision, i.e., not 0
	*/
	inline bool is_zero(double number) {
		static_assert(std::numeric_limits<double>::is_iec559, "IEEE 754 floating point not verified!");
		// 5.684341886080802e-14   <->   0 | 01111010011 | 0000000000000000000000000000000000000000000000000000
		uint64_t tmp; // silence compiler warnings (at 0 overhead)
		std::memcpy(&tmp, &number, sizeof(tmp));
		return static_cast<uint16_t>((tmp >> 52) & 0x7FF) < 0b01111010011;
	}

	inline bool is_zero(std::complex<double> number) {
		return is_zero(std::abs(number));
	}
}