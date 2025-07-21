#pragma once
#include <vector>
#include <cmath>
#include <boost/math/constants/constants.hpp>
#include "EnergyRanges.hpp"

#define _CONST_LONG_FLOATING constexpr _internal_precision

namespace DOS {
    using _internal_precision = long double;
    using abscissa_t = long double;

    _CONST_LONG_FLOATING one_half = 0.5L; // 1/2
    _CONST_LONG_FLOATING three_halves = 1.5L; // 3/2
    _CONST_LONG_FLOATING LONG_1_PI = 0.5L * boost::math::constants::two_div_pi<_internal_precision>(); // 1 / pi
    _CONST_LONG_FLOATING LONG_PI = boost::math::constants::pi<_internal_precision>(); // pi
    _CONST_LONG_FLOATING FOUR_OVER_PI_SQR = boost::math::constants::two_div_pi<_internal_precision>() * boost::math::constants::two_div_pi<_internal_precision>(); // 4 / pi^2
    _CONST_LONG_FLOATING LOG_4 = 2.L * boost::math::constants::ln_two<_internal_precision>(); // ln(4) = 2 ln(2)

    template <class RealType>
	inline RealType sqrt_1_minus_x_squared(const RealType& x) {
		return sqrt((1 - x) * (1 + x));
	}

    struct Base {
        virtual ~Base() = default;
        Base(size_t N, double min_energy, double max_energy, double fine_center, double fine_range);

        // The integration weight is contained in the returned vector
        // so that a simple summation of the vector elements suffices to compute an integral
        // Must hence be used in unison with the EnergyRanges provided in the class
        std::vector<double> const& get_dos();
        inline EnergyRanges const& get_energies() const
        {
            return this->_energies;
        }
        
        void normalize();

        inline double get_min_energy() const noexcept {
            return _energies.E_min;
        }
        inline double get_max_energy() const noexcept {
            return _energies.E_min;
        }
        inline size_t size() const noexcept {
            return _dos.size();
        }
    protected:
        static constexpr double _pi = 3.1415926535897932384626;
        virtual void compute() = 0;
        
        EnergyRanges _energies;
        std::vector<double> _dos;
        bool _computed{};
    };
}