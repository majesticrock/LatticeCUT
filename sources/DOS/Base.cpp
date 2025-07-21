#include "Base.hpp"
#include <iostream>
#include <numeric>

namespace DOS {
    Base::Base(size_t N, double min_energy, double max_energy, double fine_center, double fine_range)
        : _energies(min_energy, max_energy, fine_center - (max_energy - min_energy) * fine_range, fine_center + (max_energy - min_energy) * fine_range, N), _dos(N)
    {
    }

    std::vector<double> const& Base::get_dos()
    {
        if (!_computed) {
            this->compute();
            this->normalize();
            this->_computed = true;
        }
        return this->_dos; 
    }

    void Base::normalize() {
        const double sum{std::reduce(_dos.begin(), _dos.end(), double{})};
		std::cout << "DOS Norm = " << sum << std::endl;
        for (auto& d : _dos)
            d /= sum;
    }
}