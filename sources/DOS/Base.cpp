#include "Base.hpp"
#include <iostream>

namespace DOS {
    Base::Base(size_t N, double min_energy, double max_energy)
        : _dos(N), _min_energy{min_energy}, _max_energy{max_energy}
    { }

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
        const double dE = static_cast<double>((_max_energy - _min_energy) / _dos.size());
        double sum{};
		for (const auto& d : _dos)
			sum += d;
		sum *= dE;
		std::cout << "DOS Norm = " << sum << std::endl;
        for (auto& d : _dos)
            d /= sum;
    }
}