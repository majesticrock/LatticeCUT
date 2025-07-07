#include "Base.hpp"

namespace DOS {
    Base::Base(size_t N, double min_energy, double max_energy)
        : _dos(N), _min_energy{min_energy}, _max_energy{max_energy}
    { }

    std::vector<double> const& Base::get_dos()
    {
        if (!_computed) {
            this->compute();
            this->_computed = true;
        }
        return this->_dos; 
    }
}