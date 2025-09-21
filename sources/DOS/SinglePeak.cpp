#include "SinglePeak.hpp"
#include <numeric>

namespace DOS {
    constexpr double factor_range = 25;
    SinglePeak::SinglePeak(size_t N, double E_F, double debye, double peak_weight)
        : Base(N, -1, 1, E_F, factor_range * debye), rel_peak_weight{ peak_weight / 100.0 }
    { }

    void SinglePeak::compute()
    { 
        constexpr double peak_width = 0.01;
        const double rel_peak_height = rel_peak_weight / peak_width * (1. - peak_width) / (1. - rel_peak_weight);

        for(int k = 0; k < _dos.size(); ++k) {
            const double energy = _energies.index_to_energy(k);
            if (std::abs(energy) < peak_width) {
                _dos[k] = rel_peak_height;
            }
            else {
                _dos[k] = 1.0;
            }
        }

        const double dos_norm = std::reduce(_dos.begin(), _dos.end(), double{});
        for(auto& val : _dos) {
            val /= dos_norm;
        }
    }
}