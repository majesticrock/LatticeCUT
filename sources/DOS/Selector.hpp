#pragma once
#include <vector>
#include <string>
#include <memory>
#include <cassert>

#include "Base.hpp"

namespace DOS {
    struct Selector {
        std::vector<double> const& select_dos(const std::string& name, int N, double E_F, double debye);

        inline double get_min_energy() const {
            assert(dos_ptr && "dos_ptr is nullptr!");
            return dos_ptr->get_min_energy();
        }
        inline double get_max_energy() const {
            assert(dos_ptr && "dos_ptr is nullptr!");
            return dos_ptr->get_max_energy();
        }
        inline double get_band_width() const {
            return (get_max_energy() - get_min_energy());
        }

        inline const EnergyRanges& get_energies() const {
            assert(dos_ptr != nullptr);
            return dos_ptr->get_energies();
        }

        // We include dE in the dos-vector. Thus, if we want access to the plain dos values
        // we need to account for that. This function does exactly that and returns the plain dos values
        std::vector<double> get_raw_dos() const;

        double average_in_range(double low, double up) const;
    private:
        std::unique_ptr<Base> dos_ptr;
    };
}