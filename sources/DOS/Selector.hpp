#pragma once
#include <vector>
#include <string>
#include <memory>
#include <cassert>

#include "Base.hpp"

namespace DOS {
    struct Selector {
        std::vector<double> const& select_dos(const std::string& name, int N, double bandwidth);

        inline double get_min_energy() const {
            assert(dos_ptr && "dos_ptr is nullptr!");
            return dos_ptr->get_min_energy();
        }
        inline double get_max_energy() const {
            assert(dos_ptr && "dos_ptr is nullptr!");
            return dos_ptr->get_max_energy();
        }
    private:
        std::unique_ptr<Base> dos_ptr;
    };
}