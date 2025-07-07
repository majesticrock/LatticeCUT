#include "Selector.hpp"

#include <stdexcept>
#include <cctype>
#include <algorithm>

#include "FreeElectrons.hpp"

namespace DOS {
    std::vector<double> const& Selector::select_dos(const std::string& name, int N, double bandwidth) 
    {
        const std::string free_electrons_name = "free_electrons";

        if (dos_ptr) return dos_ptr->get_dos();

        if (name.rfind(free_electrons_name, 0) == 0) {
            std::string dimension_str = name.substr(free_electrons_name.length());
            // std::isdigit does not seem to work for some reason unbeknownst to me
            if (dimension_str.empty() || !std::all_of(dimension_str.begin(), dimension_str.end(), ::isdigit));
            
            dos_ptr = std::make_unique<FreeElectrons>(N, bandwidth, std::stoi(dimension_str));
        }
        else {
            throw std::invalid_argument("DOS '" + name + "' not recognized!");
        }

        return dos_ptr->get_dos();
    }
}