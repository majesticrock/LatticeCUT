#include "Selector.hpp"

#include <stdexcept>
#include <cctype>
#include <algorithm>
#include <numeric>
#include <iostream>

#include "FreeElectrons.hpp"
#include "SimpleCubic.hpp"
#include "FCC.hpp"
#include "BCC.hpp"
#include "HoneyComb.hpp"

constexpr double FREE_3D_MAX_ENERGY = 5.0;

namespace DOS {
    std::vector<double> const& Selector::select_dos(const std::string& name, int N, double E_F, double debye) 
    {
        const std::string free_electrons_name = "free_electrons";
        const std::string sc_name = "sc";
        const std::string fcc_name = "fcc";
        const std::string bcc_name = "bcc";
        const std::string honeycomb_name = "hc";

        if (dos_ptr) return dos_ptr->get_dos();

        if (name.rfind(free_electrons_name, 0) == 0) {
            std::string dimension_str = name.substr(free_electrons_name.length());
            // std::isdigit does not seem to work for some reason unbeknownst to me
            if (dimension_str.empty() || !std::all_of(dimension_str.begin(), dimension_str.end(), ::isdigit));
            
            dos_ptr = std::make_unique<FreeElectrons>(N, FREE_3D_MAX_ENERGY, std::stoi(dimension_str), E_F, debye);
        }
        else if (name == sc_name) {
            dos_ptr = std::make_unique<SimpleCubic>(N, E_F, debye);
        }
        else if (name == fcc_name) {
            dos_ptr = std::make_unique<FCC>(N, E_F, debye);
        }
        else if (name == bcc_name) {
            dos_ptr = std::make_unique<BCC>(N, E_F, debye);
        }
        else if (name == honeycomb_name) {
            dos_ptr = std::make_unique<HoneyComb>(N, E_F, debye);
        }
        else {
            throw std::invalid_argument("DOS '" + name + "' not recognized!");
        }

        return dos_ptr->get_dos();
    }

    std::vector<double> Selector::get_raw_dos() const
    {
        std::vector<double> ret(this->dos_ptr->get_dos());
        for(int i = 0; i < ret.size(); ++i) {
            ret[i] /= this->dos_ptr->get_energies().get_dE(i);
        }
        return ret;
    }

    double Selector::average_in_range(double low, double up) const
    {
        int N_low = dos_ptr->get_energies().energy_to_index(low);
        int N_up  = dos_ptr->get_energies().energy_to_index(up);

        const double avg{ std::reduce(dos_ptr->get_dos().begin() + N_low, dos_ptr->get_dos().begin() + N_up, double{}) / (up - low) };
        std::cout << "Average DOS in range = " << avg << std::endl;
        return avg;
    }
}