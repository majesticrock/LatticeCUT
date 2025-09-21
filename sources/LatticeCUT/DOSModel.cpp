#include "DOSModel.hpp"

#include <mrock/utility/better_to_string.hpp>

#include "../DOS/Selector.hpp"

namespace LatticeCUT {
    DOSModel::DOSModel(mrock::utility::InputFileReader& input) 
        : phonon_coupling_in{input.getDouble("phonon_coupling")},
        local_interaction{input.getDouble("local_interaction")},
        fermi_energy{ input.getDouble("fermi_energy") },
        omega_debye_in{ input.getDouble("omega_debye") },
        N{input.getInt("N")},
        selector(),
        dos_name{input.getString("dos")},
        density_of_states{selector.select_dos(dos_name, N, fermi_energy, omega_debye_in)},
        energies{selector.get_energies()},
        omega_debye{ omega_debye_in * energies.total_range },
        phonon_coupling{phonon_coupling_in / selector.average_in_range(fermi_energy - omega_debye, fermi_energy + omega_debye)},
        local_interaction_energy_units{local_interaction / selector.average_in_range(fermi_energy - omega_debye, fermi_energy + omega_debye)},
        Delta(decltype(Delta)::FromAllocator([&](int k) -> l_float {
			const l_float magnitude = (k < omega_debye_in * N || k > omega_debye_in * N) ? 0.01 : 0.1;
			if (k < N) {
				return magnitude;
			}
			return l_float{};
			}, N))
    {
    }

    void DOSModel::iteration_step(const ParameterVector& initial_values, ParameterVector& result)
    {
        static int step_num = 0;
        result.setZero();
        this->Delta.fill_with(initial_values);
        this->get_expectation_values();

        for (int k = 0; k < N; k++)
        {
            const l_float energy_k = energies.index_to_energy(k);
            l_float __part{};
            for (int l = phonon_lower_bound(energy_k); l <= phonon_upper_bound(energy_k); ++l)
            {
                __part -= _expecs[mrock::symbolic_operators::SC_Type][l] * density_of_states[l];
            }
            result(k) = phonon_coupling * __part;
            __part = l_float{};
            for (int l = 0; l < N; ++l)
            {
                __part += _expecs[mrock::symbolic_operators::SC_Type][l] * density_of_states[l];
            }
            result(k) += local_interaction_energy_units * __part;
        }

        this->Delta.fill_with(result, 0.5);
        //std::cout << step_num << ": " << delta_max() << std::endl;
		this->Delta.clear_noise(PRECISION);
		result -= initial_values;
		++step_num;
    }

    l_float DOSModel::sc_expectation_value_index(int k) const
    {
        if (is_zero(Delta[k])) return l_float{};
		return -Delta[k] / (2 * quasiparticle_energy_index(k));
    }

    l_float DOSModel::occupation_index(int k) const
    {
        const l_float eps_mu = single_particle_energy(k);
		if (is_zero(Delta[k])) {
			if (is_zero(eps_mu)) return 0.5;
			return (eps_mu < 0 ? 1 : 0);
		}
		return 0.5 * (1 - (eps_mu / sqrt(eps_mu * eps_mu + std::norm(Delta[k]))));
    }

    l_float DOSModel::compute_coefficient(mrock::symbolic_operators::Coefficient const &coeff, int first, int second) const
    {
        if (coeff.name == "\\epsilon_0")
		{
			return single_particle_energy(first);
		}
		else if (coeff.name == "g")
		{
            if (std::abs(energies.index_to_energy(first) - energies.index_to_energy(second)) > omega_debye) {
                return l_float{};
            }
            return phonon_coupling;
        }
        else if (coeff.name == "U")
        {
            return local_interaction_energy_units;
        }
        else
		{
			throw std::invalid_argument("Coefficient not recognized! " + coeff.name);
		}
    }

    const std::map<mrock::symbolic_operators::OperatorType, std::vector<l_float>> &DOSModel::get_expectation_values() const
    {
		if (_expecs.empty()) {
			_expecs.emplace(mrock::symbolic_operators::Number_Type, std::vector<l_float>(N));
			_expecs.emplace(mrock::symbolic_operators::SC_Type, std::vector<l_float>(N));
		}
		for (int k = 0; k < N; ++k) {
			_expecs.at(mrock::symbolic_operators::Number_Type)[k] = this->occupation_index(k);
			_expecs.at(mrock::symbolic_operators::SC_Type)[k] = this->sc_expectation_value_index(k);
		}

		return _expecs;
	}

    l_float DOSModel::delta_max() const
    {
        return std::abs(*std::max_element(Delta.begin(), Delta.begin() + N,
			[](const l_float& lhs, const l_float& rhs) {
				return std::abs(lhs) < std::abs(rhs);
			}
		));
    }

    std::string DOSModel::info() const
    {
        return "DOSModel: g=" + std::to_string(phonon_coupling) + "\tU=" + std::to_string(local_interaction_energy_units) + "\tDOS=" + dos_name;
    }

    std::string DOSModel::to_folder() const
    {
        auto improved_string = [](l_float number) -> std::string {
			if (std::floor(number) == number) {
				// If the number is a whole number, format it with one decimal place
				std::ostringstream out;
				out.precision(1);
				out << std::fixed << number;
				return out.str();
			}
			else {
				std::string str = mrock::utility::better_to_string(number, std::chars_format::fixed);
				// Remove trailing zeroes
				str.erase(str.find_last_not_of('0') + 1, std::string::npos);
				str.erase(str.find_last_not_of('.') + 1, std::string::npos);
				return str;
			}
			};
        
        return dos_name + "/N=" + std::to_string(N) + "/"
            + "g=" + improved_string(phonon_coupling_in) + "/"
            + "U=" + improved_string(local_interaction) + "/"
            + "E_F=" + improved_string(fermi_energy) + "/"
            + "omega_D=" + improved_string(omega_debye_in) + "/";
    }
}