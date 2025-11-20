#include "DOSModel.hpp"

#include <mrock/utility/better_to_string.hpp>
#include <boost/math/tools/roots.hpp>

#include "../DOS/Selector.hpp"
#include <limits>

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
        beta{ input.getDouble("beta") },
        chemical_potential{ fermi_energy },
        filling_at_zero_temp{ 0. },
        Delta(decltype(Delta)::FromAllocator([&](int k) -> l_float {
            if (k == N) {
                return fermi_energy;
            }
            const double index_at_ef = 0.5 * N * (fermi_energy + 1); // only used for comparison, so double is fine
            const double index_at_0 = N / 2;
            const double range = omega_debye_in * N;
            
            double magnitude = local_interaction > 0.0 ? -0.001 : 0.001;
            if (k < index_at_ef + range && k > index_at_ef - range) {
                magnitude += 0.1;
            }
            else if (k < index_at_0 + range && k > index_at_0 - range) {
                magnitude += local_interaction > 0.0 ? -0.1 : 0.1;
            }
            return magnitude;
			}, N + 1))
    {
        int i = 0;
        l_float eps{this->single_particle_energy(0)};
        
        while((eps = this->single_particle_energy(++i)) < l_float{}) {
            filling_at_zero_temp += density_of_states[i];
        }
        if (eps == l_float{}) {
            filling_at_zero_temp += 0.5 * density_of_states[i];
        }
        std::cout << "Filling at zero temperature: " << filling_at_zero_temp << std::endl;
    }

    void DOSModel::iteration_step(const ParameterVector& initial_values, ParameterVector& result)
    {
        static int step_num = 0;
        result.setZero();
        this->Delta.fill_with(initial_values);
        this->chemical_potential = initial_values(N);
        this->get_expectation_values();

        auto fit_occupation = [&](const l_float mu) -> l_float {
            return filling_at_zero_temp - compute_filling(mu);
        };
        auto bracket_root = [&]() -> std::pair<l_float, l_float> {
            constexpr l_float step = 0.005;
            l_float a{chemical_potential - step};
            l_float b{chemical_potential + step};

            l_float fa{fit_occupation(a)};
            l_float fb{fit_occupation(b)};

            l_float dfa, dfb;

            if (fa * fb < 0.0) {
                return std::pair<l_float, l_float>{a, b};
            }
            else {
                dfa = fit_occupation(a - step) - fa;
                dfb = fit_occupation(b + step) - fb;
            }

            while (fa * fb > 0.0) {
                if (fb * dfb < 0.0){
                    b += step;
                    fb = fit_occupation(b);
                    dfb = fit_occupation(b + step) - fb;
                }
                if (fa * dfa < 0.0) {
                    a -= step;
                    fa = fit_occupation(a);
                    dfa = fit_occupation(a - step) - fa;
                }
            }
            //std::cout << "Found [" << a << "," << b << "]: " << fa << "\t" << fb << std::endl;
            return std::pair<l_float, l_float>{a, b};
        };

        std::uintmax_t boost_max_it{100U};
        const auto bracket = bracket_root();
        const auto best_mu = boost::math::tools::toms748_solve(fit_occupation, 
                bracket.first, bracket.second,
                boost::math::tools::eps_tolerance<l_float>(), boost_max_it);
                
        result(N) = 0.5 * (best_mu.first + best_mu.second);
        std::cout << std::setprecision(16);
        //std::cout << "Step: " << step_num << " found mu=" << result(N) << std::endl;

#pragma omp parallel for
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
        this->Delta[N] = 0.5 * (this->chemical_potential + initial_values(N));
        this->chemical_potential = this->Delta[N];
        //std::cout << step_num << ": " << delta_max() << std::endl;
		this->Delta.clear_noise(PRECISION);
		result -= initial_values;
		++step_num;
    }

    l_float DOSModel::sc_expectation_value_index(int k) const
    {
        if (is_zero(Delta[k])) return l_float{};
        if (beta < 0)
		    return -Delta[k] / (2 * quasiparticle_energy_index(k));
        else
            return -Delta[k] / (2 * quasiparticle_energy_index(k)) * std::tanh(0.5 * beta * quasiparticle_energy_index(k));
    }

    l_float DOSModel::occupation_index(int k) const
    {
        const l_float eps_mu = single_particle_energy(k);
		if (is_zero(Delta[k])) {
			return fermi_function(eps_mu, beta);
		}
        const l_float E = quasiparticle_energy_index(k);
        if (beta < 0)
            return 0.5 * (1 - (eps_mu / E));
        else
		    return 0.5 * (1 - (eps_mu / E) * std::tanh(0.5 * beta * E));
    }

    l_float DOSModel::compute_filling(const l_float mu) const
    {
        auto number_operator = [&](int i) -> l_float {
            const l_float eps = energies.index_to_energy(i) - mu;
            const l_float E = sqrt(eps*eps + Delta[i]*Delta[i]);
            if (is_zero(Delta[i]))
			    return fermi_function(eps, beta);
            if (beta < 0)
                return 0.5 * (1 - (eps / E));
            else
		        return 0.5 * (1 - (eps / E) * std::tanh(0.5 * beta * E));
        };
        l_float n{};
        for(int i = 0U; i < N; ++i) {
            
            n += density_of_states[i] * number_operator(i);
        }
        return n;
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
#pragma omp parallel for
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

    l_float DOSModel::true_gap() const
    {
        l_float current_min = 1000.;
        for (int k = 0; k < N; ++k) {
            const l_float E = quasiparticle_energy_index(k);
            if (E < current_min) {
                current_min = E;
            }
        }
        return current_min;
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