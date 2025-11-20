#include "T_C.hpp"
#include <mrock/utility/Selfconsistency/IterativeSolver.hpp>
#include <mrock/utility/Selfconsistency/BroydenSolver.hpp>
#include <mrock/utility/better_to_string.hpp>
#include <algorithm>

constexpr double TARGET_DT = 1e-5;
constexpr double INITIAL_DT = 5e-3;
constexpr double ZERO_EPS = 1e-10;
constexpr double DELTA_F_EPS = 1e-4;

namespace LatticeCUT {
    template<class T>
    void permute(std::vector<T>& input, const std::vector<size_t>& indices)
    {
        std::vector<bool> done(input.size());
        for (size_t i = 0U; i < input.size(); ++i) {
            if (done[i]) continue;

            done[i] = true;
            size_t prev_j = i;
            size_t j = indices[i];
            while(i != j) {
                std::swap(input[prev_j], input[j]);
                done[j] = true;
                prev_j = j;
                j = indices[j];
            }
        }
    };

    T_C::T_C(mrock::utility::InputFileReader &input)
        : model(DOSModel(input))
    { 
        temperatures.reserve(500);
        finite_gaps.reserve(500);
    }

    void T_C::compute()
    {
#ifdef _iterative_selfconsistency
		auto solver = mrock::utility::Selfconsistency::make_iterative<l_float>(&model, &(model.Delta));
#else
		auto solver = mrock::utility::Selfconsistency::make_broyden<l_float>(&model, &(model.Delta), 300);
#endif
        l_float T{};
        l_float current_dT{INITIAL_DT};
        l_float delta_max{};
        l_float last_delta{};
        l_float last_delta_F{};
        bool did_last_converge{true};

        model.beta = is_zero(T) ? -1. : 1. / T;
        solver.compute(false);
        int index_at_ef = static_cast<int>(0.5 * model.N * (model.chemical_potential + 1));
        // the delta_max function uses the absolute value
        delta_max = model.delta_max();
        // use U(1) symmetry to unifiy delta_max > 0
        if (delta_max < 0) {
            for (auto& d : model.Delta) {
                d *= -1;
            }
            delta_max *= -1;
        }
        last_delta = delta_max;
        last_delta_F = model.Delta[index_at_ef];

        std::cout << "Starting T_c computation at T=" << T << " (beta=" << model.beta << ")" << std::endl;
        std::cout << "\t\tDelta_max=" << delta_max << "\tDelta_F=" << model.Delta[index_at_ef] << std::endl;

        temperatures.emplace_back(T);
        finite_gaps.emplace_back(model.Delta.as_vector(model.N));
        gaps_at_ef.emplace_back(model.Delta[index_at_ef]);
        chemical_potentials.emplace_back(model.Delta[model.N]);
        true_gaps.emplace_back(model.true_gap());

        auto increase_dT = [&]() -> bool {
            if (current_dT >= 0.5 * INITIAL_DT) return false;
            if (is_zero(last_delta)) return false;
            if (std::abs((delta_max - last_delta) / last_delta) > 0.05) return false; 
            if (std::abs(last_delta_F) < DELTA_F_EPS * delta_max) {
                return std::abs(model.Delta[index_at_ef]) < DELTA_F_EPS * delta_max;
            }
            return (std::abs(model.Delta[index_at_ef] - last_delta_F) / last_delta_F) < 0.05; 
        };
        auto decrease_dT = [&]() -> bool {
            if (current_dT < TARGET_DT) return false;
            if (delta_max < 0.85 * last_delta) return true;
            if (std::abs(model.Delta[index_at_ef]) < 0.85 * std::abs(last_delta_F)) {
                return std::abs(model.Delta[index_at_ef]) > DELTA_F_EPS * delta_max;
            }
            return false;
        };
        std::cout << std::setprecision(6) << std::endl;

        do {
            T += current_dT;
            model.beta = 1. / T;

            const auto found_iter = std::find_if(temperatures.begin(), temperatures.end(), [&T](const l_float Tvec) -> bool { 
                    return std::abs(Tvec-T) < TARGET_DT; }
                );
            if (found_iter != temperatures.end()) {
                const size_t index = std::distance(temperatures.begin(), found_iter);

                last_delta = std::abs(std::ranges::max(finite_gaps[index], [](const double lhs, const double rhs){ return std::abs(lhs) < std::abs(rhs);}));
                last_delta_F = finite_gaps[index][index_at_ef];
                continue;
            }
            std::cout << "Working... T=" << T << "    current dT=" << current_dT << std::endl;
            model.Delta.converged = false;
            solver.compute(false);
            index_at_ef = static_cast<int>(0.5 * model.N * (model.chemical_potential + 1));
            delta_max = model.delta_max();
            

            if (!model.Delta.converged) {
                std::cerr << "Self-consistency not achieved while computing T_C! Retrying... at beta=" << model.beta << std::endl;
                solver.compute(false);
                index_at_ef = static_cast<int>(0.5 * model.N * (model.chemical_potential + 1));
                delta_max = model.delta_max();
                if (!model.Delta.converged) {
		        	std::cerr << "No convergence even after retry. Skipping data point." << std::endl;
                    if (!did_last_converge) {
                        break;
                    }
                    else {
                        did_last_converge = false;
                        continue;
                    }
		        }
            }
            std::cout << "\t\tDelta_max=" << delta_max 
                << "\tDelta_F=" << model.Delta[index_at_ef] 
                << "\tmu=" << model.chemical_potential << std::endl;

            if (decrease_dT()) {
                if (std::abs(delta_max) > ZERO_EPS) { // the is_zero function is sometimes too precise
                    temperatures.emplace_back(T);
                    finite_gaps.emplace_back(model.Delta.as_vector(model.N));
                    gaps_at_ef.emplace_back(model.Delta[index_at_ef]);
                    chemical_potentials.emplace_back(model.Delta[model.N]);
                    true_gaps.emplace_back(model.true_gap());
                }

                if ((std::abs(last_delta_F) > ZERO_EPS) && (std::abs(finite_gaps.back()[index_at_ef]) < ZERO_EPS)) {
                    model.Delta.fill_with(*(finite_gaps.end() - 2));
                    model.Delta[model.N] = *(chemical_potentials.end() - 2);
                    model.chemical_potential = model.Delta[model.N];
                }
                else {
                    model.Delta.fill_with(finite_gaps.back());
                    model.Delta[model.N] = chemical_potentials.back();
                    model.chemical_potential = model.Delta[model.N];
                }
                
                T -= current_dT;
                if (T < 0) T = 0.0; // circumvent rare case floating point arithmetic issues
                current_dT *= 0.2;
                // As we say in German: man muss ja nicht gleich uebertreiben...
                if (current_dT < 0.5 * TARGET_DT) current_dT = 0.5 * TARGET_DT;
            }
            else {
                if (increase_dT()) {
                    current_dT *= 2.0;
                }
                last_delta = delta_max;
                last_delta_F = model.Delta[index_at_ef];

                if (std::abs(delta_max) > ZERO_EPS) { // the is_zero function is sometimes too precise
                    temperatures.emplace_back(T);
                    finite_gaps.emplace_back(model.Delta.as_vector(model.N));
                    gaps_at_ef.emplace_back(model.Delta[index_at_ef]);
                    chemical_potentials.emplace_back(model.Delta[model.N]);
                    true_gaps.emplace_back(model.true_gap());
                }
                else {
                    model.Delta.fill_with(finite_gaps.back());
                    model.Delta[model.N] = chemical_potentials.back();
                    model.chemical_potential = model.Delta[model.N];
                }
                did_last_converge = true;
            }
        } while(std::abs(delta_max) > ZERO_EPS || current_dT >= TARGET_DT);

        std::vector<size_t> indices(temperatures.size());
        std::iota(indices.begin(), indices.end(), size_t{});
        std::sort(indices.begin(), indices.end(), [&](size_t A, size_t B) -> bool { return temperatures[A] < temperatures[B]; });
        permute(temperatures, indices);
        permute(finite_gaps, indices);
        permute(true_gaps, indices);
        permute(gaps_at_ef, indices);
        permute(chemical_potentials, indices);

        max_gaps.resize(finite_gaps.size());
        for (size_t i = 0U; i < max_gaps.size(); ++i) {
            max_gaps[i] = std::abs(std::ranges::max(finite_gaps[i], [](const l_float A, const l_float B) -> bool { 
                    return std::abs(A) < std::abs(B);
                }));
        }

        std::cout << "Finished T_C computation at T_C = " << T << " (beta=" << 1./(T > 0.0 ? T : -1.) << ") in " << temperatures.size() << "iterations." << std::endl;
    }

    std::string T_C::to_folder() const
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
        
        return "T_C/" + model.dos_name + "/N=" + std::to_string(model.N) + "/"
            + "g=" + improved_string(model.phonon_coupling_in) + "/"
            + "U=" + improved_string(model.local_interaction) + "/"
            + "E_F=" + improved_string(model.fermi_energy) + "/"
            + "omega_D=" + improved_string(model.omega_debye_in) + "/";
    }
}
