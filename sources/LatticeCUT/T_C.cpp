#include "T_C.hpp"
#include <mrock/utility/Selfconsistency/IterativeSolver.hpp>
#include <mrock/utility/Selfconsistency/BroydenSolver.hpp>
#include <mrock/utility/better_to_string.hpp>
#include <algorithm>

constexpr double TARGET_DT = 1e-5;
constexpr double INITIAL_DT = 5e-3;
constexpr double ZERO_EPS = 1e-10;

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

        model.beta = -1.;
        solver.compute(false);
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

        temperatures.emplace_back(T);
        finite_gaps.emplace_back(model.Delta.as_vector());

        do {
            T += current_dT;
            model.beta = 1. / T;

            if (std::find_if(temperatures.begin(), temperatures.end(), [&T](const l_float Tvec) -> bool { return is_zero(Tvec-T); }) != temperatures.end()) {
                continue;
            }
            std::cout << "Working... T=" << T << std::endl;
            model.Delta.converged = false;
            solver.compute(false);
            delta_max = model.delta_max();
            std::cout << "\t\tDelta_max=" << delta_max << std::endl;

            if (!model.Delta.converged) {
                std::cerr << "Self-consistency not achieved while computing T_C! at beta=" << model.beta << std::endl;
                break;
            }

            if (std::abs(delta_max) > ZERO_EPS) { // the is_zero function is sometimes too precise
                temperatures.emplace_back(T);
                finite_gaps.emplace_back(model.Delta.as_vector());
            }
            if (delta_max < 0.85 * last_delta && current_dT >= TARGET_DT) {
                T -= current_dT;
                if (T < 0) T = 0.0; // circumvent rare case floating point arithmetic issues
                current_dT *= 0.2;
            }
            if (std::abs(delta_max) > ZERO_EPS) {
                last_delta = delta_max;
            }
            else {
                model.Delta.fill_with(finite_gaps.back());
            }

        } while(std::abs(delta_max) > ZERO_EPS || current_dT >= TARGET_DT);

        std::vector<size_t> indices(temperatures.size());
        std::iota(indices.begin(), indices.end(), size_t{});
        std::sort(indices.begin(), indices.end(), [&](size_t A, size_t B) -> bool { return temperatures[A] < temperatures[B]; });
        permute(temperatures, indices);
        permute(finite_gaps, indices);

        max_gaps.resize(finite_gaps.size());
        for (size_t i = 0U; i < max_gaps.size(); ++i) {
            max_gaps[i] = std::abs(std::ranges::max(finite_gaps[i], [](const l_float A, const l_float B) -> bool { 
                    return std::abs(A) < std::abs(B);
                }));
        }

        std::cout << "Finished T_C computation at T_C = " << T << " (beta=" << model.beta << ") in " << temperatures.size() << "iterations." << std::endl;
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
