#include "T_C.hpp"
#include <mrock/utility/Selfconsistency/IterativeSolver.hpp>
#include <mrock/utility/Selfconsistency/BroydenSolver.hpp>
#include <mrock/utility/better_to_string.hpp>

constexpr double TARGET_DT = 1e-4;
constexpr double INITIAL_DT = 0.01;

namespace LatticeCUT {
    T_C::T_C(mrock::utility::InputFileReader &input)
        : model(DOSModel(input))
    { 
        temperatures.reserve(200);
        finite_gaps.reserve(200);
    }

    void T_C::compute()
    {
#ifdef _iterative_selfconsistency
		auto solver = mrock::utility::Selfconsistency::make_iterative<l_float>(&model, &(model.Delta));
#else
		auto solver = mrock::utility::Selfconsistency::make_broyden<l_float>(&model, &(model.Delta), 200);
#endif
        l_float T{};
        l_float current_dT{INITIAL_DT};
        l_float delta_max{};
        model.beta = -1.;
        do {
            std::cout << "Working... T=" << T << std::endl;
            model.Delta.converged = false;
            solver.compute(false);
            delta_max = model.delta_max();
            if (!model.Delta.converged) {
                std::cerr << "Self-consistency not achieved while computing T_C! at beta=" << model.beta << std::endl;
            }
            else if (!is_zero(delta_max)) {
                temperatures.push_back(T);
                finite_gaps.push_back(model.Delta.as_vector());
                T += current_dT;
            }
            else {
                model.Delta.fill_with(finite_gaps.back());
                T -= current_dT;
                current_dT /= 5.;
                T += current_dT;
            }
            model.beta = 1. / T;
        } while((!is_zero(delta_max) || current_dT >= TARGET_DT) && model.Delta.converged);

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
