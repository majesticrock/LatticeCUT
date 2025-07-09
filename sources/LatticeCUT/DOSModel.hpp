#pragma once
#include "GlobalDefinitions.hpp"
#include "ModelAttributes.hpp"
#include "../DOS/Selector.hpp"

#include <string>
#include <map>

#include <mrock/utility/InputFileReader.hpp>
#include <mrock/symbolic_operators/WickTerm.hpp>

namespace LatticeCUT {
    struct DOSModel {
        typedef Eigen::VectorXd ParameterVector;

        const l_float phonon_coupling; ///< g_in
        const l_float local_interaction; ///< in units of the hopping constant
        const l_float fermi_energy; ///< in units of the hopping constant
        const l_float omega_debye_in; ///< fraction of discretization points (N) in which the phonon_coupling is active
        const int N; ///< Number of discretization points
        const int omega_debye; ///< in units of Delta epsilon := epislon_(i+1) - epsilon_i
        
        DOS::Selector selector;
        const std::string dos_name;
        const std::vector<l_float>& density_of_states;
        const l_float min_energy; ///< in units of the hopping constant
        const l_float delta_epsilon; ///< in units of the hopping constant
           
        ModelAttributes<l_float> Delta;

        DOSModel(mrock::utility::InputFileReader& input);

        void iteration_step(const ParameterVector& initial_values, ParameterVector& result);
        
        l_float sc_expectation_value_index(int k) const;
        l_float occupation_index(int k) const;

        inline l_float single_particle_energy(int k) const {
            return k * delta_epsilon - fermi_energy + min_energy;
        }

        inline l_float quasiparticle_energy_index(int k) const
        {
            return sqrt(single_particle_energy(k) * single_particle_energy(k) + std::norm(Delta[k]));
        }

        inline int phonon_lower_bound(int k) const {
            return std::max(k - omega_debye, int{});
        }
        inline int phonon_upper_bound(int k) const {
            return std::min(k + omega_debye, N - 1);
        }

        l_float compute_coefficient(mrock::symbolic_operators::Coefficient const& coeff, int first, int second) const;

        const std::map<mrock::symbolic_operators::OperatorType, std::vector<l_float>>& get_expectation_values() const;

        l_float delta_max() const;

        std::string info() const;
        std::string to_folder() const;

    private:
		mutable std::map<mrock::symbolic_operators::OperatorType, std::vector<l_float>> _expecs;
    };
}