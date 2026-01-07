#pragma once
#include "GlobalDefinitions.hpp"
#include "ModelAttributes.hpp"
#include "../DOS/Selector.hpp"
#include "../DOS/EnergyRanges.hpp"

#include <string>
#include <map>

#include <mrock/utility/InputFileReader.hpp>
#include <mrock/symbolic_operators/WickTerm.hpp>

namespace LatticeCUT {
    struct DOSModel {
        typedef Eigen::VectorXd ParameterVector;

        const l_float phonon_coupling_in; ///< g_in
        const l_float local_interaction; ///< in units of W
        const l_float fermi_energy; ///< in units of W
        const l_float omega_debye_in; ///< fraction of the bandwidth in which the phonon_coupling is active
        const int N; ///< Number of discretization points
        
        DOS::Selector selector;
        const std::string dos_name;
        const std::vector<l_float>& density_of_states;
        const DOS::EnergyRanges& energies;
        const l_float omega_debye; ///< used for calculations
        const l_float phonon_coupling; ///< g_in / \int_(E_F-omega_D)^(E_F+omega_D) rho(E) dE
        const l_float local_interaction_energy_units; ///< g_in / \int_(E_F-omega_D)^(E_F+omega_D) rho(E) dE
        l_float beta; ///< inverse temperature in units of W
        l_float chemical_potential; ///< chemical potential in units of W
        l_float filling_at_zero_temp; ///< filling at zero temperature

        ModelAttributes<l_float> Delta;

        DOSModel(mrock::utility::InputFileReader& input);

        void iteration_step(const ParameterVector& initial_values, ParameterVector& result);
        
        l_float sc_expectation_value_index(int k) const;
        l_float occupation_index(int k) const;

        inline l_float single_particle_energy(int k) const {
            return energies.index_to_energy(k) - chemical_potential;
        }

        inline l_float quasiparticle_energy_index(int k) const
        {
            return sqrt(single_particle_energy(k) * single_particle_energy(k) + std::norm(Delta[k]));
        }

        inline int phonon_lower_bound(l_float eps) const {
            return std::max(energies.energy_to_index(eps - omega_debye), int{});
        }
        inline int phonon_upper_bound(l_float eps) const {
            return std::min(energies.energy_to_index(eps + omega_debye), N - 1);
        }

        l_float compute_filling(const l_float mu) const;
        
        l_float compute_coefficient(mrock::symbolic_operators::Coefficient const& coeff, int first, int second) const;

        const std::map<mrock::symbolic_operators::OperatorType, std::vector<l_float>>& get_expectation_values() const;

        l_float delta_max() const;

        l_float true_gap() const;

        std::string info() const;
        std::string to_folder() const;

    private:
		mutable std::map<mrock::symbolic_operators::OperatorType, std::vector<l_float>> _expecs;
        const bool guaranteed_E_F{};
    };
}