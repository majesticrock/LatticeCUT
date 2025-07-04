#pragma once
#include "GlobalDefinitions.hpp"
#include "ModelAttributes.hpp"

#include <string>

#include <mrock/utility/InputFileReader.hpp>


namespace LatticeCUT {
    struct DOSModel {
        typedef Eigen::VectorXd ParameterVector;

        l_float phonon_coupling;
        l_float local_interaction;
        l_float omega_debye;
        std::string dos_name;
        std::vector<l_float> density_of_states;

        ModelAttributes<l_float> Delta;

        DOSModel(mrock::utility::InputFileReader& input);

        void iteration_step(const ParameterVector& initial_values, ParameterVector& result);
        std::string info() const;

        l_float sc_expectation_value_index(int k) const;
        l_float occupation_index(int k) const;
    };
}