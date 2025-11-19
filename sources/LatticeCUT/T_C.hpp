#pragma once

#include "GlobalDefinitions.hpp"
#include "DOSModel.hpp"

#include <mrock/utility/InputFileReader.hpp>
#include <vector>

namespace LatticeCUT {
    struct T_C {
        DOSModel model;
        std::vector<l_float> temperatures;
        std::vector<std::vector<l_float>> finite_gaps;
        std::vector<l_float> max_gaps;
        std::vector<l_float> true_gaps;
        std::vector<l_float> gaps_at_ef;
        std::vector<l_float> chemical_potentials;

        T_C(mrock::utility::InputFileReader& input);

        void compute();

        std::string to_folder() const;
    };
}