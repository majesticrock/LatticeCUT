#pragma once

#include "DOSModel.hpp"
#include "GlobalDefinitions.hpp"

#include <mrock/utility/InputFileReader.hpp>

#include <string>
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
}  // namespace LatticeCUT