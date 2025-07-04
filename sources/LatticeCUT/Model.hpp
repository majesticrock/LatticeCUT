#pragma once
#include "GlobalDefinitions.hpp"
#include <mrock/utility/InputFileReader.hpp>

namespace LatticeCUT {
    struct Model {
        l_float phonon_interaction;
        l_float local_interaction;
        std::vector<l_float> density_of_states;

        Model(mrock::utility::InputFileReader& input);
    };
}