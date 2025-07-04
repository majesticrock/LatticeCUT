#include "DOSModel.hpp"

namespace LatticeCUT {
    DOSModel::DOSModel(mrock::utility::InputFileReader& input) 
        : phonon_coupling{input.getDouble("phonon_coupling")},
        local_interaction{input.getDouble("local_interaction")},
        omega_debye{input.getDouble("omega_debye")}
    { 
        // TODO: Fill DOS
    }

    void DOSModel::iteration_step(const ParameterVector& initial_values, ParameterVector& result)
    {

    }

    std::string DOSModel::info() const
    {
        return "DOSModel: g=" + std::to_string(phonon_coupling) + "\tU=" + std::to_string(local_interaction) + "\tDOS=" + dos_name;
    }
}