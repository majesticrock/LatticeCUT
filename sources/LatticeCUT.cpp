#include <nlohmann/json.hpp>
#include <mrock/info.h>
#include <mrock/utility/OutputConvenience.hpp>
#include <mrock/utility/info_to_json.hpp>

#include "LatticeCUT/ModeHelper.hpp"
#include "../build_header/info.h"

using namespace LatticeCUT;
const std::string BASE_FOLDER = "../../data/lattice_cut/";

int main(int argc, char** argv) {
    if (argc < 2) {
		std::cerr << "Invalid number of arguments: Use <path_to_executable> <configfile>" << std::endl;
		return -1;
	}
    mrock::utility::InputFileReader input(argv[1]);

    ModeHelper modes(input);
    const std::string output_folder = BASE_FOLDER + input.getString("output_folder") + "/" + modes.getModel().to_folder();
	std::filesystem::create_directories(output_folder);

    nlohmann::json comments = {
		{ "time", 				mrock::utility::time_stamp() },
        { "dos_name",           modes.getModel().dos_name },
		{ "g", 					modes.getModel().phonon_coupling },
		{ "U",	                modes.getModel().local_interaction },
		{ "band_width",			modes.getModel().selector.get_band_width() },
        { "E_F", 				modes.getModel().fermi_energy },
		{ "omega_D", 			modes.getModel().omega_debye },
        { "N",              	modes.getModel().N },
        { "Delta_max", 			modes.getModel().delta_max() },
	};

    nlohmann::json info_json = mrock::utility::generate_json<LatticeCUT::info>("lattice_cut_");
	info_json.update(mrock::utility::generate_json<mrock::info>("mrock_"));
	
	std::cout << "Saving data to " << output_folder << std::endl;
	mrock::utility::saveString(info_json.dump(4), output_folder + "metadata.json.gz");

    /*
	* Compute and output gap data
	*/
	nlohmann::json jDelta = {
		{ "energies",	modes.getModel().energies.get_abscissa() },
#ifndef UNIFORM_DISCRETIZATION
		{ "inner_min",  modes.getModel().energies.inner_min },
		{ "inner_max",  modes.getModel().energies.inner_max },
#endif
        { "dos",        modes.getModel().selector.get_raw_dos() },
		{ "Delta", 	    modes.getModel().Delta.as_vector() }
	};
	jDelta.merge_patch(comments);
	mrock::utility::saveString(jDelta.dump(4), output_folder + "gap.json.gz");

    auto resolvents = modes.compute_collective_modes(150);
	if (!resolvents.empty()) {
		nlohmann::json jResolvents = {
			{ "resolvents", resolvents },
			{ "continuum_boundaries", modes.continuum_boundaries() }
		};
		jResolvents.merge_patch(comments);
		mrock::utility::saveString(jResolvents.dump(4),output_folder + "resolvents.json.gz");
	}

    return 0;
}
