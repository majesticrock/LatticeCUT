#include <nlohmann/json.hpp>
#include <mrock/info.h>
#include <mrock/utility/OutputConvenience.hpp>
#include <mrock/utility/info_to_json.hpp>
#include <mrock/utility/FunctionTime.hpp>

#include "LatticeCUT/T_C.hpp"
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

	// Compute T_C
	if (input.getBool("T_C")) {
		std::cout << "Starting T_C computation..." << std::endl;
		T_C tc(input);
		mrock::utility::member_function_time_ms(tc, &T_C::compute);

		nlohmann::json comments = {
			{ "time", 				   	mrock::utility::time_stamp() },
	        { "dos_name",              	tc.model.dos_name },
			{ "g", 					   	tc.model.phonon_coupling_in },
			{ "U",	                   	tc.model.local_interaction },
	        { "E_F", 				   	tc.model.fermi_energy },
			{ "omega_D", 			   	tc.model.omega_debye_in },
	        { "N",              	   	tc.model.N },
		    { "filling_at_zero_temp",	tc.model.filling_at_zero_temp },
		};
		
		nlohmann::json info_json = mrock::utility::generate_json<LatticeCUT::info>("lattice_cut_");
		info_json.update(mrock::utility::generate_json<mrock::info>("mrock_"));

		const std::string output_folder = BASE_FOLDER + input.getString("output_folder") + "/" + tc.to_folder();
		std::filesystem::create_directories(output_folder);
		std::cout << "Saving data to " << output_folder << std::endl;
		mrock::utility::saveString(info_json.dump(4), output_folder + "metadata.json.gz");

		nlohmann::json jT_C = {
			{ "temperatures",			tc.temperatures },
			{ "max_gaps", 				tc.max_gaps },
			{ "true_gaps", 				tc.true_gaps },
			{ "gaps_at_ef", 			tc.gaps_at_ef },
			{ "chemical_potentials",	tc.chemical_potentials }
		};
		// All gaps are seperated because the data is rarely needed and takes a lot of time loading
		nlohmann::json jAllGaps = { { "finite_gaps",		tc.finite_gaps } };

		jT_C.merge_patch(comments);
		mrock::utility::saveString(jT_C.dump(4), output_folder + "T_C.json.gz");

		jAllGaps.merge_patch(comments);
		mrock::utility::saveString(jAllGaps.dump(4), output_folder + "all_gaps.json.gz");
		std::cout << "T_C computation finished." << std::endl;
	}

	// Study collective excitations
	if (input.getBool("collective_modes")) {
	    ModeHelper modes(input);
	    
	    nlohmann::json comments = {
			{ "time", 				   	mrock::utility::time_stamp() },
	        { "dos_name",              	modes.getModel().dos_name },
			{ "g", 					   	modes.getModel().phonon_coupling_in },
			{ "U",	                   	modes.getModel().local_interaction },
	        { "E_F", 				   	modes.getModel().fermi_energy },
			{ "omega_D", 			   	modes.getModel().omega_debye_in },
	        { "N",              	   	modes.getModel().N },
	        { "Delta_max", 			   	modes.getModel().delta_max() },
			{ "beta", 				   	modes.getModel().beta },
			{ "investigated_operator", 	static_cast<int>(modes.investigated_operator) },
			{ "filling_at_zero_temp",	modes.getModel().filling_at_zero_temp },
			{ "chemical_potential", 	modes.getModel().chemical_potential }
		};

	    nlohmann::json info_json = mrock::utility::generate_json<LatticeCUT::info>("lattice_cut_");
		info_json.update(mrock::utility::generate_json<mrock::info>("mrock_"));

		std::string mode_dir;
		if (input.getInt("investigated_operator") == 1) {
			mode_dir = "confined/";
		}
		else if (input.getInt("investigated_operator") == 2) {
			mode_dir = "double_confined/";
		}
		else {
			mode_dir = "/";
		}
		const std::string output_folder = BASE_FOLDER + input.getString("output_folder") +
			mode_dir + modes.getModel().to_folder();
		std::filesystem::create_directories(output_folder);
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
			{ "Delta", 	    modes.getModel().Delta.as_vector(modes.getModel().N) }
		};
		jDelta.merge_patch(comments);
		mrock::utility::saveString(jDelta.dump(4), output_folder + "gap.json.gz");
		std::cout << "Order parameter data saved!" << std::endl;

#ifdef LATTICE_CUT_RESIDUALS
	    auto [resolvents, residual_infos] = modes.compute_collective_modes_with_residuals(400);
#else
		auto resolvents = modes.compute_collective_modes(400);
#endif
		if (!resolvents.empty()) {
			nlohmann::json jResolvents = {
				{ "resolvents", resolvents },
				{ "continuum_boundaries", modes.continuum_boundaries() }
			};
			jResolvents.merge_patch(comments);
			mrock::utility::saveString(jResolvents.dump(4), output_folder + "resolvents.json.gz");
#ifdef LATTICE_CUT_RESIDUALS
			nlohmann::json jResiduals = {
				{ "phase", residual_infos.front() },
				{ "amplitude", residual_infos.back() }
			};
			jResiduals.merge_patch(comments);
			mrock::utility::saveString(jResiduals.dump(4), output_folder + "residuals.json.gz");
		}
#endif
		/* auto [phase_data, amplitude_data] = modes.full_diagonalization();
		nlohmann::json jFullDiag = {
			{ "phase", phase_data },
			{ "amplitude", amplitude_data },
			{ "continuum_boundaries", modes.continuum_boundaries() }
		};
		jFullDiag.merge_patch(comments);
		mrock::utility::saveString(jFullDiag.dump(4), output_folder + "full_diagonalization.json.gz"); */
	}

    return 0;
}
