#include "ModeHelper.hpp"
#include <mrock/utility/Selfconsistency/IterativeSolver.hpp>
#include <mrock/utility/Selfconsistency/BroydenSolver.hpp>

#define ieom_diag model->density_of_states[k]
#define ieom_offdiag model->density_of_states[k] * model->density_of_states[l]

#ifdef _complex
#define __conj(z) std::conj(z)
#else
#define __conj(z) z
#endif

#define loop_offset 0//model->energies.lower_discretization
#define loop_end model->N//(model->energies.inner_discretization + 1)
#define loop_end_init input.getInt("N")//(input.getInt("N") * DOS::EnergyRanges::inner_fraction + 1)

namespace LatticeCUT {
    int ModeHelper::select_epsilon(mrock::symbolic_operators::Momentum const &momentum, int k, int l, int q) const
    {
        assert(momentum.size() == 1U && "There should not be an addition of momenta, so that a momentum can be easily associated with an energy!");
        assert(std::abs(momentum.front().factor) == 1 && "The momentum should either come with + or - and nothing more!");
        switch(static_cast<char>(momentum.front().name)) {
            case 'k':
                return k;
            case 'l':
                return l;
            case 'q':
                return q;
            default:
				throw std::runtime_error("Momentum not recognized! " + momentum.front().name);
        }
    }

    l_float ModeHelper::get_expectation_value(mrock::symbolic_operators::WickOperator const& op, int k) const
    {
        if (op.type == mrock::symbolic_operators::Number_Type) {
			return model->occupation_index(k);
		}
		else if (op.type == mrock::symbolic_operators::SC_Type) {
			if (op.is_daggered) {
				return __conj(model->sc_expectation_value_index(k));
			}
			return model->sc_expectation_value_index(k);
		}
		throw std::runtime_error("Expectation value not recognized!");
    }

	l_float ModeHelper::compute_phonon_sum(const mrock::symbolic_operators::WickTerm& term, int k) const
    {
        const int q_dependend = term.which_operator_depends_on('q');
		mrock::symbolic_operators::WickOperator const* const summed_op = &(term.operators[q_dependend]);
		mrock::symbolic_operators::WickOperator const* const other_op = term.is_bilinear() ? nullptr : &(term.operators[q_dependend == 0]);
		l_float value{};

		const l_float energy_k = model->energies.index_to_energy(k);
        for (int q = model->phonon_lower_bound(energy_k); q <= model->phonon_upper_bound(energy_k); ++q) {
            value += this->get_expectation_value(*summed_op, q) * model->density_of_states[q];
        }
		if (other_op) {
			value *= this->get_expectation_value(*other_op, k);
		}
		value *= static_cast<l_float>(term.multiplicity) * model->phonon_coupling;
		return value;
    }

	l_float ModeHelper::compute_local_sum(const mrock::symbolic_operators::WickTerm& term, int k) const
    {
        const int q_dependend = term.which_operator_depends_on('q');
		mrock::symbolic_operators::WickOperator const* const summed_op = &(term.operators[q_dependend]);
		mrock::symbolic_operators::WickOperator const* const other_op = term.is_bilinear() ? nullptr : &(term.operators[q_dependend == 0]);
		l_float value{};

        for (int q = 0; q < model->N; ++q) {
            value += this->get_expectation_value(*summed_op, q) * model->density_of_states[q];
        }
		if (other_op) {
			value *= this->get_expectation_value(*other_op, k);
		}
		value *= static_cast<l_float>(term.multiplicity) * model->local_interaction_energy_units;
		return value;
    }

	void ModeHelper::createStartingStates()
    {
		starting_states.push_back({ _parent::Vector::Zero(antihermitian_discretization), _parent::Vector::Zero(hermitian_discretization), "SC" });
		int __b, __e;
		if (investigated_operator == InvestigatedOperator::Full) {
			__b = 0;
			__e = loop_end;
		}
		else if (investigated_operator == InvestigatedOperator::NearZero) {
			const int tenth = model->N / 10;
			__b = model->N / 2 - tenth;
			__e = model->N / 2 + tenth;
		}
		else if (investigated_operator == InvestigatedOperator::NearDoublePeaks) {
			const int tenth = model->N / 10;
			__b = model->N / 4 - tenth;
			__e = model->N / 4 + tenth;
		}
		else {
			throw std::invalid_argument("investigated_operator not recognized!");
		}
		for (int k = __b; k < __e; ++k) {
			starting_states[0][0](k) = 1.;
			starting_states[0][1](k) = 1.;
		}

		if (investigated_operator == InvestigatedOperator::NearDoublePeaks) {
			const int tenth = model->N / 10;
			__b = 3 * model->N / 4 - tenth;
			__e = 3 * model->N / 4 + tenth;

			for (int k = __b; k < __e; ++k) {
				starting_states[0][0](k) = 1.;
				starting_states[0][1](k) = 1.;
			}
		}
    }

	void ModeHelper::fillMatrices()
    {
		K_plus.setZero(hermitian_discretization, hermitian_discretization);
		K_minus.setZero(antihermitian_discretization, antihermitian_discretization);
		L.setZero(hermitian_discretization, antihermitian_discretization);

#pragma omp parallel for
		for (int i = 0; i < number_of_basis_terms; ++i)
		{
			for (int j = 0; j < number_of_basis_terms; ++j)
			{
				// Ignore the offdiagonal blocks as they are 0
				if ((i < hermitian_size && j < hermitian_size) || (j >= hermitian_size && i >= hermitian_size)) {
					fill_block_M(i, j);
				}
				// N only contains offdiagonal blocks
				else if (i < hermitian_size && j >= hermitian_size) {
					fill_block_N(i, j);
				}
			}
		}
		std::cout << "||K_+ - K_+^+|| = " << (K_plus - K_plus.adjoint()).norm() << std::endl;
		std::cout << "||K_- - K_-^+|| = " << (K_minus - K_minus.adjoint()).norm() << std::endl;
    }

	void ModeHelper::fill_M()
    {
        throw std::runtime_error("fill_M not implemented!");
    }
    
	void ModeHelper::fill_block_M(int i, int j)
    {
        for (int k = 0; k < loop_end; ++k) {
			l_float diag_buffer{};
			for (const auto& term : wicks.M[number_of_basis_terms * j + i]) {
				if (!term.delta_momenta.empty()) {
					if (term.sums.momenta.empty()) {
						if (term.coefficients.front().name == "g" || term.coefficients.front().name == "U") {
							// These kinds of terms scale as 1/V -> 0
							continue;
						}
					}
					diag_buffer += computeTerm(term, k + loop_offset, k + loop_offset);
				}
				else {
					for (int l = 0; l < loop_end; ++l) {
						M(i * loop_end + k, j * loop_end + l) += ieom_offdiag * computeTerm(term, k + loop_offset, l + loop_offset);
					}
				}
			}
			M(i * loop_end + k, j * loop_end + k) += ieom_diag * diag_buffer;
		}
    }

	void ModeHelper::fill_block_N(int i, int j)
    {
        for (int k = 0; k < loop_end; ++k) {
			for (const auto& term : wicks.N[number_of_basis_terms * j + i]) {
				if (!term.delta_momenta.empty()) {
					// only k=l and k=-l should occur. Additionally, only the magntitude should matter
					N(i * loop_end + k, j * loop_end + k) += computeTerm(term, k + loop_offset, k + loop_offset);
				}
				else {
					throw std::runtime_error("Offdiagonal term in N!");
				}
			}
			N(i * loop_end + k, j * loop_end + k) *= ieom_diag;
		}
    }

	l_float ModeHelper::computeTerm(const mrock::symbolic_operators::WickTerm& term, int k, int l) const
    {
        if (term.sums.momenta.empty()) {
			l_float value{ static_cast<l_float>(term.multiplicity) };

			if (!term.coefficients.empty()) {
				const mrock::symbolic_operators::Coefficient* coeff_ptr = &term.coefficients.front();
				if (coeff_ptr->momenta.size() == 2U) {
					value *= model->compute_coefficient(*coeff_ptr, select_epsilon(coeff_ptr->momenta[0], k, l), select_epsilon(coeff_ptr->momenta[1], k, l));
				}
				else if (coeff_ptr->momenta.size() == 1U || coeff_ptr->momenta.empty()) {
					value *= model->compute_coefficient(*coeff_ptr, k, l);
				}
				else {
					throw std::runtime_error("Number of momenta of coefficient is not handled! " + std::to_string(coeff_ptr->momenta.size()));
				}
			}

			if (term.operators.empty()) return value;
			for (const auto& op : term.operators) {
				value *= this->get_expectation_value(op, select_epsilon(op.momentum, k, l));
			}
			return value * static_cast<l_float>(term.sums.spins.size() + 1U);
		}
		assert(term.coefficients.size() == 1U);

		if (term.coefficients.front().name == "g")
		{
			return compute_phonon_sum(term, k);
		}
		else if (term.coefficients.front().name == "U") {
			return compute_local_sum(term, k);
		}
		throw std::runtime_error("Something went wrong while computing terms...");
    }

	std::vector<l_float> ModeHelper::continuum_boundaries() const
    {
        l_float prev_min { std::numeric_limits<l_float>::max() };
		l_float min { std::numeric_limits<l_float>::max() };
        for(int k = 0; k < model->N; ++k) {
            if (model->quasiparticle_energy_index(k) < prev_min) {
                prev_min = model->quasiparticle_energy_index(k);
                min = model->single_particle_energy(k);
            }
        }
        std::cout << "Found minimum at epsilon=" << min 
            << "   E(min)=" << prev_min
            << "   E(0)=" << model->quasiparticle_energy_index(model->energies.energy_to_index(model->fermi_energy)) << std::endl;
        return { 2 * prev_min,
			2 * std::max(model->quasiparticle_energy_index(loop_offset), model->quasiparticle_energy_index(loop_end - 1)) };
    }

	ModeHelper::ModeHelper(mrock::utility::InputFileReader& input)
        : _parent(this, SQRT_PRECISION, 
            static_cast<int>(loop_end_init * hermitian_size    ),
            static_cast<int>(loop_end_init * antihermitian_size), 
            false, false),
        model(std::make_unique<DOSModel>(input)),
		hermitian_discretization{loop_end * hermitian_size },
		antihermitian_discretization{loop_end * antihermitian_size },
		total_matrix_size{loop_end * number_of_basis_terms },
		investigated_operator{ static_cast<InvestigatedOperator>(input.getInt("investigated_operator")) }
	{
		std::cout << "Working on " << model->info() << std::endl;
		wicks.load("../commutators/lattice_cut/", true, number_of_basis_terms, 0);

#ifdef _iterative_selfconsistency
		auto solver = mrock::utility::Selfconsistency::make_iterative<l_float>(model.get(), &model->Delta);
#else
		auto solver = mrock::utility::Selfconsistency::make_broyden<l_float>(model.get(), &model->Delta, 200);
#endif
		solver.compute(true);
		std::cout << "\n########################################\n  ---  Delta_max = " << model->delta_max() << "  ---\n########################################\n" << std::endl;
	}
}