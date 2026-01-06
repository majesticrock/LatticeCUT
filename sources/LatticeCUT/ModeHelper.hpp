#pragma once
#include "GlobalDefinitions.hpp"
#include <mrock/symbolic_operators/TermLoader.hpp>
#include <mrock/symbolic_operators/WickTerm.hpp>
#include <mrock/utility/better_to_string.hpp>
#include "DOSModel.hpp"
#include <memory>

#ifndef _complex
#define _XP
#endif

#define IEOM_CLASS_TYPE mrock::utility::Numerics::iEoM::XPResolvent<ModeHelper, l_float, 10>

#include <mrock/utility/Numerics/iEoM/XPResolvent.hpp>

namespace LatticeCUT {
	enum class InvestigatedOperator { Full = 0, NearZero = 1, NearDoublePeaks = 2 };

	class ModeHelper : public IEOM_CLASS_TYPE
	{
	private:
		friend struct IEOM_CLASS_TYPE;
		using _parent = IEOM_CLASS_TYPE;

        int select_epsilon(mrock::symbolic_operators::Momentum const& momentum, int k, int l, int q=int{}) const;
        l_float get_expectation_value(mrock::symbolic_operators::WickOperator const& op, int k) const;

		l_float compute_phonon_sum(const mrock::symbolic_operators::WickTerm& term, int k) const;
		l_float compute_local_sum(const mrock::symbolic_operators::WickTerm& term, int k) const;
	protected:
		mrock::symbolic_operators::TermLoader wicks;
		//size_t TOTAL_BASIS{};
		constexpr static int hermitian_size = 2;
		constexpr static int antihermitian_size = 1;
		constexpr static int number_of_basis_terms = hermitian_size + antihermitian_size;

		std::unique_ptr<DOSModel> model;

		const int hermitian_discretization;
		const int antihermitian_discretization;
		const int total_matrix_size;		

		void createStartingStates();
		void fillMatrices();
		void fill_M();

		void fill_block_M(int i, int j);
		void fill_block_N(int i, int j);

		l_float computeTerm(const mrock::symbolic_operators::WickTerm& term, int k, int l) const;
	public:
		const InvestigatedOperator investigated_operator;
		std::vector<l_float> continuum_boundaries() const;

		inline DOSModel& getModel() {
			return *model;
		};
		inline const DOSModel& getModel() const {
			return *model;
		};
		ModeHelper(mrock::utility::InputFileReader& input);
	};
}

#undef IEOM_CLASS_TYPE