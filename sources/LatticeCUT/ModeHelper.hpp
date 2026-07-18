#pragma once
#include "DOSModel.hpp"
#include "GlobalDefinitions.hpp"

#include <mrock/iEoM/XPResolvent.hpp>
#include <mrock/symbolic_operators/TermLoader.hpp>
#include <mrock/symbolic_operators/WickOperator.hpp>
#include <mrock/symbolic_operators/WickTerm.hpp>
#include <mrock/utility/better_to_string.hpp>

#include <memory>
#include <vector>

namespace LatticeCUT {
enum class InvestigatedOperator { Full = 0, NearZero = 1, NearDoublePeaks = 2, Preset = 255 };

typedef mrock::iEoM::XPResolvent<l_float, 12, false> iEoM_algorithm;

class ModeHelper : public iEoM_algorithm {
private:
    constexpr static int hermitian_size = 2;
    constexpr static int antihermitian_size = 1;
    constexpr static int number_of_basis_terms = hermitian_size + antihermitian_size;

    mrock::symbolic_operators::TermLoader wicks;
    std::unique_ptr<DOSModel> model;

    const int hermitian_discretization;
    const int antihermitian_discretization;
    const int total_matrix_size;

    void create_starting_states();
    void fill_matrices();
    void fill_M();

    void fill_block_M(int i, int j);
    void fill_block_N(int i, int j);

    l_float computeTerm(const mrock::symbolic_operators::WickTerm& term, int k, int l) const;

    int select_epsilon(mrock::symbolic_operators::Momentum const& momentum, int k, int l, int q = int{}) const;
    l_float get_expectation_value(mrock::symbolic_operators::WickOperator const& op, int k) const;

    l_float compute_phonon_sum(const mrock::symbolic_operators::WickTerm& term, int k) const;
    l_float compute_local_sum(const mrock::symbolic_operators::WickTerm& term, int k) const;

public:
    const InvestigatedOperator investigated_operator;
    std::vector<l_float> continuum_boundaries() const;

    inline DOSModel& getModel() { return *model; };
    inline const DOSModel& getModel() const { return *model; };
    ModeHelper(mrock::utility::InputFileReader& input);
};
}  // namespace LatticeCUT

#undef iEoM_algorithm