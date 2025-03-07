/*++
Copyright (c) 2006 Microsoft Corporation

Module Name:

    smt_params.cpp

Abstract:

    <abstract>

Author:

    Leonardo de Moura (leonardo) 2008-02-20.

Revision History:

--*/
#include "smt/params/smt_params.h"
#include "smt/params/smt_params_helper.hpp"
#include "util/gparams.h"

void smt_params::updt_local_params(params_ref const & _p) {
    smt_params_helper p(_p);
    m_auto_config = p.auto_config() && gparams::get_value("auto_config") == "true"; // auto-config is not scoped by smt in gparams.
    m_random_seed = p.random_seed();
    m_use_ls=p.use_ls();
    // sampling
    m_sampling=p.sampling();
    m_seed=p.seed();
    m_input_file=p.input_file().str();
    m_output_file=p.output_file().str();
    m_ls_time=p.ls_time();
    m_relevancy_lvl = p.relevancy();
    m_ematching   = p.ematching();
    m_induction   = p.induction();
    m_clause_proof = p.clause_proof();
    m_phase_selection = static_cast<phase_selection>(p.phase_selection());
    if (m_phase_selection > PS_THEORY) throw default_exception("illegal phase selection numeral");
    m_phase_caching_on = p.phase_caching_on();
    m_phase_caching_off = p.phase_caching_off();
    m_restart_strategy = static_cast<restart_strategy>(p.restart_strategy());
    if (m_restart_strategy > RS_ARITHMETIC) throw default_exception("illegal restart strategy numeral");
    m_restart_factor = p.restart_factor();
    m_case_split_strategy = static_cast<case_split_strategy>(p.case_split());
    m_theory_case_split = p.theory_case_split();
    m_theory_aware_branching = p.theory_aware_branching();
    m_delay_units = p.delay_units();
    m_delay_units_threshold = p.delay_units_threshold();
    m_preprocess = _p.get_bool("preprocess", true); // hidden parameter
    m_max_conflicts = p.max_conflicts();
    m_restart_max   = p.restart_max();
    m_cube_depth    = p.cube_depth();
    m_threads       = p.threads();
    m_threads_max_conflicts  = p.threads_max_conflicts();
    m_threads_cube_frequency = p.threads_cube_frequency();
    m_core_validate = p.core_validate();
    m_logic = _p.get_sym("logic", m_logic);
    m_string_solver = p.string_solver();
    validate_string_solver(m_string_solver);
    if (_p.get_bool("arith.greatest_error_pivot", false))
        m_arith_pivot_strategy = arith_pivot_strategy::ARITH_PIVOT_GREATEST_ERROR;
    else if (_p.get_bool("arith.least_error_pivot", false))
        m_arith_pivot_strategy = arith_pivot_strategy::ARITH_PIVOT_LEAST_ERROR;
    theory_array_params::updt_params(_p);
    m_dump_benchmarks = false;
    m_dump_min_time = 0.5;
    m_dump_recheck = false;
}

void smt_params::updt_params(params_ref const & p) {
    preprocessor_params::updt_params(p);
    dyn_ack_params::updt_params(p);
    qi_params::updt_params(p);
    theory_arith_params::updt_params(p);
    theory_bv_params::updt_params(p);
    theory_pb_params::updt_params(p);
    // theory_array_params::updt_params(p);
    theory_datatype_params::updt_params(p);
    theory_str_params::updt_params(p);
    updt_local_params(p);
}

void smt_params::updt_params(context_params const & p) {
    m_auto_config    = p.m_auto_config;
    m_model          = p.m_model;
}

#define DISPLAY_PARAM(X) out << #X"=" << X << std::endl;

void smt_params::display(std::ostream & out) const {
    preprocessor_params::display(out);
    dyn_ack_params::display(out);
    qi_params::display(out);
    theory_arith_params::display(out);
    theory_array_params::display(out);
    theory_bv_params::display(out);
    theory_pb_params::display(out);
    theory_datatype_params::display(out);
    theory_str_params::display(out);

    DISPLAY_PARAM(m_display_proof);
    DISPLAY_PARAM(m_display_dot_proof);
    DISPLAY_PARAM(m_display_unsat_core);
    DISPLAY_PARAM(m_check_proof);
    DISPLAY_PARAM(m_eq_propagation);
    DISPLAY_PARAM(m_binary_clause_opt);
    DISPLAY_PARAM(m_relevancy_lvl);
    DISPLAY_PARAM(m_relevancy_lemma);
    DISPLAY_PARAM(m_random_seed);
    DISPLAY_PARAM(m_random_var_freq);
    DISPLAY_PARAM(m_inv_decay);
    DISPLAY_PARAM(m_clause_decay);
    DISPLAY_PARAM(m_random_initial_activity);
    DISPLAY_PARAM(m_phase_selection);
    DISPLAY_PARAM(m_phase_caching_on);
    DISPLAY_PARAM(m_phase_caching_off);
    DISPLAY_PARAM(m_minimize_lemmas);
    DISPLAY_PARAM(m_max_conflicts);
    DISPLAY_PARAM(m_cube_depth);
    DISPLAY_PARAM(m_threads);
    DISPLAY_PARAM(m_threads_max_conflicts);
    DISPLAY_PARAM(m_threads_cube_frequency);
    DISPLAY_PARAM(m_simplify_clauses);
    DISPLAY_PARAM(m_tick);
    DISPLAY_PARAM(m_display_features);
    DISPLAY_PARAM(m_new_core2th_eq);
    DISPLAY_PARAM(m_ematching);
    DISPLAY_PARAM(m_induction);
    DISPLAY_PARAM(m_clause_proof);

    DISPLAY_PARAM(m_case_split_strategy);
    DISPLAY_PARAM(m_rel_case_split_order);
    DISPLAY_PARAM(m_lookahead_diseq);

    DISPLAY_PARAM(m_delay_units);
    DISPLAY_PARAM(m_delay_units_threshold);

    DISPLAY_PARAM(m_theory_resolve);

    DISPLAY_PARAM(m_restart_strategy);
    DISPLAY_PARAM(m_restart_initial);
    DISPLAY_PARAM(m_restart_factor);
    DISPLAY_PARAM(m_restart_adaptive);
    DISPLAY_PARAM(m_agility_factor);
    DISPLAY_PARAM(m_restart_agility_threshold);

    DISPLAY_PARAM(m_lemma_gc_strategy);
    DISPLAY_PARAM(m_lemma_gc_half);
    DISPLAY_PARAM(m_recent_lemmas_size);
    DISPLAY_PARAM(m_lemma_gc_initial);
    DISPLAY_PARAM(m_lemma_gc_factor);
    DISPLAY_PARAM(m_new_old_ratio);
    DISPLAY_PARAM(m_new_clause_activity);
    DISPLAY_PARAM(m_old_clause_activity);
    DISPLAY_PARAM(m_new_clause_relevancy);
    DISPLAY_PARAM(m_old_clause_relevancy);
    DISPLAY_PARAM(m_inv_clause_decay);

    DISPLAY_PARAM(m_smtlib_dump_lemmas);
    DISPLAY_PARAM(m_logic);
    DISPLAY_PARAM(m_string_solver);

    DISPLAY_PARAM(m_profile_res_sub);
    DISPLAY_PARAM(m_display_bool_var2expr);
    DISPLAY_PARAM(m_display_ll_bool_var2expr);

    DISPLAY_PARAM(m_model);
    DISPLAY_PARAM(m_model_on_timeout);
    DISPLAY_PARAM(m_model_on_final_check);

    DISPLAY_PARAM(m_progress_sampling_freq);

    DISPLAY_PARAM(m_core_validate);

    DISPLAY_PARAM(m_preprocess);
    DISPLAY_PARAM(m_user_theory_preprocess_axioms);
    DISPLAY_PARAM(m_user_theory_persist_axioms);
    DISPLAY_PARAM(m_at_labels_cex);
    DISPLAY_PARAM(m_check_at_labels);
    DISPLAY_PARAM(m_dump_goal_as_smt);
    DISPLAY_PARAM(m_auto_config);
    DISPLAY_PARAM(m_use_ls);
    DISPLAY_PARAM(m_sampling);
    DISPLAY_PARAM(m_seed);
    DISPLAY_PARAM(m_ls_time);
    DISPLAY_PARAM(m_input_file);
    DISPLAY_PARAM(m_output_file);
}

void smt_params::validate_string_solver(symbol const& s) const {
    if (s == "z3str3" || s == "seq" || s == "empty" || s == "auto" || s == "none")
        return;
    throw default_exception("Invalid string solver value. Legal values are z3str3, seq, empty, auto, none");
}
