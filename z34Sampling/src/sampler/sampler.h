#pragma once
#include <signal.h>
#include <string.h>
#include <sys/types.h>
#include <time.h>
#include <unistd.h>

#include <algorithm>
#include <chrono>
#include <cstdio>
#include <cstdlib>
#include <exception>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <random>
#include <set>
#include <stack>
#include <string>
#include <unordered_set>

#include "lia_Array.h"
#include "util/debug.h"
#include "util/trace.h"

// #define EQ2INEQ_BUILD_TABLE
// #define PROB_GUIDED
#define INTERVAL_MOVE
#define BUILD_OCCS_CLOSURES
// #define BIT_GUIDED
// #define CALC_DIVERSE_SAMPLE

// #define DEBUG
// #define VERBOSE
// #define PRINT_PAR_SAMPLES
// #define PRINT_INTERVAL
// #define PRINT_OCCS
// #define PRINT_UNCOV_BITS
// #define PRINT_BIT_GAIN

// #define PRINT_PROGRESS
// #define SHOW_PROGRESS_BAR

#define SAMPLER_TRACE(CODE) TRACE("sampler", CODE)
#define SAMPLER_CTRACE(COND, CODE) CTRACE("sampler", COND, CODE)

namespace sampler {

typedef std::vector<__int128_t> SampleBitSet;

struct Hash {
    std::size_t operator()(__int128_t x) const {
        return static_cast<std::size_t>(x ^ (x >> 64));
    }
};

const __int128_t max_int = __int128_t(INT64_MAX) * __int128_t(INT64_MAX);
const int64_t pos_inf_64 = INT64_MAX;
const int64_t neg_inf_64 = INT64_MIN;
const int32_t pos_inf_32 = INT32_MAX;
const int32_t neg_inf_32 = INT32_MIN;
const int16_t pos_inf_16 = INT16_MAX;
const int16_t neg_inf_16 = INT16_MIN;

struct lit {
    std::vector<int> pos_coff_var_idx;  // index of positive coefficient variables
    std::vector<__int128_t> pos_coff;   // coefficients of positive coefficient variables
    std::vector<int> neg_coff_var_idx;  // index of negative coefficient variables
    std::vector<__int128_t> neg_coff;   // coefficients of variables with negative coefficients
    __int128_t key;                     // constant term after conversion to standard form
    int lits_index;                     // literal index
    __int128_t delta;                   // for the lia lit, the current value of left side; for the boolean lit, the var's index in _resolution_vars
    bool is_equal = false;              // true means a-b-k==0, else a-b-k<=0
    int equal_pair = -1;
    bool is_lia_lit = false;  // true means this is a lia lit
};

struct variable {
    // literals[i] -- literal_clause[i] -- literal_coff[i] one-to-one correspondence
    std::vector<int> literals;             // literals[i]=l means the ith literal of the var is the pos(neg) of lth of _lits, it can be negative
                                           // duplicate elements may exist
    std::vector<int> literal_clause;       // literal_clause[i]=c means the ith literal containing the var is in cth clause.
                                           // there may be repetition, i.e. a clause containing repeated lia literal
    std::vector<__int128_t> literal_coff;  // literal_coff[i] denotes the coff of the var in corresponding literal, it can be negative

    std::vector<uint64_t> clause_idxs;  // clause index, implies which clauses the current variable is in
    std::string var_name;               // variable's name
    __int128_t low_bound = -max_int;    // lower bound
    __int128_t upper_bound = max_int;   // upper bound
    bool is_lia;                        // true means lia variable(int), else boolean variable
    bool is_delete = false;             // true means var is deleted
    __int128_t score;                   // if it is a bool var, then the score is calculated beforehand
    int up_bool = 0;                    // the bool value of variables deleted(1 true -1 false)
    bool is_in_equal = false;
    int64_t s_lower_bound = neg_inf_64;
    int64_t s_upper_bound = pos_inf_64;
    int occs = 0;
    bool init_to_zero = false;
};

struct clause {
    std::vector<int> literals;       // literals[i]=l means the ith literal of the clause if the pos(neg) of the _lits, it can be negative
    std::vector<int> lia_literals;   // linear integer arithmetic literals
    std::vector<int> bool_literals;  // boolean literals
    int weight = 1;                  // clause weight
    int sat_count;                   // the number of literals in the clause that have been satisfied ?
    __int128_t min_delta;            // a positive value, the distance from sat, delta for pos lit, 1-delta for neg lit
    int min_delta_lit_index;         // the lit index with the min_delta
    bool is_delete = false;          // true means var is deleted
};

class ls_sampler {
   private:
    // strategy
    bool bit_guided = false;

    // params
    unsigned _random_seed;
    std::string _input_path;
    std::string _output_path;
    int equal_cnt = 0;

    /* statistic */
    // variables
    uint64_t _num_vars;
    uint64_t _num_lia_vars;
    // literals
    uint64_t _num_lits;  // the number of literals read in from Z3
    int _num_lia_lits;   // the number of lia literals
    int _num_bool_lits;  // the number of boolean literals
    // clauses
    uint64_t _num_clauses;  // the number of clauses
    uint64_t _num_opt = 0;  // the number of vars in all literals, which is the max number of operations

    /* sampling */
    // out put
    std::ofstream samples_file;
    std::ofstream verify_file;

    /* internal data structure */
    // variables
    int _lia_var_idx_with_most_lits;
    std::vector<variable> _resolution_vars;  // or vars/if vars/ boolean vars
    std::vector<int> _lia_var_vec;           // lia vars index in _resolution_vars
    std::vector<int> _bool_var_vec;          // boolean vars index in _resolution_vars
    std::vector<variable> _vars;             // ??
    std::vector<variable> _tmp_vars;         // LIA vars
    // literals
    std::vector<lit> _lits;
    std::vector<int> _bound_lits;   // record the index of bounded lits
    Array* _lit_occur;              // the lit containing the lia var in one single clause
    std::vector<bool> _lit_appear;  // ??
    // clauses
    std::vector<clause> _clauses;
    Array* _unsat_clauses;                  // the set of unsat clauses
    Array* _sat_clause_with_false_literal;  // clauses with 0<sat_num<literal_num, from which swap operation are choosen
    Array* _contain_bool_unsat_clauses;     // unsat clause with at least one boolean var
    // diff logic
    Array* pair_x;  // x-y-->z
    Array* pair_y;
    std::vector<__int128_t> pair_x_value;  // x - y
    std::vector<__int128_t> pair_y_value;
    std::vector<std::pair<int, int>> equal_table;

    /* auxiliary variables */
    const uint64_t _additional_len;

    /* clause weighting */
    uint64_t _total_clause_weight;  // total weight of all clauses

    /* control */
    // for overflow
    bool is_overflow = false;
    bool update_sampling_interval = false;
    int _last_flip_lia_lit = -1;  // 记录上一步操作是反转的哪个文字
    // for random

    std::mt19937 mt;  // random number generator
    // step
    uint64_t _step;              // the number of steps executed by the algorithm
    uint64_t _outer_layer_step;  // the steps of the outer loop
    const uint64_t _max_step;    // the maximum allowed number of steps
    // data structure for clause weighting
    const uint64_t smooth_probability;  // the smoothing probability used in the clause weighting algorithm
    uint64_t _swt_threshold;            // the threshold for switching or modifying the clause weight
    float _swt_p;                       // w=w*p+ave_w*(1-p)
    uint64_t total_clause_weight;       // total weight of all clauses
    int _lit_in_unsat_clause_num;       // the number of literals in unsatisfied clauses
    int _bool_lit_in_unsat_clause_num;  // the number of Boolean literals in unsatisfied clauses
    // boolean variable
    bool use_swap_from_from_small_weight;  // ??
    bool use_pbs = false;
    bool is_pb = true;
    bool is_idl = true;  // if it is the IDL mode
    // map
    std::map<std::string, uint64_t> _name2resolution_var;  // map the name of a resolution variable to its index
    std::map<std::string, uint64_t> _name2var;             // map the name of a variable to its index
    std::map<std::string, uint64_t> _name2tmp_var;  // map from temporary variable names to their indices
    // cc and tabu
    int CC_mode;                      // controlling the variable selection mode (integer mode or boolean mode)
    std::vector<int> _CClist;         // CClist[2*var]==1 and [2*var+1]==1 means that var is allowed to move forward or backward
    std::vector<uint64_t> _tabulist;  // tabulist[2*var] and [2*var+1] means that var are forbidden to move forward or backward until then
                                      // Move forward: flips the variable from false (0) to true (1)
                                      // Move backward: flips the variable from true (1) to false (0)
    // vector
    std::vector<int> _operation_var_idx_vec;              // indices of variables involved in operations
    std::vector<__int128_t> _operation_change_value_vec;  // the change values for the variables involved in operations
    std::vector<int> _operation_lit_idx_vec;
    std::vector<int> _operation_var_idx_bool_vec;  // the indices of Boolean variables involved in operations
    std::vector<uint64_t> _last_move;              // the last move step for each variable. (上次移动是在第几步)
                                                   // for lia operations, _last_move[2*var+1] means change_value > 0; _last_move[2*var+0] means change_value <= 0;
    std::vector<bool> _is_chosen_bool_var;         // indicating which Boolean variables are chosen for operations
    std::vector<int> _pre_value_1;                 // the 1st pre-set value of a var, if the var is in the form of (a==0 OR a==1)
    std::vector<int> _pre_value_2;                 // the 2nd pre-set value of a var
    // time
    // double _total_time;
    double _best_cost_time;  // the time when the best solution was found
    std::chrono::steady_clock::time_point _total_start;
    double _cutoff;  // the cutoff time ??

    /* solution */
    // solution
    // int _solution_size;
    std::vector<__int128_t> _solution;      // current solution
    std::vector<__int128_t> _best_solutin;  // optimal solution
    std::vector<__int128_t> _final_solution;
    std::stack<clause> _reconstruct_stack;      // 被归结的子句
    bool is_in_bool_search = false;             // ??

    // cost
    int best_found_this_restart;          // cost of the best solution found in the current restart process
    int _best_found_hard_cost_this_bool;  // ??
    int _best_found_hard_cost_this_lia;   // ??
    int _no_improve_cnt_bool = 0;         // the number of times an improved solution was not found during the search for Boolean variables
    int _no_improve_cnt_lia = 0;

   public:
    int _best_found_cost;  // cost of finding the optimal solution
    /* parse input */
    void split_string(std::string& in_string, std::vector<std::string>& str_vec, std::string pattern);
    void build_lits(std::string& in_string);
    void build_instance(std::vector<std::vector<int>>& clause_vec);
    uint64_t transfer_name_to_resolution_var(std::string& name, bool is_lia, bool in_equal);
    uint64_t transfer_name_to_tmp_var(std::string& name, bool in_equal);                   // lia var is first inserted into _tmp_var when build lit,
                                                                                           // then inserted into _resolution_var when reduce var(x-y->z)
    uint64_t transfer_name_to_reduced_var(std::string& name, bool is_lia, bool in_equal);  // after the resolution vars are inserted into _vars

    /* initialize */
    ls_sampler(unsigned seed, std::string input, std::string output)
        : _random_seed(seed),
          _input_path(input),
          _output_path(output),
          _additional_len(10),
          CC_mode(-1),
          _max_step(UINT64_MAX),
          smooth_probability(3),
          _swt_p(0.3),
          _swt_threshold(50),
          _cutoff(1200) {
        std::random_device rd;
        mt.seed(rd());
    }
    void make_space();
    void free_space();
    void make_lits_space(uint64_t num_lits) {
        _num_lits = num_lits;
        _lits.resize(num_lits + _additional_len);
    };
    void initialize();
    void initialize_variable_datas();
    void initialize_lit_datas();
    void initialize_clause_datas();
    void build_neighbor();
    void unit_prop();
    void resolution();
    void reduce_clause();
    void set_pre_value();
    void reduce_vars();

    // random walk
    void update_clause_weight();
    void smooth_clause_weight();
    void random_walk();

    /* basic operation */
    // literals
    void invert_lit(lit& l);
    __int128_t delta_lit(lit& l);
    // move
    void modify_CC(uint64_t var_idx, int direction);
    int pick_critical_move(__int128_t& best_value);
    int flip_lia_lit(lit* l, __int128_t& best_value, bool to_ture);
    int pick_critical_move_bool();
    void critical_move(uint64_t var_idx, __int128_t change_value);
    void insert_operation(int var_idx, __int128_t change_value, int& operation_idx);
    void insert_operation(int var_idx, __int128_t change_value, int& operation_idx, int const& flip_lit);
    __int128_t devide(__int128_t a, __int128_t b);
    void add_swap_operation(int& operation_idx);
    void interval_move(uint64_t var_idx, __int128_t& change_value);

    // clauses
    inline void unsat_a_clause(uint64_t clause_idx) {
        _unsat_clauses->insert_element((int)clause_idx);
        if (_clauses[clause_idx].bool_literals.size() > 0)
            _contain_bool_unsat_clauses->insert_element((int)clause_idx);
    };
    inline void sat_a_clause(uint64_t clause_idx) {
        _unsat_clauses->delete_element((int)clause_idx);
        _contain_bool_unsat_clauses->delete_element((int)clause_idx);
    };
    // distance/score
    void convert_to_pos_delta(__int128_t& delta, int l_idx);  // Convert the value of delta to a positive value
    __int128_t critical_score(uint64_t var_idx, __int128_t change_value);
    __int128_t critical_subscore(uint64_t var_idx, __int128_t change_value);  // calculate dscore
    void critical_score_subscore(uint64_t var_idx, __int128_t change_value);
    void critical_score_subscore(uint64_t var_idx);  // dedicated for boolean var
    // clean
    void clear_prev_data();
    void construct_solution_score();  // construct the solution based on score
    bool update_best_solution();

    // time
    double TimeElapsed_total();

    /* ls_sampling */
    void ls_sampling();
    void shrinkSampleInterval(lit* l);

    /* search */
    bool search();
    void swap_from_small_weight_clause();
    void enter_lia_mode();
    void enter_bool_mode();
    bool update_outer_best_solution();
    bool update_inner_best_solution();
    void up_bool_vars();
    void choose_value_for_pair();
    void build_equal_table();
    void build_occs_closures();
    /* random */
    __int128_t random_int128_in_range(__int128_t low, __int128_t up);

    int64_t random_int64_in_range(int64_t left, int64_t right);

    // check
    bool check_solution();

    /* prety print*/
    void print_var_solution(std::string& var_name, std::string& var_value);
    std::string print_128(__int128 n);
    void print_vars(std::ostream& out);
    inline __int128_t abs_128(__int128_t n) { return n >= 0 ? n : -n; }
    void print_literal(std::ostream& out, lit& l);
    void print_formula(std::ostream& out);
    void print_interal_data_strcture(std::ostream& out);
    void print_move_status(std::ostream& out, int v_idx, __int128_t change_val);
    void print_result();
};
}  // namespace sampler