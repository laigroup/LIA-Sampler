#include "liasampler.h"
#include <filesystem>  // C++17 引入的库
#include <fstream>
#include <regex>

namespace sampler {

std::string processNegNumber(const std::string& input) {
    // 定义正则表达式以匹配类似 "(- 4294967281)" 的模式
    std::regex pattern(R"(\(\s*(-\s*\d+)\s*\))");
    std::string result;

    // 使用 std::sregex_iterator 进行匹配和替换
    std::sregex_iterator begin(input.begin(), input.end(), pattern);
    std::sregex_iterator end;

    size_t lastPos = 0;
    for (auto it = begin; it != end; ++it) {
        // 追加匹配前的部分
        result += input.substr(lastPos, it->position() - lastPos);
        // 去掉捕获组中的空格
        std::string number = (*it)[1].str();
        number.erase(std::remove_if(number.begin(), number.end(), ::isspace), number.end());
        // 追加处理后的部分
        result += number;
        // 更新位置
        lastPos = it->position() + it->length();
    }

    // 追加剩余部分
    result += input.substr(lastPos);

    return result;
}

std::string extract_filename(const std::string& path) {
    // 找到最后一个 '/' 的位置
    size_t pos = path.find_last_of('/');
    if (pos != std::string::npos) {
        return path.substr(pos + 1);
    }
    return path;  // 如果没有 '/'，假设整个路径就是文件名
}

void LiaSampler::print_statistic() {
    std::cout << "--------------------- After sampling: statistic ---------------------\n";
    std::cout << "Sampling time: " << TimeElapsed() << "\n";
    std::cout << "Total samples number: " << num_samples << "\n";
}

double LiaSampler::TimeElapsed() {
    std::chrono::steady_clock::time_point finish = std::chrono::steady_clock::now();
    std::chrono::duration<double> duration = finish - time_sampling_start;
    return duration.count();
}

z3::tactic LiaSampler::mk_preamble_tactic(z3::context& ctx) {
    z3::params pull_ite_p(ctx);

    pull_ite_p.set("pull_cheap_ite", true);
    pull_ite_p.set("push_ite_arith", false);
    pull_ite_p.set("local_ctx", true);
    pull_ite_p.set("local_ctx_limit", 10000000u);
    pull_ite_p.set("hoist_ite", true);

    z3::params ctx_simp_p(ctx);
    ctx_simp_p.set("max_depth", 30u);
    ctx_simp_p.set("max_steps", 5000000u);

    z3::params solver_eqs_p(ctx);
    solver_eqs_p.set("solve_eqs_max_occs", 5U);

    z3::params lhs_p(ctx);
    lhs_p.set("arith_lhs", true);

    return z3::tactic(ctx, "simplify") &
           z3::tactic(ctx, "propagate-values") &
           z3::with(z3::tactic(ctx, "ctx-simplify"), ctx_simp_p) &
           z3::with(z3::tactic(ctx, "simplify"), pull_ite_p) &
           z3::with(z3::tactic(ctx, "solve-eqs"), solver_eqs_p) &
           z3::tactic(ctx, "elim-uncnstr") &
           z3::with(z3::tactic(ctx, "simplify"), lhs_p);
}

void LiaSampler::sampling() {
    time_sampling_start = std::chrono::steady_clock::now();

    z3::solver ls_solver(c);
    ls_solver.set("sampling", true);
    ls_solver.set("input_file", smtFilePath.c_str());
    ls_solver.set("output_file", samplesFileDir.c_str());
    ls_solver.set("seed", seed);
    ls_solver.set("random_seed", seed);

    parseSmtFile();

    z3::goal g(c);
    g.add(original_formula);
    z3::tactic preamble_tactic = mk_preamble_tactic(c);
    z3::apply_result simp_ar = preamble_tactic(g);

    assert(simp_ar.size() == 1);
    z3::goal subgoal = simp_ar[0];

    for (unsigned i = 0; i < subgoal.size(); i++) {
        ls_solver.add(subgoal[i]);
    }

    // ls_solver.add(simp_formula);

    std::string samplesFileName = samplesFileDir + "/" + extract_filename(smtFilePath) + ".samples";
    std::ofstream samplesFile(samplesFileName);  // 打开文件
    if (!samplesFile) {
        std::cerr << "Unable to open file " << samplesFileName << std::endl;
        return;
    }

    for (size_t i = 0; i < maxNumSamples; ++i) {
        z3::check_result check_res = ls_solver.check();
        if (z3::sat != check_res) {
            std::cout << "Unsat or known case!\n";
            return;
        }
        z3::model m = ls_solver.get_model();
        m = subgoal.convert_model(m);
        samplesFile << i << ": ";
        for (size_t j = 0; j < m.size(); ++j) {
            if (m[j].is_const()) {
                samplesFile << m[j].name() << ":" << processNegNumber(m.get_const_interp(m[j]).to_string()) << ";";
            }
        }
        samplesFile << "\n";
        num_samples ++;
#if defined(VERBOSE) && defined(PRINT_PROGRESS)
        std::cout << "The " << i + 1 << " sample is being generated ..." << std::endl;
#endif
        if (TimeElapsed() > maxTimeLimit){
            break;
        }
    }

    samplesFile.close();

    print_statistic();
}
}  // namespace sampler