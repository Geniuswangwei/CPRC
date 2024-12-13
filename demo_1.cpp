#include "solver.hpp"
#include <iostream>
#include <iomanip>

void print_complex(const std::string& name, const std::complex<double>& value) {
    std::cout << std::setw(15) << name << ": "
              << std::setprecision(4) << std::fixed 
              << value.real() << " + j" << value.imag() << std::endl;
}

int main() {
    // 定义基本参数
    double epsilon_b = 4.0;
    std::complex<double> epsilon_a(12.0, -1.0);
    double mu_b = 1.0;
    std::complex<double> mu_a(3.0, -1.0);

    // 创建求解器实例
    HoneycombSolver solver(epsilon_b, epsilon_a, mu_b, mu_a);

    // 结构参数
    double t = 0.5;  // mm
    double x = 5.0;  // mm
    double d = (x - t) * 0.05;  // 使用d/(x-t) = 0.05计算d

    // 计算等效参数
    auto epsilon_r = solver.calculate_epsilon_r(d, x, t, epsilon_a);
    auto epsilon_p = solver.calculate_epsilon_p(d, x, t);

    auto epsilon_parallel = solver.calculate_epsilon_parallel(x, t, epsilon_b, epsilon_p);
    auto epsilon_perpendicular = solver.calculate_epsilon_perpendicular(t, x, epsilon_r, epsilon_b);

    auto mu_r = solver.calculate_epsilon_r(d, x, t, mu_a);
    auto mu_p = solver.calculate_mu_p(d, x, t);

    auto mu_parallel = solver.calculate_mu_parallel(x, t, mu_b, mu_p);
    auto mu_perpendicular = solver.calculate_mu_perpendicular(t, x, mu_r, mu_b);

    // 打印结果
    std::cout << "计算结果 (d/(x-t) = 0.05):" << std::endl;
    std::cout << "----------------------------------------" << std::endl;
    print_complex("epsilon //", epsilon_parallel);
    print_complex("epsilon ⊥", epsilon_perpendicular);
    print_complex("mu //", mu_parallel);
    print_complex("mu ⊥", mu_perpendicular);
    
    // 参考值
    std::cout << "\n参考值 (d/(x-t) = 0.05):" << std::endl;
    std::cout << "----------------------------------------" << std::endl;
    std::cout << "epsilon //: 2.4387 - j0.0790" << std::endl;
    std::cout << "epsilon ⊥:  1.4827 - j0.0036" << std::endl;
    std::cout << "mu //:      1.1580 - j0.0790" << std::endl;
    std::cout << "mu ⊥:       1.0890 - j0.0219" << std::endl;

    return 0;
}