#pragma once
#include <iostream>
#include <vector>
#include <complex>
#include <cmath>

class HoneycombSolver {
public:
    // Physical constants
    static constexpr double EPSILON_0 = 8.854e-12;  // vacuum permittivity (F/m)
    static constexpr double MU_0 = 4 * M_PI * 1e-7; // vacuum permeability (H/m)
    static constexpr double C = 3e8;                // speed of light (m/s)

    // Constructor with electromagnetic parameters
    HoneycombSolver(double epsilon_b, std::complex<double> epsilon_a,
                   double mu_b, std::complex<double> mu_a)
        : epsilon_b_(epsilon_b), epsilon_a_(epsilon_a),
          mu_b_(mu_b), mu_a_(mu_a) {}

    // Calculate epsilon_p
    std::complex<double> calculate_epsilon_p(double d, double x, double t) const {
        double ratio = d / (x - t);
        auto term1 = ratio * (2 - ratio) * epsilon_a_;
        auto term2 = std::pow(1 - ratio, 2);
        return term1 + term2;
    }

    // Calculate mu_p
    std::complex<double> calculate_mu_p(double d, double x, double t) const {
        double ratio = d / (x - t);
        auto term1 = ratio * (2 - ratio) * mu_a_;
        auto term2 = std::pow(1 - ratio, 2);
        return term1 + term2;
    }

    // Calculate epsilon_r or mu_r
    std::complex<double> calculate_epsilon_r(double d, double x, double t, 
                                          std::complex<double> epsilon_a) const {
        double ratio = d / (x - t);
        auto term1 = std::pow(1 - ratio, 2) - ratio * (2 - ratio);
        auto term2 = term1 * (1.0 - epsilon_a);
        auto term3 = std::sqrt(std::pow(term1, 2) * std::pow(1.0 - epsilon_a, 2) + 4.0 * epsilon_a);
        return 0.5 * (term2 + term3);
    }

    // Calculate parallel component (epsilon or mu)
    std::complex<double> calculate_epsilon_parallel(double x, double t,
                                          std::complex<double> epsilon_b,
                                          std::complex<double> epsilon_p) const {
        double ratio = t / x;
        auto term1 = ratio * (2 - ratio) * epsilon_b;
        auto term2 = std::pow(1 - ratio, 2) * epsilon_p;
        return term1 + term2;
    }

    std::complex<double> calculate_mu_parallel(double x, double t,
                                          std::complex<double> mu_b,
                                          std::complex<double> mu_p) const {
        double ratio = t / x;
        auto term1 = ratio * (2 - ratio) * mu_b;
        auto term2 = std::pow(1 - ratio, 2) * mu_p;
        return term1 + term2;
    }

    // Calculate perpendicular component (epsilon or mu)
    std::complex<double> calculate_epsilon_perpendicular(double t, double x,
                                               std::complex<double> epsilon_r,
                                               std::complex<double> epsilon_b) const {
        double ratio = t / x;
        auto term1 = std::pow(1 - ratio, 2) - ratio * (2 - ratio);
        auto term2 = term1 * (epsilon_r - epsilon_b);
        auto term3 = std::sqrt(std::pow(term1, 2) * std::pow(epsilon_r - epsilon_b, 2) + 
                              4.0 * epsilon_r * epsilon_b);
        return 0.5 * (term2 + term3);
    }

    std::complex<double> calculate_mu_perpendicular(double t, double x,
                                               std::complex<double> mu_r,
                                               std::complex<double> mu_b) const {
        double ratio = t / x;
        auto term1 = std::pow(1 - ratio, 2) - ratio * (2 - ratio);
        auto term2 = term1 * (mu_r - mu_b);
        auto term3 = std::sqrt(std::pow(term1, 2) * std::pow(mu_r - mu_b, 2) + 
                              4.0 * mu_r * mu_b);
        return 0.5 * (term2 + term3);
    }

    // Calculate wave number k0
    double calculate_k0(double f) const {
        double omega = 2 * M_PI * f;
        return omega * std::sqrt(EPSILON_0 * MU_0);
    }

    // Calculate reflection coefficient
    double calculate_reflection_coefficient(std::complex<double> epsilon_perp,
                                         std::complex<double> mu_perp,
                                         double f, double h) const {
        double k0 = calculate_k0(f);
        
        auto sqrt_epsilon = std::sqrt(epsilon_perp);
        auto sqrt_mu = std::sqrt(mu_perp);
        
        auto gamma = k0 * std::sqrt(mu_perp * epsilon_perp);
        auto exp_term = std::exp(std::complex<double>(0, 2) * gamma * h);
        
        auto numerator = (sqrt_epsilon + sqrt_mu) + (sqrt_epsilon - sqrt_mu) * exp_term;
        auto denominator = (sqrt_epsilon - sqrt_mu) + (sqrt_epsilon + sqrt_mu) * exp_term;
        
        return 20 * std::log10(std::abs(numerator / denominator));
    }

    // Getters for electromagnetic parameters
    double get_epsilon_b() const { return epsilon_b_; }
    std::complex<double> get_epsilon_a() const { return epsilon_a_; }
    double get_mu_b() const { return mu_b_; }
    std::complex<double> get_mu_a() const { return mu_a_; }

private:
    double epsilon_b_;              // relative permittivity of frame
    std::complex<double> epsilon_a_; // relative permittivity of absorbing layer
    double mu_b_;                   // relative permeability of frame
    std::complex<double> mu_a_;     // relative permeability of absorbing layer
};