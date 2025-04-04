#include "CosineLaser.hpp"
#include <cmath>

namespace FHHG::Laser {
    CosineLaser::CosineLaser(h_float photon_energy, h_float E_0, int max_frequency, int n_subdivision, int n_laser_cycles)
        : Laser(photon_energy, E_0, max_frequency, n_subdivision), n_laser_cycles(n_laser_cycles) 
    {}

    void CosineLaser::laser_convolution(const State &current_state, State &new_state, const h_float m_x, const h_float m_z) const
    {
        auto prefactor = [this](h_float omega) -> h_complex {
            return (omega / sqrt_two_pi) * imaginary_unit * (std::polar(1.0, 2.0 * pi * n_laser_cycles * omega) - 1.0);
        };
        auto main_kernel = [this, &prefactor](h_float omega) -> h_complex {
            if (is_zero(std::abs(omega) - 1.0)) {
                return (-pi / sqrt_two_pi) * n_subdivisions;
            }
            else if (is_zero(std::abs(omega) - (1. + 1. / n_laser_cycles))) {
                return (-0.5 * pi / sqrt_two_pi) * n_subdivisions / (1. + 1. / n_laser_cycles);
            }
            else if (is_zero(std::abs(omega) - (1. - 1. / n_laser_cycles))) {
                return (-0.5 * pi / sqrt_two_pi) * n_subdivisions / (1. - 1. / n_laser_cycles);
            }

            const h_float omega_squared = omega * omega;
            return prefactor(omega) * (
                1.   / (1. - omega_squared)
                -0.5 / ((1. + 1. / n_laser_cycles) * (1. + 1. / n_laser_cycles) - omega_squared)
                -0.5 / ((1. - 1. / n_laser_cycles) * (1. - 1. / n_laser_cycles) - omega_squared)
            );
        };

        h_float omega, integration_variable;
        h_complex kernel_minus, kernel_plus;
        for (int i = 0; i < N; ++i) {
            omega = frequency(i);
            for (int j = 0; j < N; ++j) {
                integration_variable = frequency(j);

                kernel_minus = main_kernel(omega - integration_variable);
                kernel_plus = main_kernel(omega + integration_variable);

                new_state(i, 0) += (m_z / integration_variable) * (kernel_minus * current_state(j, 1) + kernel_plus * std::conj(current_state(j, 1)));
                new_state(i, 1) -= (m_z / integration_variable) * (kernel_minus * current_state(j, 0) + kernel_plus * std::conj(current_state(j, 0)));
                new_state(i, 1) += (m_x / integration_variable) * (kernel_minus * current_state(j, 2) + kernel_plus * std::conj(current_state(j, 2)));
                new_state(i, 1) -= (m_x / integration_variable) * (kernel_minus * current_state(j, 1) + kernel_plus * std::conj(current_state(j, 1)));
            }
        }
        new_state /= N; // normalization factor
    }
}