#include "DiracSystem.hpp"

namespace FHHG {
    DiracSystem::DiracSystem(h_float temperature, h_float _E_F, h_float _v_F, h_float _band_width, h_float _photon_energy)
        : beta { is_zero(temperature) ? std::numeric_limits<h_float>::infinity() : 1. / (k_B * temperature * _photon_energy) },
        E_F{ _E_F / _photon_energy }, 
        v_F{ _v_F * ((1e12 * hbar) / _photon_energy) }, // 1e12 for conversion to pm; T_L = hbar / _photon_energy
        band_width{ _band_width },
        max_k { band_width }, // in units of omega_L / v_F
        max_kappa_compare { band_width * band_width }  // in units of (omega_L / v_F)^2
    {  }

    void DiracSystem::compute_single_current_density(nd_vector& rhos, Laser::Laser const * const laser, h_float k_z, h_float kappa) const
    {
        State current_state(laser->n_subdivisions * laser->max_frequency);
        State new_state(laser->n_subdivisions * laser->max_frequency);

        for (int i = 0; i < current_state.n_frequencies; ++i) {
            current_state(i) = { h_complex{}, h_complex{}, h_complex{ 2. / static_cast<h_float>(i + 1), 0.} };
        }

        const h_float magnitude_k = norm(k_z, kappa);
        const h_float m_x = laser->momentum_amplitude * 2.0 * v_F * kappa / magnitude_k;
        const h_float m_z = -2.0 * laser->momentum_amplitude * magnitude_k;

        const h_complex diagonal_factor = 2.0 * imaginary_unit * magnitude_k / sqrt(2. * pi);

        bool go_on = true;
        while(go_on) {
            for (size_t i = 0U; i < current_state.n_frequencies; ++i) {
                new_state(i, 0) = (diagonal_factor / laser->frequency(i)) * current_state(i, 1);
                new_state(i, 1) = -(diagonal_factor / laser->frequency(i)) * current_state(i, 0);
                new_state(i, 2) = h_complex{};
            }
            laser->laser_convolution(current_state, new_state, m_x, m_z);

            for (size_t i = 0U; i < current_state.n_frequencies; ++i) {
                State::state_type diff = new_state(i) - current_state(i);
                if (diff.norm() / new_state(i).norm() > 1e-4) {
                    break;
                }
                if (i == current_state.n_frequencies - 1) {
                    go_on = false;
                }
            }
            current_state = new_state;
        }
    }

    std::string DiracSystem::info() const
    {
        return "DiracSystem\nE_F=" + std::to_string(E_F) + " * hbar omega_L"
            + "\nv_F=" + std::to_string(v_F) + " * pm / T_L"
            + "\nband_width=" + std::to_string(band_width)
            + "\nmax_k=" + std::to_string(max_k) + " pm"
            + "\nmax_kappa_compare=" + std::to_string(max_kappa_compare) + " pm^2";
    }
}