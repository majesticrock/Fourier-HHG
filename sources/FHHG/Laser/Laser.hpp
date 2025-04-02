#pragma once
#include "../GlobalDefinitions.hpp"
#include "../State.hpp"

namespace FHHG::Laser {
    struct Laser {
        const h_float momentum_amplitude{}; // e E_0 / (hbar omega_L)
        const int max_frequency{}; // in units of hbar omega_L
        const int n_subdivisions{}; // number of subdivisions per hbar omega_L

        /**
         * @param photon_energy \f$ \hbar \omega_L \f$ in meV
         * @param E_0 peak electric field strength in MV / cm
         */
        Laser(h_float photon_energy, h_float E_0, int max_frequency, int n_subdivisions);

        inline h_float frequency(int i) const noexcept {
            return (static_cast<h_float>(i) + 0.5) / static_cast<h_float>(n_subdivisions);
        }

        virtual void laser_convolution(const State& current_state, State& new_state, const h_float m_x, const h_float m_z) const = 0;
    };
}