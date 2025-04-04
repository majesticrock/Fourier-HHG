#include "ContinuousLaser.hpp"
#include <cassert>

namespace FHHG::Laser {
    ContinuousLaser::ContinuousLaser(h_float photon_energy, h_float E_0, int max_frequency, int n_subdivision)
        : Laser(photon_energy, E_0, max_frequency, n_subdivision) {}

    void ContinuousLaser::laser_convolution(const State& current_state, State& new_state, const h_float m_x, const h_float m_z) const {
        constexpr h_float factor = sqrt_pi_over_two;
        for(size_t i = n_subdivisions; i < current_state.n_frequencies; ++i) {
            new_state(i, 0) += (factor / frequency(i - n_subdivisions)) * m_z * current_state(i - n_subdivisions, 1);
            new_state(i, 1) += (factor / frequency(i - n_subdivisions)) * (m_x * current_state(i - n_subdivisions, 2) - m_z * current_state(i - n_subdivisions, 0));
            new_state(i, 2) -= (factor / frequency(i - n_subdivisions)) * m_x * current_state(i - n_subdivisions, 1);
        }

        for(size_t i = 0U; i < current_state.n_frequencies - n_subdivisions; ++i) {
            new_state(i, 0) += (factor / frequency(i + n_subdivisions)) * m_z * current_state(i + n_subdivisions, 1);
            new_state(i, 1) += (factor / frequency(i + n_subdivisions)) * (m_x * current_state(i + n_subdivisions, 2) - m_z * current_state(i + n_subdivisions, 0));
            new_state(i, 2) -= (factor / frequency(i + n_subdivisions)) * m_x * current_state(i + n_subdivisions, 1);
        }

        for(size_t i = 0U; i < n_subdivisions; ++i) {
            new_state(i, 0) += (factor / frequency(n_subdivisions - 1 - i)) * m_z * current_state(n_subdivisions - 1 - i, 1);
            new_state(i, 1) += (factor / frequency(n_subdivisions - 1 - i)) * (m_x * current_state(n_subdivisions - 1 - i, 2) - m_z * current_state(n_subdivisions - 1 - i, 0));
            new_state(i, 2) -= (factor / frequency(n_subdivisions - 1 - i)) * m_x * current_state(n_subdivisions - 1 - i, 1);
        }
    }
}