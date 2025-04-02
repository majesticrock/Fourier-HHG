#pragma once

#include "Laser.hpp"

namespace FHHG::Laser {
    struct ContinuousLaser : public Laser {
        ContinuousLaser(h_float photon_energy, h_float E_0, int max_frequency, int n_subdivision);

        void laser_convolution(const State& current_state, State& new_state, const h_float m_x, const h_float m_z) const final;
    };
}