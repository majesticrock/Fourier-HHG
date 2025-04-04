#pragma once
#include "Laser.hpp"

namespace FHHG::Laser {
    struct CosineLaser : public Laser {
        int n_laser_cycles{};

        CosineLaser(h_float photon_energy, h_float E_0, int max_frequency, int n_subdivision, int n_laser_cycles);

        void laser_convolution(const State& current_state, State& new_state, const h_float m_x, const h_float m_z) const final;
    };
}