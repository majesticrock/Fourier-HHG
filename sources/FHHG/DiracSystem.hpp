#pragma once
#include "GlobalDefinitions.hpp"
#include "Laser/Laser.hpp"
#include <string>

namespace FHHG {
    class DiracSystem {
    public:
        DiracSystem() = delete;
        /**
         * @param _E_F Fermi energy in meV
         * @param _v_F Fermi velocity in m/s
         * @param _band_width in multiples of the photon energy
         * @param _photon_energy hbar omega_L in meV
         */
        DiracSystem(h_float temperature, h_float _E_F, h_float _v_F, h_float _band_width, h_float _photon_energy);

        h_float dispersion(h_float k_z, h_float kappa) const;

        void compute_single_current_density(ncd_vector& j_buffer, Laser::Laser const * const laser, h_float k_z, h_float kappa) const;

        std::vector<h_complex> compute_current_density(Laser::Laser const * const laser, const int n_z, const int n_kappa = 20,
            const h_float kappa_threshold = 1e-3, std::string const& debug_dir = "") const;

        std::string info() const;

    private:
        inline h_float z_integration_upper_limit() const noexcept { return max_k; }
        h_float kappa_integration_upper_limit(h_float k_z) const;

        const h_float beta{}; ///< in units of the 1 / photon energy
        const h_float E_F{}; ///< in units of the photon energy
        const h_float v_F{}; ///< in units of pm / T_L, where T_L = 1 / omega_L
        const h_float band_width{}; ///< in units of the photon energy
        const h_float max_k{}; ///< in units of omega_L / v_F
        const h_float max_kappa_compare{}; ///< in units of (omega_L / v_F)^2
    };
}