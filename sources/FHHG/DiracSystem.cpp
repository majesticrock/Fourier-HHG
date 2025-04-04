#include "DiracSystem.hpp"

#include <mrock/utility/Numerics/Integration/AdaptiveTrapezoidalRule.hpp>
#include <mrock/utility/Numerics/Integration/AdaptiveTrapezoidalRule.hpp>
#include <mrock/utility/Numerics/ErrorFunctors.hpp>
#include <mrock/utility/progress_bar.hpp>

#include <iostream>
#include <omp.h>

#include <nlohmann/json.hpp>
#include <mrock/utility/OutputConvenience.hpp>

#pragma omp declare reduction(vec_plus : std::vector<FHHG::h_complex> : \
    std::transform(omp_out.begin(), omp_out.end(), omp_in.begin(), omp_out.begin(), std::plus<FHHG::h_complex>())) \
    initializer(omp_priv = decltype(omp_orig)(omp_orig.size()))

namespace FHHG {
    DiracSystem::DiracSystem(h_float temperature, h_float _E_F, h_float _v_F, h_float _band_width, h_float _photon_energy)
        : beta { is_zero(temperature) ? std::numeric_limits<h_float>::infinity() : 1. / (k_B * temperature * _photon_energy) },
        E_F{ _E_F / _photon_energy }, 
        v_F{ _v_F * ((1e12 * hbar) / _photon_energy) }, // 1e12 for conversion to pm; T_L = hbar / _photon_energy
        band_width{ _band_width },
        max_k { band_width }, // in units of omega_L / v_F
        max_kappa_compare { band_width * band_width }  // in units of (omega_L / v_F)^2
    {  }

    h_float DiracSystem::dispersion(h_float k_z, h_float kappa) const
    {
        return norm(kappa, k_z); // v_F is already contained within the k values
    }

    void DiracSystem::compute_single_current_density(ncd_vector& j_buffer, Laser::Laser const * const laser, h_float k_z, h_float kappa) const
    {
        State current_state(laser->N);
        State new_state(current_state.n_frequencies);

        for (int i = 0; i < current_state.n_frequencies; ++i) {
            current_state(i) = { h_complex{1e-4, 1e-4}, h_complex{1e-4, -1e-4}, h_complex{ 2. / static_cast<h_float>(i + 1), 0.1} };
        }

        const h_float magnitude_k = norm(k_z, kappa);
        const h_float m_x =  2.0 * laser->momentum_amplitude * v_F * kappa / magnitude_k;
        const h_float m_z = -2.0 * laser->momentum_amplitude * v_F * k_z / magnitude_k;

        const h_float diagonal_factor = magnitude_k * 2.0 * sqrt_two_pi;

        while (!current_state.is_close_relative(new_state, 1e-8)) {
            new_state.setZero();
            laser->laser_convolution(current_state, new_state, m_x, m_z);
            for (size_t i = 0U; i < current_state.n_frequencies; ++i) {
                new_state(i, 0) += (diagonal_factor / laser->frequency(i)) * current_state(i, 1);
                new_state(i, 1) -= (diagonal_factor / laser->frequency(i)) * current_state(i, 0);
            }
            current_state = -imaginary_unit * new_state;
        }

        const h_float alpha_0 = fermi_function(E_F + dispersion(k_z, kappa), beta);
        const h_float beta_0 = fermi_function(E_F - dispersion(k_z, kappa), beta);
        const h_float rho_x = alpha_0 * beta_0;
        const h_float rho_z = alpha_0 * alpha_0 - beta_0 * beta_0;

        for (int i = 0; i < current_state.n_frequencies; ++i) {
            //std::cout << current_state(i, 0) << " " << current_state(i, 1) << " " << current_state(i, 2) << std::endl;
            j_buffer[i] = current_state(i, 0) * rho_x + current_state(i, 2) * rho_z; // factor of 2 does not matter
        }
    }

    std::vector<h_complex> DiracSystem::compute_current_density(Laser::Laser const * const laser, const int n_z, const int n_kappa/* = 20 */, 
        const h_float kappa_threshold /* = 1e-3 */, std::string const& debug_dir/* = "" */) const
    {
        constexpr mrock::utility::Numerics::Integration::adapative_trapezoidal_rule_print_policy m_policy{false, false, false};
        mrock::utility::Numerics::Integration::adapative_trapezoidal_rule<h_float, m_policy> integrator;

        const int N = laser->n_subdivisions * laser->max_frequency;
        ncd_vector j_buffer = nd_vector::Zero(N);
        std::vector<h_complex> current_density(N, h_complex{});
        const auto delta_z = 2.0 * z_integration_upper_limit() / n_z;

        auto integration_weight = [](h_float k_z, h_float kappa) {
            return k_z * kappa / norm(k_z, kappa);
        };

        std::vector<int> progresses(omp_get_max_threads(), int{});

        /* Debug output config */
        std::vector<h_float> k_zs(n_z);
        std::vector<h_float> j_k_0(n_z);
        std::vector<h_float> j_k_1(n_z);
        std::vector<h_float> j_k_2(n_z);
        std::vector<h_float> j_k_3(n_z);
        k_zs.front() = - z_integration_upper_limit();
        k_zs.back() = z_integration_upper_limit();
        /* End debug config */

#pragma omp parallel for firstprivate(j_buffer) reduction(vec_plus:current_density) schedule(dynamic)
        for (int z = 1; z < n_z; ++z) { // f(|z| = z_max) = 0
            ++(progresses[omp_get_thread_num()]);
            if (omp_get_thread_num() == 0) {
                mrock::utility::progress_bar( 
                    static_cast<float>(std::reduce(progresses.begin(), progresses.end())) / static_cast<float>(n_z)
                );
            }
            const auto k_z = (z - n_z / 2) * delta_z;

            auto kappa_integrand = [&](h_float kappa) -> const ncd_vector& {
                if (is_zero(kappa)) {
                    j_buffer.setZero();
                }
                else {
                    compute_single_current_density(j_buffer, laser, k_z, kappa);
                    j_buffer *= integration_weight(k_z, kappa);
                }
                return j_buffer;
            };

            auto kappa_result = integrator.split_integrate<10>(kappa_integrand, h_float{}, kappa_integration_upper_limit(k_z), 
                n_kappa, kappa_threshold, mrock::utility::Numerics::vector_elementwise_error<ncd_vector, ncd_vector::value_type, false>(), ncd_vector::Zero(N));

            k_zs[z] = k_z;
            j_k_0[z] = std::abs(kappa_result(0));
            j_k_1[z] = std::abs(kappa_result(N / 5));
            j_k_2[z] = std::abs(kappa_result(N / 4));
            j_k_3[z] = std::abs(kappa_result(N / 2));
            
            std::transform(current_density.begin(), current_density.end(), kappa_result.begin(), current_density.begin(), std::plus<>());
        }
        std::cout << std::endl;
        for (int i = 0; i < N; ++i) {
            current_density[i] *= delta_z;
        }

        if (debug_dir != "") {
            nlohmann::json debug_json = {
                { "time",   mrock::utility::time_stamp() },
                { "k_zs",   k_zs },
                { "j_k_0", j_k_0 },
                { "j_k_1", j_k_1 },
                { "j_k_2", j_k_2 },
                { "j_k_3", j_k_3 },
            };
            std::filesystem::create_directories(debug_dir);
            mrock::utility::saveString(debug_json.dump(4), debug_dir + "debug.json.gz");
        }

        return current_density;
    }

    std::string DiracSystem::info() const
    {
        return "DiracSystem\nE_F=" + std::to_string(E_F) + " * hbar omega_L"
            + "\nv_F=" + std::to_string(v_F) + " * pm / T_L"
            + "\nband_width=" + std::to_string(band_width)
            + "\nmax_k=" + std::to_string(max_k) + " pm"
            + "\nmax_kappa_compare=" + std::to_string(max_kappa_compare) + " pm^2";
    }

    h_float DiracSystem::kappa_integration_upper_limit(h_float k_z) const
    {
        const h_float k_z_squared = k_z * k_z;
        assert(max_kappa_compare >= k_z_squared);
        return sqrt(max_kappa_compare - k_z_squared);
    }
}