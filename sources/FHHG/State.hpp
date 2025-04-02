#pragma once
#include "GlobalDefinitions.hpp"
#include <vector>

namespace FHHG {
    struct State {
        static constexpr size_t D = 3; ///< Dimension of the state vector
        const size_t n_frequencies{}; ///< Number of frequencies

        using state_type = Eigen::Vector<h_complex, D>; ///< Type of the state vector
        std::vector<state_type> state_vector{}; ///< State vectors per frequency

        State(size_t n_frequencies);

        inline size_t size() const noexcept { return D * n_frequencies; }
        inline h_complex& operator[](size_t i) noexcept { return state_vector[i / D](i % D); }
        inline const h_complex& operator[](size_t i) const noexcept { return state_vector[i / D](i % D); }

        inline h_complex& operator()(size_t i, size_t j) noexcept { return state_vector[i](j); }
        inline const h_complex& operator()(size_t i, size_t j) const noexcept { return state_vector[i](j); }

        inline state_type& operator()(size_t i) noexcept { return state_vector[i]; }
        inline const state_type& operator()(size_t i) const noexcept { return state_vector[i]; }
    };
}