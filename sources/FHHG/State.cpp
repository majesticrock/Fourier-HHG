#include "State.hpp"
#include <numeric>

namespace FHHG {
    State::State(size_t n_frequencies)
        : n_frequencies(n_frequencies), state_vector(n_frequencies, state_type::Zero()) 
    { }
    h_float State::relative_difference(const State &other) const
    {
        h_float diff{};
        for(size_t i = 0; i < n_frequencies; ++i) {
            diff += (state_vector[i] - other.state_vector[i]).norm() / (is_zero(state_vector[i].norm()) ? 1.0 : state_vector[i].norm());
        }
        return diff;
    }

    h_float State::absolut_difference(const State &other) const
    {
        h_float diff{};
        for(size_t i = 0; i < n_frequencies; ++i) {
            diff += (state_vector[i] - other.state_vector[i]).norm();
        }
        return diff;
    }

    bool State::is_close_relative(const State &other, h_float tolerance) const
    {
        h_float diff_norm{};
        for (size_t i = 0U; i < n_frequencies; ++i) {
            diff_norm = (other.state_vector[i] - state_vector[i]).norm();
            if (is_zero(other(i).norm())) {
                if (!is_zero(diff_norm)) return false;
            }
            else if (diff_norm > tolerance * other.state_vector[i].norm()) {
                return false;
            }
        }
        return true;
    }

    bool State::is_close_absolute(const State &other, h_float tolerance) const
    {
        return std::abs(absolut_difference(other)) < tolerance;
    }
    
    void State::setZero() noexcept
    {
        for (auto& state : state_vector) {
            state.setZero();
        }
    }
}