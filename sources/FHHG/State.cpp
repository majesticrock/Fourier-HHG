#include "State.hpp"

namespace FHHG {
    State::State(size_t n_frequencies)
        : n_frequencies(n_frequencies), state_vector(n_frequencies, state_type::Zero()) 
    { }
}