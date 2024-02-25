#include "impl/utility.hpp"

#include <cstddef>

namespace tnt::dsp::impl {

// Returns true if x is a power of 2
bool is_power_of_2(std::size_t N)
{
    return N && (!(N & (N - 1)));
}

}