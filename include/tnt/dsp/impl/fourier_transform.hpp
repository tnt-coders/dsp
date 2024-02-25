#pragma once

#include "utility.hpp"

#include <cmath>
#include <complex>
#include <concepts>
#include <cstddef>

namespace tnt::dsp::impl {

template <std::floating_point T>
void in_place_fourier_transform(std::complex<T>* x, std::size_t N)
{
    assert(impl::is_power_of_2(N));

    // Bit-reversal permutation
    std::size_t log2_N = static_cast<std::size_t>(std::log2(N));
    for (std::size_t i = 0; i < N; ++i)
    {
        std::size_t j = 0;
        for (std::size_t k = 0; k < log2_N; ++k)
        {
            j |= ((i >> k) & 1) << (log2_N - 1 - k);
        }
        if (j > i)
        {
            std::swap(x[i], x[j]);
        }
    }

    // Cooley-Tukey iteration
    for (std::size_t size = 2; size <= N; size *= 2)
    {
        T angle = -2 * static_cast<T>(M_PI) / size;
        std::complex<T> w(1.0, 0.0);
        std::complex<T> wn(std::cos(angle), std::sin(angle));

        for (std::size_t j = 0; j < N; j += size)
        {
            std::complex<T> w_current(1.0, 0.0);
            for (std::size_t k = 0; k < size / 2; ++k)
            {
                std::complex<T> t = w_current * x[j + k + size / 2];
                std::complex<T> u = x[j + k];
                x[j + k] = u + t;
                x[j + k + size / 2] = u - t;
                w_current *= wn;
            }
            w *= wn;
        }
    }
}

}