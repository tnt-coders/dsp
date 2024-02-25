#pragma once

#include "impl/fourier_transform.hpp"
#include "impl/utility.hpp"

#include <algorithm>
#include <cmath>
#include <complex>
#include <concepts>
#include <cstddef>

// TODO: Assert input size is a power of 2

namespace tnt::dsp {

template <std::floating_point T>
void fourier_transform(const T* x, std::complex<T>* X, std::size_t N)
{
    assert(impl::is_power_of_2(N));

    // Taking advantage of symmetry the FFT of a real signal can be computed
    // using a single N/2-point complex FFT. Split the input signal into its
    // even and odd components and load the data into a single complex vector.
    for (std::size_t n = 0; n < N / 2; ++n)
    {
        X[n] = { x[n * 2], x[n * 2 + 1] };
    }

    // Perform an in-place FFT
    impl::in_place_fourier_transform(X, N / 2);

    // Save a copy of the original X[0]
    // It is needed later, but will be overwritten by the in-place algorithm
    std::complex<T> X_0 = X[0];

    // The FFT is periodic so it is valid append X[0] to the end. This is
    // required to avoid a buffer overflow in the next section.
    X[N / 2] = X_0;

    // Extract the real FFT from the output of the complex FFT
    for (std::size_t m = 0; m < N / 2; ++m)
    {
        T X_r_first = (X[m].real() + X[N / 2 - m].real()) / 2;
        T X_r_second = (X[m].real() - X[N / 2 - m].real()) / 2;

        T X_i_first = (X[m].imag() + X[N / 2 - m].imag()) / 2;
        T X_i_second = (X[m].imag() - X[N / 2 - m].imag()) / 2;

        T a = std::cos(static_cast<T>(M_PI) * m / (N / 2));
        T b = std::sin(static_cast<T>(M_PI) * m / (N / 2));

        T real = X_r_first + a * X_i_first - b * X_r_second;
        T imag = X_i_second - b * X_i_first - a * X_r_second;

        // Store in 2nd half of the array to avoid clobbering data
        X[m + N / 2] = { real, imag };
    }

    // Copy the 2nd half of the array into the first half
    std::copy(X + N / 2, X + N, X);

    // X[m] = X[N - m] where 1 <= m <= N / 2 - 1
    for (std::size_t m = 1; m < N / 2; ++m)
    {
        X[N - m] = std::conj(X[m]);  // (>*.*)> symmetry! <(*.*<)
    }

    // X[N / 2] is a special case
    X[N / 2] = X_0.real() - X_0.imag();
}

template <std::floating_point T>
void fourier_transform(const std::complex<T>* x, std::complex<T>* X, std::size_t N)
{
    assert(impl::is_power_of_2(N));

    // Copy the input data into the output array
    std::copy(x, x + N, X);

    // Perform an in-place FFT
    impl::in_place_fourier_transform(X, N);
}

}