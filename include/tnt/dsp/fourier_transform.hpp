#pragma once

#include <algorithm>
#include <cmath>
#include <complex>
#include <concepts>

// TODO: Assert input size is a power of 2

namespace tnt::dsp {

template <std::floating_point T>
void in_place_fourier_transform(std::complex<T>* x, std::size_t N);

template <std::floating_point T>
void fourier_transform(const T* x, std::complex<T>* X, std::size_t N)
{
    // Taking advantage of symmetry the FFT of a real signal can be computed
    // using a single N/2-point complex FFT. Split the input signal into its
    // even and odd components and load the data into a single complex vector.
    for (std::size_t n = 0; n < N / 2; ++n)
    {
        X[n] = { x[n * 2], x[n * 2 + 1] };
    }

    // Perform an in-place FFT
    in_place_fourier_transform(X, N / 2);

    // The FFT is periodic so it is valid append X[0] to the end. This is
    // required to avoid a buffer overflow in the next section.
    X[N / 2] = X[0];

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

    std::complex<T> temp = X[0];

    // Copy the 2nd half of the array into the first half
    std::copy(X + N / 2, X + N, X);

    // X[m] = X[N - m] where 1 <= m <= N / 2 - 1
    for (std::size_t m = 1; m < N / 2; ++m)
    {
        X[N - m] = std::conj(X[m]);  // (>*.*)> symmetry! <(*.*<)
    }

    // X[N / 2] is a special case
    X[N / 2] = temp.real() - temp.imag();
}

template <std::floating_point T>
void fourier_transform(const std::complex<T>* x, std::complex<T>* X, std::size_t N)
{
    // Copy the input data into the output array
    std::copy(x, x + N, X);

    // Perform an in-place FFT
    in_place_fourier_transform(X, N);
}

template <std::floating_point T>
void in_place_fourier_transform(std::complex<T>* x, std::size_t N) {
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