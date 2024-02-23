#include <catch2/catch_approx.hpp>
#include <catch2/catch_template_test_macros.hpp>
#include <tnt/dsp/fourier_transform.hpp>
#include <tnt/dsp/signal_generator.hpp>

#include <complex>
#include <concepts>
#include <valarray>

// TODO: Test fourier_transform/inverse_fourier_transform of larger sizes for speed

using namespace tnt;

// TODO: add double
TEMPLATE_TEST_CASE("fourier_transform", "[fourier_transform]", float)
{
    const dsp::signal_generator<TestType> g(44100, 8);

    SECTION("fourier transform of a complex signal")
    {
        std::vector<TestType> x = g.cosine(100);
        std::vector<TestType> x_imag = g.sine(100);

        /*std::vector<std::complex<TestType>> x(128);
        for (std::size_t n = 0; n < x.size(); ++n)
        {
            x[n] = { x_real[n], x_imag[n] };
        }*/

        std::vector<std::complex<TestType>> X1(x.size());
        std::vector<std::complex<TestType>> X2(x.size());

        dsp::fourier_transform(x.data(), X1.data(), x.size());
        dft(x.data(), X2.data(), x.size());

        for (auto m = 0; m < x.size(); ++m)
        {
            CHECK(X1[m].real() == Catch::Approx(X2[m].real()).scale(1));
            CHECK(X1[m].imag() == Catch::Approx(X2[m].imag()).scale(1));
        }
    }
}

// Implementation of slow Fourier transform to compare against
template <std::floating_point T>
void dft(const T* x, std::complex<T>* X, std::size_t N)
{
    // Take advantage of DFT symmetry when dealing with real input signals
    // Only the first N/2 + 1 outputs are unique
    for (size_t k = 0; k < N / 2 + 1; ++k)
    {
        for (size_t n = 0; n < N; ++n)
        {
            X[k] += x[n] * std::polar(static_cast<T>(1), -2 * static_cast<T>(M_PI) * n * k / N);
        }

        // X(N-k) = X(k)* for k = 1 -> N/2
        if (k != 0)
        {
            X[N - k] = std::conj(X[k]);
        }
    }
}

// Implementation of slow fourier transform to compare against
template <std::floating_point T>
void dft(const std::complex<T>* x, std::complex<T>* X, std::size_t N)
{
    for (size_t k = 0; k < N; ++k)
    {
        for (size_t n = 0; n < N; ++n)
        {
            X[k] += x[n] * std::polar(static_cast<T>(1), -2 * static_cast<T>(M_PI) * n * k / N);
        }
    }
}
