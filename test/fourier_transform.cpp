#include <catch2/catch_template_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <tnt/dsp/fourier_transform.hpp>
#include <tnt/dsp/signal_generator.hpp>

#include <complex>
#include <concepts>
#include <valarray>

#include <fftw3.h>

// TODO: Test fourier_transform/inverse_fourier_transform of larger sizes for speed

using namespace tnt;
using namespace Catch::Matchers;

// TODO: add double
TEMPLATE_TEST_CASE("fourier_transform", "[fourier_transform]", float)
{
    const std::size_t N = 128;
    const dsp::signal_generator<TestType> g(44100, N);

    SECTION("fourier transform of a complex signal")
    {
        std::vector<TestType> x_real = g.cosine(100);
        std::vector<TestType> x_imag = g.sine(100);

        fftw_complex* input = (fftw_complex*)fftwf_malloc(sizeof(fftw_complex) * N);
        fftw_complex* output = (fftw_complex*)fftwf_malloc(sizeof(fftw_complex) * N);
        fftw_plan plan = fftw_plan_dft_1d(N, input, output, FFTW_FORWARD, FFTW_ESTIMATE);

        std::vector<std::complex<TestType>> x(g.size());
        for (std::size_t n = 0; n < x.size(); ++n)
        {
            x[n] = { x_real[n], x_imag[n] };
            input[n][0] = x_real[n];
            input[n][1] = x_imag[n];
        }

        fftw_execute(plan);

        std::vector<std::complex<TestType>> X1(x.size());
        std::vector<std::complex<TestType>> X2(x.size());

        dsp::fourier_transform(x.data(), X1.data(), x.size());
        dft(x.data(), X2.data(), x.size());

        for (auto m = 0; m < x.size(); ++m)
        {
            //CHECK_THAT(X1[m].real(), WithinULP(X2[m].real(), 100));
            //CHECK_THAT(X1[m].imag(), WithinULP(X2[m].imag(), 100));
            CHECK_THAT(output[m][0], WithinULP(X1[m].real(), 100));
            CHECK_THAT(output[m][1], WithinULP(X1[m].imag(), 100));
        }

        fftw_destroy_plan(plan);
        fftw_free(input);
        fftw_free(output);
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