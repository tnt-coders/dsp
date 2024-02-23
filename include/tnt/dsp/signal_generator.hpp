#include <cmath>
#include <complex>
#include <concepts>
#include <vector>

namespace tnt::dsp
{

template <std::floating_point T>
class signal_generator final
{
public:
    signal_generator(std::size_t sample_rate, std::size_t size)
        : m_sample_rate(sample_rate)
        , m_sample_interval(1 / static_cast<T>(m_sample_rate))
        , m_size(size)
    {}

    std::vector<T> cosine(T frequency, T amplitude = 1, T phase_shift = 0, T vertical_shift = 0) const
    {
        std::vector<T> signal(m_size);
        for (size_t n = 0; n < signal.size(); ++n)
        {
            signal[n] = amplitude * std::cos(2 * static_cast<T>(M_PI) * frequency * n * m_sample_interval + phase_shift) + vertical_shift;
        }

        return signal;
    }

    std::vector<T> sine(T frequency, T amplitude = 1, T phase_shift = 0, T vertical_shift = 0) const
    {
        std::vector<T> signal(m_size);
        for (size_t n = 0; n < signal.size(); ++n)
        {
            signal[n] = amplitude * std::sin(2 * static_cast<T>(M_PI) * frequency * n * m_sample_interval + phase_shift) + vertical_shift;
        }

        return signal;
    }

private:
    std::size_t m_sample_rate;
    T m_sample_interval;
    std::size_t m_size;
};

}
