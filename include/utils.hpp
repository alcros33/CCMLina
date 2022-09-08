#include <chrono>

class PerfCounter
{
public:
    inline void tick() { m_then = std::chrono::high_resolution_clock::now(); }
    inline auto tock()
    {
        std::chrono::duration<double> diff = std::chrono::high_resolution_clock::now() - m_then;
        return diff.count();
    }

private:
    decltype(std::chrono::high_resolution_clock::now()) m_then;
};