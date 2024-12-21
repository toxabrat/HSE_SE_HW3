
#include <random>
#include <limits>
#include <type_traits>
#include "Fixed.h"

static std::mt19937_64 rnd(1337);

class RandomGenerator {
    
    public:
        template <size_t N, size_t K, bool fast>
        static Fixed<N, K> random01(const Fixed<N, K, fast>&) {
            static std::uniform_int_distribution<> distrib(0, (1LL << K) - 1);
            return Fixed<N, K>::from_raw(distrib(rnd));
        }

        static float random01(const float&) {  
            static std::uniform_real_distribution<float> distrib(0.0f, 1.0f);
            return distrib(rnd);
        }

        static double random01(const double&) { 
            static std::uniform_real_distribution<double> distrib(0.0, 1.0);
            return distrib(rnd);
        }

        template <typename T>
        static T random01(const T&) {
            static std::uniform_real_distribution<T> distrib(0.0, 1.0);
            return distrib(rnd);
        }

        template <typename P>
        static auto random01() {
            return random01(P{});
        }

};