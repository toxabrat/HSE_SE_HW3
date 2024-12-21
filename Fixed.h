#pragma once
#include <iostream>
#include <ostream>
#include <cassert>
#include <random>

template<size_t N, size_t K, bool fast = false>
struct Fixed {
    using v_type = std::conditional_t<
        fast,
        std::conditional_t<
            N <= 8,  int_fast8_t,
            std::conditional_t<
                N <= 16, int_fast16_t,
                std::conditional_t<
                    N <= 32, int_fast32_t,
                    std::conditional_t<
                        N <= 64, int_fast64_t,
                        void 
                    >
                >
            >
        >,
        int64_t 
    >;

    static_assert(N <= std::numeric_limits<int64_t>::digits, "N must be less than or equal to 63");
    static_assert(K < N, "K must be less than N");

    static constexpr size_t type_bits = std::numeric_limits<v_type>::digits;

    static constexpr int64_t mask_N = (1LL << N) - 1;

    static constexpr int64_t sign_extension_shift_N = type_bits - N;

    static constexpr int64_t trim_to_n_bits(int64_t x) {
        return (x & mask_N) << sign_extension_shift_N >> sign_extension_shift_N;
    }

    static constexpr u_int64_t mask_K = (1LL << K) - 1;

    static constexpr u_int64_t trim_to_Uk_bits(u_int64_t x) {
        return (x & mask_K);
    }

    constexpr Fixed(int x) : v(trim_to_n_bits(static_cast<v_type>(x) << K)) {}
    constexpr Fixed(float f) : v(trim_to_n_bits(static_cast<v_type>(f * (1LL << K)))) {}
    constexpr Fixed(double f) : v(trim_to_n_bits(static_cast<v_type>(f * (1LL << K)))) {}
    constexpr Fixed() : v(0) {}
    template<typename T>
    constexpr Fixed(T f) : Fixed((double) f) {}

    static constexpr Fixed from_raw(v_type x) {
        Fixed ret;
        ret.v = trim_to_n_bits(x);
        return ret;
    }

    static constexpr Fixed get_inf() {
        return from_raw(mask_N); 
    }

    static constexpr Fixed get_eps() {
        return Fixed::from_raw(1);
    }

    v_type v;

    operator double() const { return static_cast<double>(v) / (1LL << K); }

    auto operator<=>(const Fixed& other) const {
        return v <=> other.v;
    }

    template <typename T>
    auto operator<=>(const T& other) const {
        return static_cast<double>(*this) <=> static_cast<double>(other);
    }

    bool operator==(const Fixed& other) const {
        return v == other.v;
    }

    template <typename T>
    bool operator==(const T& other) const {
        return static_cast<double>(*this) == static_cast<double>(other);
    }
};

template <std::size_t N, std::size_t K>
using FastFixed = Fixed<N, K, true>;

template<size_t N, size_t K, bool fast>
Fixed<N, K, fast> operator+(Fixed<N, K, fast> a, Fixed<N, K, fast> b) {
    return Fixed<N, K, fast>::from_raw(a.v + b.v);
}

template<size_t N, size_t K, bool fast>
Fixed<N, K, fast> operator-(Fixed<N, K, fast> a, Fixed<N, K, fast> b) {
    return Fixed<N, K, fast>::from_raw(a.v - b.v);
}

template<size_t N, size_t K, bool fast>
Fixed<N, K, fast> operator*(Fixed<N, K, fast> a, Fixed<N, K, fast> b) {
    return Fixed<N, K, fast>::from_raw((a.v * b.v) >> K);
}

template<size_t N, size_t K, bool fast>
Fixed<N, K, fast> operator/(Fixed<N, K, fast> a, Fixed<N, K, fast> b) {
    return Fixed<N, K, fast>::from_raw((a.v << K) / b.v);
}

template<size_t N, size_t K, bool fast>
Fixed<N, K, fast>& operator+=(Fixed<N, K, fast>& a, Fixed<N, K, fast> b) {
    return a = a + b;
}

template<size_t N, size_t K, bool fast>
Fixed<N, K, fast>& operator-=(Fixed<N, K, fast>& a, Fixed<N, K, fast> b) {
    return a = a - b;
}

template<size_t N, size_t K, bool fast>
Fixed<N, K, fast>& operator*=(Fixed<N, K, fast>& a, Fixed<N, K, fast> b) {
    return a = a * b;
}

template<size_t N, size_t K, bool fast>
Fixed<N, K, fast>& operator/=(Fixed<N, K, fast>& a, Fixed<N, K, fast> b) {
    return a = a / b;
}

template<size_t N, size_t K, bool fast>
Fixed<N, K, fast> operator-(Fixed<N, K, fast> x) {
    return Fixed<N, K, fast>::from_raw(-x.v);
}


template<typename T, size_t N, size_t K, bool fast>
Fixed<N, K, fast> operator+(Fixed<N, K, fast> a, T b) {
    return a + Fixed<N, K, fast>(static_cast<double>(b));
}

template<typename T, size_t N, size_t K, bool fast>
Fixed<N, K, fast> operator-(Fixed<N, K, fast> a, T b) {
    return a - Fixed<N, K, fast>(static_cast<double>(b));
}

template<typename T, size_t N, size_t K, bool fast>
Fixed<N, K, fast> operator*(Fixed<N, K, fast> a, T b) {
    return a * Fixed<N, K, fast>(static_cast<double>(b));
}

template<typename T, size_t N, size_t K, bool fast>
Fixed<N, K, fast> operator/(Fixed<N, K, fast> a, T b) {
    return a / Fixed<N, K, fast>(static_cast<double>(b));
}

template<typename T, size_t N, size_t K, bool fast>
Fixed<N, K, fast>& operator+=(Fixed<N, K, fast>& a, T b) {
    return a = a + b;
}

template<typename T, size_t N, size_t K, bool fast>
Fixed<N, K, fast>& operator-=(Fixed<N, K, fast>& a, T b) {
    return a = a - b;
}

template<typename T, size_t N, size_t K, bool fast>
Fixed<N, K, fast>& operator*=(Fixed<N, K, fast>& a, T b) {
    return a = a * b;
}

template<typename T, size_t N, size_t K, bool fast>
Fixed<N, K, fast>& operator/=(Fixed<N, K, fast>& a, T b) {
    return a = a / b;
}

template<size_t N, size_t K, bool fast>
Fixed<N, K, fast> abs(Fixed<N, K, fast> x) {
    return Fixed<N, K, fast>::from_raw(x.v < 0 ? -x.v : x.v);
}

template<size_t N, size_t K, bool fast>
std::ostream& operator<<(std::ostream& out, Fixed<N, K, fast> x) {
    return out << static_cast<double>(x);
}
