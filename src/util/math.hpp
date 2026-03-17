#pragma once

#include <array>
#include <cmath>
#include <span>
#include <type_traits>
#include "../raytracer.hpp"
#include "check.h"

constexpr Float ShadowEpsilon = 0.0001f;
constexpr Float Pi            = 3.14159265358979323846;
constexpr Float InvPi         = 0.31830988618379067154;
constexpr Float Inv2Pi        = 0.15915494309189533577;
constexpr Float Inv4Pi        = 0.07957747154594766788;
constexpr Float PiOver2       = 1.57079632679489661923;
constexpr Float PiOver4       = 0.78539816339744830961;
constexpr Float Sqrt2         = 1.41421356237309504880;
constexpr Float infinity      = std::numeric_limits<Float>::infinity();

template <typename T>
inline constexpr T sqr(T x) { return x*x; }

template <typename T, typename U, typename V>
inline constexpr T clamp(T val, U low, V high) {
    if (val < low)
        return T(low);
    else if (val > high)
        return T(high);
    else
        return val;
}

constexpr inline Float radians(Float deg) { return (Pi / 180) * deg; }
constexpr inline Float degrees(Float rad) { return (180 / Pi) * rad; }

inline constexpr Float safeASin(Float x) { return std::asin(clamp(x, -1, 1)); }
inline constexpr Float safeACos(Float x) { return std::acos(clamp(x, -1, 1)); }

template <typename Ta, typename Tb, typename Tc, typename Td>
inline auto differenceOfProducts(Ta a, Tb b, Tc c, Td d) {
    auto cd = c * d;
    auto differenceOfProducts = std::fma(a, b, -cd);
    auto error = std::fma(-c, d, cd);
    return differenceOfProducts + error;
}

template <typename Ta, typename Tb, typename Tc, typename Td>
inline auto sumOfProducts(Ta a, Tb b, Tc c, Td d) {
    auto cd = c * d;
    auto sumOfProducts = std::fma(a, b, cd);
    auto error = std::fma(c, d, -cd);
    return sumOfProducts + error;
}

inline Float innerProduct(Float a, Float b) {
    return a * b;
}

template <typename... T>
inline Float innerProduct(Float a, Float b, T... terms) {
    Float ab = a * b;
    Float tp = innerProduct(terms...);
    Float sum = ab + tp;
    return sum;
}

inline Float lerp(Float x, Float a, Float b) {
    return (1 - x) * a + x * b;
}

template <int n>
using SquareMatrixType = std::array<std::array<Float, n>, n>;

template <int N, class Mat>
inline void init(Mat&, int, int) {}

template <int N, class Mat, class V, class... Args>
inline void init(Mat& m, int i, int j, V v, Args... args) {
    static_assert(std::is_convertible_v<V, Float>,
                  "SquareMatrix init values must be convertible to Float");

    m[i][j] = static_cast<Float>(v);
    if (++j == N) { j = 0; ++i; }
    init<N>(m, i, j, args...);
}

template <typename Predicate>
inline size_t findInterval(size_t sz, const Predicate &pred) {
    using ssize_t = std::make_signed_t<size_t>;
    ssize_t size = (ssize_t)sz - 2, first = 1;
    while (size > 0) {
        size_t half = (size_t)size >> 1, middle = first + half;
        bool predResult = pred(middle);
        first = predResult ? middle + 1 : first;
        size = predResult ? size - (half + 1) : half;
    }
    return (size_t)clamp((ssize_t)first - 1, 0, sz - 2);
}

template <int N>
inline void initDiag(SquareMatrixType<N > &m, int i) {}

template <int N, typename... Args>
inline void initDiag(SquareMatrixType<N> &m, int i, Float v, Args... args) {
    m[i][i] = v;
    initDiag<N>(m, i + 1, args...);
}

template <int N>
class SquareMatrix {
public:
    SquareMatrix() {
        for (int i = 0; i < N; ++i)
            for (int j = 0; j < N; ++j)
                m[i][j] = (i == j) ? Float(1) : Float(0);
    }

    explicit SquareMatrix(const Float mat[N][N]) {
        for (int i = 0; i < N; ++i)
            for (int j = 0; j < N; ++j)
                m[i][j] = mat[i][j];
    }

    explicit SquareMatrix(const SquareMatrixType<N> &mat) : m(mat) {}

    template <typename... Args>
    requires (sizeof...(Args) == N * N) && (std::conjunction_v<std::is_convertible<Args, Float>...>)
    explicit SquareMatrix(Args... args) {
        init<N>(m, 0, 0, args...);
    }

    static SquareMatrix zero() {
        SquareMatrix r;
        for (int i = 0; i < N; ++i)
            for (int j = 0; j < N; ++j)
                r.m[i][j] = Float(0);
        return r;
    }

    SquareMatrix operator+(const SquareMatrix& o) const {
        SquareMatrix r = *this;
        for (int i = 0; i < N; ++i)
            for (int j = 0; j < N; ++j)
                r.m[i][j] += o.m[i][j];
        return r;
    }

    SquareMatrix& operator+=(const SquareMatrix& o) {
        for (int i = 0; i < N; ++i)
            for (int j = 0; j < N; ++j)
                m[i][j] += o.m[i][j];
        return *this;
    }

    SquareMatrix operator*(Float s) const {
        SquareMatrix r = *this;
        for (int i = 0; i < N; ++i)
            for (int j = 0; j < N; ++j)
                r.m[i][j] *= s;
        return r;
    }

    SquareMatrix& operator*=(Float s) {
        for (int i = 0; i < N; ++i)
            for (int j = 0; j < N; ++j)
                m[i][j] *= s;
        return *this;
    }

    SquareMatrix operator/(Float s) const {
        DCHECK_NE(s, 0);
        SquareMatrix r = *this;
        for (int i = 0; i < N; ++i)
            for (int j = 0; j < N; ++j)
                r.m[i][j] /= s;
        return r;
    }

    SquareMatrix& operator/=(Float s) {
        DCHECK_NE(s, 0);
        for (int i = 0; i < N; ++i)
            for (int j = 0; j < N; ++j)
                m[i][j] /= s;
        return *this;
    }

    bool operator==(const SquareMatrix& o) const {
        for (int i = 0; i < N; ++i)
            for (int j = 0; j < N; ++j)
                if (m[i][j] != o.m[i][j])
                    return false;
        return true;
    }

    bool operator!=(const SquareMatrix& o) const { return !(*this == o); }

    bool operator<(const SquareMatrix& o) const {
        for (int i = 0; i < N; ++i)
            for (int j = 0; j < N; ++j) {
                if (m[i][j] < o.m[i][j]) return true;
                if (m[i][j] > o.m[i][j]) return false;
            }
        return false;
    }

    std::span<const Float, N> operator[](int i) const {
        DCHECK(i >= 0 && i < N);
        return std::span<const Float, N>(m[i].data(), N);
    }

    std::span<Float, N> operator[](int i) {
        DCHECK(i >= 0 && i < N);
        return std::span<Float, N>(m[i].data(), N);
    }

    template <typename... Args>
    static SquareMatrix diag(Float v, Args... args) {
        CHECK(1 + sizeof...(Args) == N);
        SquareMatrix m;
        initDiag<N>(m.m, 0, v, args...);
        return m;
    }

    const SquareMatrixType<N>& getMatrix() const { return m; }
    SquareMatrixType<N>& getMatrix() { return m; }

    bool isIdentity() const;

    std::string toString() const {
        std::string s;
        s += "[\n";
        for (int i = 0; i < N; ++i) {
            s += "  [ ";
            for (int j = 0; j < N; ++j) {
                s += std::format("{}", m[i][j]);
                if (j + 1 < N) s += ", ";
            }
            s += " ]";
            if (i + 1 < N) s += ",";
            s += "\n";
        }
        s += "]";
        return s;
    }

private:
    SquareMatrixType<N> m{};
};


template <typename Tresult, int N, typename T>
inline Tresult mul(const SquareMatrix<N> &m, const T &v) {
    Tresult result;
    for (int i = 0; i < N; ++i) {
        result[i] = 0;
        for (int j = 0; j < N; ++j)
            result[i] += m[i][j] * v[j];
    }
    return result;
}


template <int N, typename T>
inline T operator*(const SquareMatrix<N> &m, const T &v) {
    return mul<T>(m, v);
}

template <int N>
inline SquareMatrix<N> operator*(Float s, const SquareMatrix<N> &M) {
    return M * s;
}

inline SquareMatrix<4> operator*(const SquareMatrix<4> &m1,
                                 const SquareMatrix<4> &m2) {
    SquareMatrix<4> r;
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j)
            r[i][j] = innerProduct(m1[i][0], m2[0][j], m1[i][1], m2[1][j], m1[i][2],
                                   m2[2][j], m1[i][3], m2[3][j]);
    return r;
}

inline SquareMatrix<3> operator*(const SquareMatrix<3> &m1,
                                 const SquareMatrix<3> &m2) {
    SquareMatrix<3> r;
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            r[i][j] =
                innerProduct(m1[i][0], m2[0][j], m1[i][1], m2[1][j], m1[i][2], m2[2][j]);
    return r;
}

template <int N>
inline SquareMatrix<N> operator*(const SquareMatrix<N> &m1,
                                 const SquareMatrix<N> &m2) {
    SquareMatrix<N> r;
    for (int i = 0; i < N; ++i)
        for (int j = 0; j < N; ++j) {
            r[i][j] = 0;
            for (int k = 0; k < N; ++k)
                r[i][j] = FMA(m1[i][k], m2[k][j], r[i][j]);
        }
    return r;
}

template <int N>
inline bool SquareMatrix<N>::isIdentity() const {
    for (int i = 0; i < N; ++i)
        for (int j = 0; j < N; ++j) {
            if (i == j) {
                if (m[i][j] != 1)
                    return false;
            } else if (m[i][j] != 0)
                return false;
        }
    return true;
}

inline Float determinant(const SquareMatrix<3> &m) {
    Float minor12 = differenceOfProducts(m[1][1], m[2][2], m[1][2], m[2][1]);
    Float minor02 = differenceOfProducts(m[1][0], m[2][2], m[1][2], m[2][0]);
    Float minor01 = differenceOfProducts(m[1][0], m[2][1], m[1][1], m[2][0]);
    return std::fma(m[0][2], minor01,
               differenceOfProducts(m[0][0], minor12, m[0][1], minor02));
}

inline std::optional<SquareMatrix<3>> inverse(const SquareMatrix<3> &m) {
    Float det = determinant(m);
    if (det == 0)
        return {};
    Float invDet = 1 / det;

    SquareMatrix<3> r;

    r[0][0] = invDet * differenceOfProducts(m[1][1], m[2][2], m[1][2], m[2][1]);
    r[1][0] = invDet * differenceOfProducts(m[1][2], m[2][0], m[1][0], m[2][2]);
    r[2][0] = invDet * differenceOfProducts(m[1][0], m[2][1], m[1][1], m[2][0]);
    r[0][1] = invDet * differenceOfProducts(m[0][2], m[2][1], m[0][1], m[2][2]);
    r[1][1] = invDet * differenceOfProducts(m[0][0], m[2][2], m[0][2], m[2][0]);
    r[2][1] = invDet * differenceOfProducts(m[0][1], m[2][0], m[0][0], m[2][1]);
    r[0][2] = invDet * differenceOfProducts(m[0][1], m[1][2], m[0][2], m[1][1]);
    r[1][2] = invDet * differenceOfProducts(m[0][2], m[1][0], m[0][0], m[1][2]);
    r[2][2] = invDet * differenceOfProducts(m[0][0], m[1][1], m[0][1], m[1][0]);

    return r;
}
// ---- 4x4 inverse (your implementation, adapted to this SquareMatrix) ----
inline std::optional<SquareMatrix<4>> inverse(const SquareMatrix<4>& m) {
    // Via: https://github.com/google/ion/blob/master/ion/math/matrixutils.cc,
    // (c) Google, Apache license.

    Float s0 = differenceOfProducts(m[0][0], m[1][1], m[1][0], m[0][1]);
    Float s1 = differenceOfProducts(m[0][0], m[1][2], m[1][0], m[0][2]);
    Float s2 = differenceOfProducts(m[0][0], m[1][3], m[1][0], m[0][3]);

    Float s3 = differenceOfProducts(m[0][1], m[1][2], m[1][1], m[0][2]);
    Float s4 = differenceOfProducts(m[0][1], m[1][3], m[1][1], m[0][3]);
    Float s5 = differenceOfProducts(m[0][2], m[1][3], m[1][2], m[0][3]);

    Float c0 = differenceOfProducts(m[2][0], m[3][1], m[3][0], m[2][1]);
    Float c1 = differenceOfProducts(m[2][0], m[3][2], m[3][0], m[2][2]);
    Float c2 = differenceOfProducts(m[2][0], m[3][3], m[3][0], m[2][3]);

    Float c3 = differenceOfProducts(m[2][1], m[3][2], m[3][1], m[2][2]);
    Float c4 = differenceOfProducts(m[2][1], m[3][3], m[3][1], m[2][3]);
    Float c5 = differenceOfProducts(m[2][2], m[3][3], m[3][2], m[2][3]);

    Float determinant =
        innerProduct(s0, c5, -s1, c4, s2, c3, s3, c2, s5, c0, -s4, c1);

    if (determinant == 0)
        return {};

    Float s = 1 / determinant;

    Float inv[4][4] = {
        {s * innerProduct(m[1][1], c5, m[1][3], c3, -m[1][2], c4),
         s * innerProduct(-m[0][1], c5, m[0][2], c4, -m[0][3], c3),
         s * innerProduct(m[3][1], s5, m[3][3], s3, -m[3][2], s4),
         s * innerProduct(-m[2][1], s5, m[2][2], s4, -m[2][3], s3)},

        {s * innerProduct(-m[1][0], c5, m[1][2], c2, -m[1][3], c1),
         s * innerProduct(m[0][0], c5, m[0][3], c1, -m[0][2], c2),
         s * innerProduct(-m[3][0], s5, m[3][2], s2, -m[3][3], s1),
         s * innerProduct(m[2][0], s5, m[2][3], s1, -m[2][2], s2)},

        {s * innerProduct(m[1][0], c4, m[1][3], c0, -m[1][1], c2),
         s * innerProduct(-m[0][0], c4, m[0][1], c2, -m[0][3], c0),
         s * innerProduct(m[3][0], s4, m[3][3], s0, -m[3][1], s2),
         s * innerProduct(-m[2][0], s4, m[2][1], s2, -m[2][3], s0)},

        {s * innerProduct(-m[1][0], c3, m[1][1], c1, -m[1][2], c0),
         s * innerProduct(m[0][0], c3, m[0][2], c0, -m[0][1], c1),
         s * innerProduct(-m[3][0], s3, m[3][1], s1, -m[3][2], s0),
         s * innerProduct(m[2][0], s3, m[2][2], s0, -m[2][1], s1)}
    };

    return SquareMatrix<4>(inv);
}

template <int N>
inline SquareMatrix<N> transpose(const SquareMatrix<N> &m) {
    SquareMatrix<N> r;
    for (int i = 0; i < N; ++i)
        for (int j = 0; j < N; ++j)
            r[i][j] = m[j][i];
    return r;
}

template <int N>
SquareMatrix<N> invertOrExit(const SquareMatrix<N> &m) {
    std::optional<SquareMatrix<N>> inv = inverse(m);
    CHECK(inv.has_value());
    return *inv;
}
