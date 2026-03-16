#pragma once

#include "../interval.h"
#include "util/math.hpp"
#include "vecmath.hpp"
#include "profiler.hpp"
#include "../raytracer.hpp"

#include <OpenImageIO/imageio.h>
#include <cmath>


using color = Vector3f;

inline Float linear_to_gamma(Float linear_component) {
    if (linear_component > 0) return std::sqrt(linear_component);
    return 0;
}

inline void write_color(std::vector<float> &out, const color &pixel_color, int index) {
    PROFILE_SCOPE("write_color");
    auto r = pixel_color[0];
    auto g = pixel_color[1];
    auto b = pixel_color[2];

    r = linear_to_gamma(r);
    g = linear_to_gamma(g);
    b = linear_to_gamma(b);

    static const interval intensity(0.000, 0.999);
    float rbyte = intensity.clamp(r);
    float gbyte = intensity.clamp(g);
    float bbyte = intensity.clamp(b);

    out[index + 0] = rbyte;
    out[index + 1] = gbyte;
    out[index + 2] = bbyte;
}


/************************************************************************************/
// new implementation


constexpr Float Lambda_min = 360, Lambda_max = 830;
static constexpr int nSpectrumSamples = 32;
static constexpr Float CIE_Y_integral = 106.856895;

class SampledSpectrum;
class SampledWavelengths;

// TODO: replace OOP approach
class Spectrum {
public:
    virtual ~Spectrum() = default;
    virtual Float operator()(Float lambda) const = 0;
    virtual Float MaxValue() const = 0;
    virtual SampledSpectrum Sample(const SampledWavelengths &lambda) const;
};

class ConstantSpectrum : public Spectrum {
public:
    ConstantSpectrum(Float c) : c(c) {}
    Float operator()(Float lambda) const { return c; }
    SampledSpectrum Sample(const SampledWavelengths &) const;
    Float MaxValue() const { return c; }

private:
    Float c;
};

class DenselySampledSpectrum {
public:
    DenselySampledSpectrum(int lambda_min = Lambda_min, int lambda_max = Lambda_max,
                           Allocator alloc = {})
        : lambda_min(lambda_min),
          lambda_max(lambda_max),
          values(lambda_max - lambda_min + 1, alloc) {}
private:
    int lambda_min, lambda_max;
    std::vector<Float> values;
};

class SampledWavelengths {
public:
    static SampledWavelengths sampleUniform(Float u, Float lambda_min = Lambda_min,
                                            Float lambda_max = Lambda_max) {
        SampledWavelengths swl;
        swl.lambda[0] = lerp(u, lambda_min, lambda_max);

        Float delta = (lambda_max - lambda_min) / nSpectrumSamples;
        for (int i = 1; i < nSpectrumSamples; ++i) {
            swl.lambda[i] = swl.lambda[i - 1] + delta;
            if (swl.lambda[i] > lambda_max)
                swl.lambda[i] = lambda_min + (swl.lambda[i] - lambda_max);
        }

        for (int i = 0; i < nSpectrumSamples; ++i)
            swl.pdf[i] = 1 / (lambda_max - lambda_min);

        return swl;
    }

    std::string toString() const;

    Float operator[](int i) const { return lambda[i]; }
    Float &operator[](int i) { return lambda[i]; }

private:
    std::array<Float, nSpectrumSamples> lambda, pdf;
};

class SampledSpectrum {
public:
    SampledSpectrum() = default;

    explicit SampledSpectrum(Float c) { values.fill(c); }

    SampledSpectrum(std::span<const Float> v) {
        DCHECK_EQ(NSpectrumSamples, v.size());
        for (int i = 0; i < nSpectrumSamples; ++i)
            values[i] = v[i];
    }

    Float operator[](int i) const {
        DCHECK(i >= 0 && i < nSpectrumSamples);
        return values[i];
    }
    Float &operator[](int i) {
        DCHECK(i >= 0 && i < nSpectrumSamples);
        return values[i];
    }

    explicit operator bool() const {
        for (int i = 0; i < nSpectrumSamples; ++i)
            if (values[i] != 0)
                return true;
        return false;
    }

private:
    std::array<Float, nSpectrumSamples> values;
};

class PiecewiseLinearSpectrum {
  public:
    PiecewiseLinearSpectrum() = default;

    void Scale(Float s) {
        for (Float &v : values)
            v *= s;
    }

    Float MaxValue() const;

    SampledSpectrum Sample(const SampledWavelengths &lambda) const {
        SampledSpectrum s;
        for (int i = 0; i < nSpectrumSamples; ++i)
            s[i] = (*this)(lambda[i]);
        return s;
    }
    Float operator()(Float lambda) const;

    std::string ToString() const;

    PiecewiseLinearSpectrum(std::span<const Float> lambdas,
                            std::span<const Float> values);

    static std::optional<Spectrum> Read(const std::string &filename);

    static PiecewiseLinearSpectrum *FromInterleaved(std::span<const Float> samples,
                                                    bool normalize=true);

  private:
    std::vector<Float> lambdas, values;
};



class RGB {
public:
    RGB() = default;
    RGB(Float r, Float g, Float b) : r(r), g(g), b(b) {}

    Float operator[](int c) const {
        DCHECK(c >= 0 && c < 3);
        if (c == 0)
            return r;
        else if (c == 1)
            return g;
        return b;
    }

    Float &operator[](int c) {
        DCHECK(c >= 0 && c < 3);
        if (c == 0)
            return r;
        else if (c == 1)
            return g;
        return b;
    }

    Float r = 0, g = 0, b = 0;
};


class XYZ {
public:
    XYZ() = default;
    XYZ(Float X, Float Y, Float Z) : X(X), Y(Y), Z(Z) {}

    Float X = 0, Y = 0, Z = 0;
};

class RGBSigmoidPolynomial {
public:
    RGBSigmoidPolynomial(Float c0, Float c1, Float c2) 
    : c0(c0), c1(c1), c2(c2) {}

    Float operator()(Float lambda) const {
        return s(c0*sqr(lambda) + c1*lambda + c2);
    }

private:
    static  Float s(Float x) {
        if (std::isinf(x)) return x > 0 ? 1 : 0;
        return .5 + x / (2 * std::sqrt(1+sqr(x)));
    }

    Float c0, c1, c2;
};


class RGBToSpectrumTable {
  public:
    static constexpr int res = 64;

    using CoefficientArray = float[3][res][res][res][3];

    RGBToSpectrumTable(const float *zNodes, const CoefficientArray *coeffs)
        : zNodes(zNodes), coeffs(coeffs) {}

    RGBSigmoidPolynomial operator()(RGB rgb) const;

    static const RGBToSpectrumTable *DCI_P3;

    std::string ToString() const;

  private:
    // RGBToSpectrumTable Private Members
    const float *zNodes;
    const CoefficientArray *coeffs;
};


template <typename Predicate>
inline size_t FindInterval(size_t sz, const Predicate &pred) {
    using ssize_t = std::make_signed_t<size_t>;
    ssize_t size = (ssize_t)sz - 2, first = 1;
    while (size > 0) {
        // Evaluate predicate at midpoint and update _first_ and _size_
        size_t half = (size_t)size >> 1, middle = first + half;
        bool predResult = pred(middle);
        first = predResult ? middle + 1 : first;
        size = predResult ? size - (half + 1) : half;
    }
    return (size_t)clamp((ssize_t)first - 1, 0, sz - 2);
}

RGBSigmoidPolynomial inline RGBToSpectrumTable::operator()(RGB rgb) const {
    DCHECK(rgb[0] >= 0.f && rgb[1] >= 0.f && rgb[2] >= 0.f && rgb[0] <= 1.f &&
           rgb[1] <= 1.f && rgb[2] <= 1.f);

    // Handle uniform _rgb_ values
    if (rgb[0] == rgb[1] && rgb[1] == rgb[2])
        return RGBSigmoidPolynomial(0, 0,
                                    (rgb[0] - .5f) / std::sqrt(rgb[0] * (1 - rgb[0])));

    // Find maximum component and compute remapped component values
    int maxc =
        (rgb[0] > rgb[1]) ? ((rgb[0] > rgb[2]) ? 0 : 2) : ((rgb[1] > rgb[2]) ? 1 : 2);
    float z = rgb[maxc];
    float x = rgb[(maxc + 1) % 3] * (res - 1) / z;
    float y = rgb[(maxc + 2) % 3] * (res - 1) / z;

    // Compute integer indices and offsets for coefficient interpolation
    int xi = std::min((int)x, res - 2), yi = std::min((int)y, res - 2),
        zi = FindInterval(res, [&](int i) { return zNodes[i] < z; });
    Float dx = x - xi, dy = y - yi, dz = (z - zNodes[zi]) / (zNodes[zi + 1] - zNodes[zi]);

    // Trilinearly interpolate sigmoid polynomial coefficients _c_
    std::array<Float, 3> c;
    for (int i = 0; i < 3; ++i) {
        // Define _co_ lambda for looking up sigmoid polynomial coefficients
        auto co = [&](int dx, int dy, int dz) {
            return (*coeffs)[maxc][zi + dz][yi + dy][xi + dx][i];
        };

        c[i] = lerp(dz,
                    lerp(dy, lerp(dx, co(0, 0, 0), co(1, 0, 0)),
                         lerp(dx, co(0, 1, 0), co(1, 1, 0))),
                    lerp(dy, lerp(dx, co(0, 0, 1), co(1, 0, 1)),
                         lerp(dx, co(0, 1, 1), co(1, 1, 1))));
    }

    return RGBSigmoidPolynomial(c[0], c[1], c[2]);
}


namespace Spectra {

DenselySampledSpectrum *x, *y, *z;

inline const DenselySampledSpectrum &X() {
    extern DenselySampledSpectrum *x;
    return *x;
}

inline const DenselySampledSpectrum &Y() {
    extern DenselySampledSpectrum *y;
    return *y;
}

inline const DenselySampledSpectrum &Z() {
    extern DenselySampledSpectrum *z;
    return *z;
}

}
