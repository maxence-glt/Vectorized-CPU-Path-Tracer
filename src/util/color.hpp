#pragma once

#include "../interval.h"
#include "util/math.hpp"
#include "vecmath.hpp"
#include "profiler.hpp"
#include "../raytracer.hpp"

#include <OpenImageIO/imageio.h>
#include <cmath>

void spectrumInit();

using color = Vector3f;

inline Float linear_to_gamma(Float linear_component) {
    if (linear_component > 0) return std::sqrt(linear_component);
    return 0;
}

inline void write_color(std::vector<float> &out, const color &pixel_color, int index) {
    PROFILE_SCOPE("write_color");

    static const interval intensity(0.000, 0.999);

    out[index + 0] = intensity.clamp(pixel_color[0]);
    out[index + 1] = intensity.clamp(pixel_color[1]);
    out[index + 2] = intensity.clamp(pixel_color[2]);
}


/************************************************************************************/
// new implementation


constexpr Float Lambda_min = 360, Lambda_max = 830;
static constexpr int nSpectrumSamples = 32;
static constexpr Float CIE_Y_integral = 106.856895;

class SampledSpectrum;
class SampledWavelengths;

// TODO: replace OOP approach

class SampledSpectrum {
public:
    SampledSpectrum() = default;

    explicit SampledSpectrum(Float c) { values.fill(c); }

    SampledSpectrum(std::span<const Float> v) {
        DCHECK_EQ(nSpectrumSamples, v.size());
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

    SampledSpectrum &operator*=(const SampledSpectrum &s) {
        for (int i = 0; i < nSpectrumSamples; ++i)
            values[i] *= s.values[i];
        return *this;
    }
    SampledSpectrum operator*(const SampledSpectrum &s) const {
        SampledSpectrum ret = *this;
        return ret *= s;
    }

private:
    std::array<Float, nSpectrumSamples> values;
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

    SampledSpectrum PDF() const { return SampledSpectrum(pdf); }

    std::string toString() const;

    Float operator[](int i) const { return lambda[i]; }
    Float &operator[](int i) { return lambda[i]; }

private:
    std::array<Float, nSpectrumSamples> lambda, pdf;
};


class Spectrum {
public:
    virtual ~Spectrum() = default;
    virtual Float operator()(Float lambda) const = 0;
    virtual Float maxValue() const = 0;
    virtual SampledSpectrum sample(const SampledWavelengths &lambda) const = 0;
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

class PiecewiseLinearSpectrum : public Spectrum {
public:
    PiecewiseLinearSpectrum() = default;
    PiecewiseLinearSpectrum(std::span<const Float> l,
                            std::span<const Float> v)
    : lambdas(l.begin(), l.end()), values(v.begin(), v.end()) 
    {/* empty ctor */}

    PiecewiseLinearSpectrum &operator=(const PiecewiseLinearSpectrum &other) = default;

    void Scale(Float s) {
        for (Float &v : values)
            v *= s;
    }

    Float operator()(Float lambda) const override {
        if (lambdas.empty() || lambda < lambdas.front() || lambda > lambdas.back())
            return 0;

        int o = findInterval(lambdas.size(), [&](int i) { return lambdas[i] <= lambda; });
        DCHECK(lambda >= lambdas[o] && lambda <= lambdas[o + 1]);
        Float t = (lambda - lambdas[o]) / (lambdas[o + 1] - lambdas[o]);
        return lerp(t, values[o], values[o + 1]);
    }

    SampledSpectrum sample(const SampledWavelengths &lambda) const override {
        SampledSpectrum s;
        for (int i = 0; i < nSpectrumSamples; ++i)
            s[i] = (*this)(lambda[i]);
        return s;
    }

    Float maxValue() const override {
        if (values.empty())
            return 0;
        return *std::max_element(values.begin(), values.end());
    }

    std::string ToString() const;

    static std::optional<Spectrum> Read(const std::string &filename);

    static PiecewiseLinearSpectrum *fromInterleaved(std::span<const Float> samples,
                                                    bool normalize=true);

private:
    std::vector<Float> lambdas, values;
};

class DenselySampledSpectrum : public Spectrum {
public:
    DenselySampledSpectrum() = default;

    DenselySampledSpectrum(const PiecewiseLinearSpectrum &spec) 
    : values(int(Lambda_max - Lambda_min + 1)) {
        for (int lambda = Lambda_min; lambda <= Lambda_max; ++lambda)
            values[lambda - Lambda_min] = spec(lambda);
    }

    Float operator()(Float lambda) const override {
        DCHECK_GT(lambda, 0);
        int offset = std::lround(lambda) - Lambda_min;
        if (offset < 0 || offset >= values.size())
            return 0;
        return values[offset];
    }

    SampledSpectrum sample(const SampledWavelengths &lambda) const override {
        SampledSpectrum s;
        for (int i = 0; i < nSpectrumSamples; ++i)
            s[i] = operator()(lambda[i]);
        return s;
    }

    DenselySampledSpectrum(Spectrum *spec)
        : values(Lambda_max - Lambda_min + 1) {
        if (spec)
            for (int lambda = Lambda_min; lambda <= Lambda_max; ++lambda)
                values[lambda - Lambda_min] = (*spec)(lambda);
    }

    Float maxValue() const override{ return *std::max_element(values.begin(), values.end()); }

private:
    std::vector<Float> values;
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

    RGB &operator/=(Float a) {
        DCHECK(!IsNaN(a));
        DCHECK_NE(a, 0);
        r /= a;
        g /= a;
        b /= a;
        return *this;
    }
    RGB operator/(Float a) const {
        RGB ret = *this;
        return ret /= a;
    }

    Float r = 0, g = 0, b = 0;
};


class XYZ {
public:
    XYZ() = default;
    XYZ(Float X, Float Y, Float Z) : X(X), Y(Y), Z(Z) {}

    Point2f xy() const { return Point2f(X / (X + Y + Z), Y / (X + Y + Z)); }

    static XYZ fromxyY(Point2f xy, Float Y = 1) {
        if (xy.y == 0)
            return XYZ(0, 0, 0);
        return XYZ(xy.x * Y / xy.y, Y, (1 - xy.x - xy.y) * Y / xy.y);
    }

    Float &operator[](int c) {
        DCHECK(c >= 0 && c < 3);
        if (c == 0)
            return X;
        else if (c == 1)
            return Y;
        return Z;
    }

    Float operator[](int c) const {
        DCHECK(c >= 0 && c < 3);
        if (c == 0)
            return X;
        else if (c == 1)
            return Y;
        return Z;
    }

    XYZ &operator/=(Float a) {
        DCHECK(!IsNaN(a));
        DCHECK_NE(a, 0);
        X /= a;
        Y /= a;
        Z /= a;
        return *this;
    }

    XYZ operator/(Float a) const {
        XYZ ret = *this;
        return ret /= a;
    }

    Float X = 0, Y = 0, Z = 0;
};


class RGBSigmoidPolynomial {
public:
    RGBSigmoidPolynomial() = default;
    RGBSigmoidPolynomial(Float c0, Float c1, Float c2) 
    : c0(c0), c1(c1), c2(c2) {}

    Float operator()(Float lambda) const {
        return s(c0*sqr(lambda) + c1*lambda + c2);
    }

    Float maxValue() const {
        Float result = std::max((*this)(360), (*this)(830));
        Float lambda = -c1 / (2 * c0);
        if (lambda >= 360 && lambda <= 830)
            result = std::max(result, (*this)(lambda));
        return result;
    }

private:
    static  Float s(Float x) {
        if (std::isinf(x)) return x > 0 ? 1 : 0;
        return .5f + x / (2 * std::sqrt(1+sqr(x)));
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


class RGBColorSpace {
public:
    RGBColorSpace(Point2f r, Point2f g, Point2f b, Spectrum *illuminant,
                  const RGBToSpectrumTable *rgbToSpectrumTable);


    RGBSigmoidPolynomial ToRGBCoeffs(RGB rgb) const {
        DCHECK(rgb.r >= 0 && rgb.g >= 0 && rgb.b >= 0);
        return (*rgbToSpectrumTable)(RGB(std::max<Float>(0, rgb.r), std::max<Float>(0, rgb.g),
                                         std::max<Float>(0, rgb.b)));
    }

    SampledSpectrum RGBToSpectrum(RGB rgb, const SampledWavelengths &lambda) const {
        RGBSigmoidPolynomial poly = ToRGBCoeffs(rgb);
        SampledSpectrum s;
        for (int i = 0; i < nSpectrumSamples; ++i)
            s[i] = poly(lambda[i]);
        return s;
    }

    static void Init();

    // RGBColorSpace Public Members
    Point2f r, g, b, w;
    DenselySampledSpectrum illuminant;
    SquareMatrix<3> XYZFromRGB, RGBFromXYZ;
    static const RGBColorSpace *DCI_P3;

    bool operator==(const RGBColorSpace &cs) const {
        return (r == cs.r && g == cs.g && b == cs.b && w == cs.w &&
                rgbToSpectrumTable == cs.rgbToSpectrumTable);
    }
    bool operator!=(const RGBColorSpace &cs) const {
        return (r != cs.r || g != cs.g || b != cs.b || w != cs.w ||
                rgbToSpectrumTable != cs.rgbToSpectrumTable);
    }

    std::string ToString() const;

    RGB LuminanceVector() const {
        return RGB(XYZFromRGB[1][0], XYZFromRGB[1][1], XYZFromRGB[1][2]);
    }

    RGB ToRGB(XYZ xyz) const { return mul<RGB>(RGBFromXYZ, xyz); }
    XYZ ToXYZ(RGB rgb) const { return mul<XYZ>(XYZFromRGB, rgb); }

    static const RGBColorSpace *GetNamed(std::string name);
    static const RGBColorSpace *Lookup(Point2f r, Point2f g, Point2f b, Point2f w);

  private:
    const RGBToSpectrumTable *rgbToSpectrumTable;
};


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

class RGBIlluminantSpectrum : public Spectrum {
  public:
    RGBIlluminantSpectrum() = default;

    RGBIlluminantSpectrum(const RGBColorSpace &cs, RGB rgb)
    : illuminant(&cs.illuminant) {
        Float m = std::max({rgb.r, rgb.g, rgb.b});
        scale = 2 * m;
        rsp = cs.ToRGBCoeffs(scale ? rgb / scale : RGB(0, 0, 0));
    }

    Float operator()(Float lambda) const override {
        if (!illuminant)
            return 0;
        return scale * rsp(lambda) * (*illuminant)(lambda);
    }

    Float maxValue() const override {
        if (!illuminant)
            return 0;
        return scale * rsp.maxValue() * illuminant->maxValue();
    }

    const DenselySampledSpectrum *Illuminant() const { return illuminant; }

    SampledSpectrum sample(const SampledWavelengths &lambda) const override {
        if (!illuminant)
            return SampledSpectrum(0);
        SampledSpectrum s;
        for (int i = 0; i < nSpectrumSamples; ++i)
            s[i] = scale * rsp(lambda[i]);
        return s * illuminant->sample(lambda);
    }

  private:
    Float scale;
    RGBSigmoidPolynomial rsp;
    const DenselySampledSpectrum *illuminant;
};


class RGBAlbedoSpectrum : public Spectrum {
public:
    Float operator()(Float lambda) const { return rsp(lambda); }
    Float MaxValue() const { return rsp.maxValue(); }

    RGBAlbedoSpectrum(const RGBColorSpace &cs, RGB rgb) {
        DCHECK_LE(std::max({rgb.r, rgb.g, rgb.b}), 1);
        DCHECK_GE(std::min({rgb.r, rgb.g, rgb.b}), 0);
        rsp = cs.ToRGBCoeffs(rgb);
    }

    SampledSpectrum Sample(const SampledWavelengths &lambda) const {
        SampledSpectrum s;
        for (int i = 0; i < nSpectrumSamples; ++i)
            s[i] = rsp(lambda[i]);
        return s;
    }

    std::string ToString() const;

  private:
    RGBSigmoidPolynomial rsp;
};

XYZ sampledSpectrumToXYZ(const SampledSpectrum& s, const SampledWavelengths& lambda);

extern PiecewiseLinearSpectrum *illumd65;

namespace Spectra {

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
