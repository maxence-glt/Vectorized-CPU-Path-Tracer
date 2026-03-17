#include "util/color.hpp"
#include "util/log.hpp"
#include "util/profiler.hpp"
#include "util/color.hpp"
#include <OpenImageIO/imageio.h>
#include <unistd.h>

SampledSpectrum SampleRGBSpectrum(const RGBColorSpace *cs,
                                  RGB rgb,
                                  const SampledWavelengths &swl) {
    RGBSigmoidPolynomial poly = cs->ToRGBCoeffs(rgb);

    SampledSpectrum s;
    for (int i = 0; i < nSpectrumSamples; ++i)
        s[i] = poly(swl[i]);

    return s;
}

Float srgbToLinear(Float c) {
    if (c <= 0.04045f) return c / 12.92f;
    return std::pow((c + 0.055f) / 1.055f, 2.4f);
}

int main() {
    LOG_VERBOSE("Starting raytracing:");
    init();
    spectrumInit();
    auto s = sample_start("main");

    //render(manyBalls());
    //SampledWavelengths lambda = SampledWavelengths::sampleUniform(0.5f);
    
    Float wavelength = 600;
    LOG_DEBUG("wavelgnth={}", wavelength);
    LOG_DEBUG("Spectral value = {}", (*illumd65)(wavelength));
    LOG_DEBUG("X value = {}", Spectra::X()(wavelength));
    LOG_DEBUG("Y value = {}", Spectra::Y()(wavelength));
    LOG_DEBUG("Z value = {}", Spectra::Z()(wavelength));

    RGB rgb(0.588f, 0.509f, 0.784f);
    RGBSigmoidPolynomial poly = RGBColorSpace::DCI_P3->ToRGBCoeffs(rgb);
    Float sDistribution = poly(wavelength);
    LOG_DEBUG("{}", sDistribution);

    SampledWavelengths swl = SampledWavelengths::sampleUniform(0.37f);
    SampledSpectrum spec = SampleRGBSpectrum(RGBColorSpace::DCI_P3, rgb, swl);
    LOG_DEBUG("swl[0]={} spec[0]={}", swl[0], spec[0]);

    XYZ xyz = RGBColorSpace::DCI_P3->ToXYZ(rgb);
    LOG_DEBUG("X value from RGB = {}", xyz.X);
    LOG_DEBUG("Y value from RGB = {}", xyz.Y);
    LOG_DEBUG("Z value from RGB = {}", xyz.Z);

    LOG_DEBUG("X value from RGB (gamma corrected) = {}", srgbToLinear(xyz.X));
    LOG_DEBUG("Y value from RGB (gamma corercted) = {}", srgbToLinear(xyz.Y));
    LOG_DEBUG("Z value from RGB (gamma corrected) = {}", srgbToLinear(xyz.Z));

    profiler.print(true);
    sample_end(s.release());
    LOG_VERBOSE("Finished render succesfully, shutting down logging\n\n******************************************************\n\n");

    return 0;
}
