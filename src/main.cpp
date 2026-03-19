#include "render/render.hpp"
#include "util/color.hpp"
#include "util/log.hpp"
#include "util/profiler.hpp"
#include "util/color.hpp"
#include "worlds/worlds.hpp"
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

    render(manyBalls());
    //SampledWavelengths lambda = SampledWavelengths::sampleUniform(0.5f);
    
    Float wavelength = 450;
    LOG_DEBUG("wavelgnth={}", wavelength);
    LOG_DEBUG("Spectral value = {}", (*illumd65)(wavelength));
    LOG_DEBUG("X value = {}", Spectra::X()(wavelength));
    LOG_DEBUG("Y value = {}", Spectra::Y()(wavelength));
    LOG_DEBUG("Z value = {}", Spectra::Z()(wavelength));



    profiler.print(true);
    sample_end(s.release());
    LOG_VERBOSE("Finished render succesfully, shutting down logging\n\n******************************************************\n\n");

    return 0;
}
