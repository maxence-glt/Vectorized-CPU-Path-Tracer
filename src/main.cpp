#include "options.hpp"
#include "render/render.hpp"
#include "util/color.hpp"
#include "util/math.hpp"
#include "worlds/worlds.hpp"
#include "util/log.hpp"
#include "util/profiler.hpp"
#include "util/color.hpp"
#include <OpenImageIO/imageio.h>
#include <unistd.h>

using namespace OIIO;

int main() {
    LOG_VERBOSE("Starting raytracing:");
    init();
    auto s = sample_start("main");

    //render(manyBalls());
    //SampledWavelengths lambda = SampledWavelengths::sampleUniform(0.5f);
    

    profiler.print(true);
    sample_end(s.release());
    LOG_VERBOSE("Finished render succesfully, shutting down logging\n\n******************************************************\n\n");

    return 0;
}
