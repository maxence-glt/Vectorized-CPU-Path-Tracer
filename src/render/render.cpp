#include "render.hpp"
#include "util/log.hpp"
#include <OpenImageIO/imageio.h>
#include <cstddef>
#include <syncstream>

using namespace OIIO;

std::vector<Bounds2<int>> makeTiles(int W, int H, int TILE) {
    std::vector<Bounds2<int>> tiles;
    int tilesX = (W + TILE - 1) / TILE;
    int tilesY = (H + TILE - 1) / TILE;

    tiles.reserve(tilesX * tilesY);

    for (int ty = 0; ty < tilesY; ++ty) {
        for (int tx = 0; tx < tilesX; ++tx) {
            int x0 = tx * TILE;
            int y0 = ty * TILE;
            int x1 = std::min(x0 + TILE, W);
            int y1 = std::min(y0 + TILE, H);
            tiles.push_back({Point2<int>(x0, y0), Point2<int>(x1, y1)});
        }
    }
    return tiles;
}

void renderThread(const Scene &scene, const Bounds2<int> &t, int W, int H, int C, std::vector<float> &film) {
    auto cam = scene.camera;

    for (int y = t[0].y; y < t[1].y; ++y) {
        for (int x = t[0].x; x < t[1].x; ++x) {
            color pixel_color(0, 0, 0);
            for (int sample = 0; sample < cam.samples_per_pixel; ++sample) {
                Ray r = cam.get_ray(x, y);
                pixel_color += cam.ray_color(r, scene.world);
            }
            write_color(film, cam.pixel_samples_scale*pixel_color, (y*W+x) * C);
        }
    }
}

void render(const Scene &scene) {
    PROFILE_SCOPE("render");

    auto cam = scene.camera;

    const char* filename = "image.exr";
    const int xres = cam.image_width, yres = cam.image_height, channels = 3;

    LOG_VERBOSE("image_height      = {}", yres);
    LOG_VERBOSE("image_width       = {}", xres);
    LOG_VERBOSE("samples_per_pixel = {}", cam.samples_per_pixel);

    std::vector<float> pixels(xres * yres * channels);

    auto tiles = makeTiles(xres, yres, 32);
    std::atomic<std::size_t> nextTile{0},tilesLeft{tiles.size()};
    std::vector<std::thread> threads;
    threads.reserve(Options->nThreads);

    auto t1 = curr_time();

    for (int t = 0; t < Options->nThreads; ++t) {
        threads.emplace_back([&](){
            while (true) {
                int i = nextTile.fetch_add(1);
                if (i >= (int)tiles.size()) break;
                renderThread(scene, tiles[i], xres, yres, channels, pixels);

                int left = tilesLeft.fetch_sub(1, std::memory_order_relaxed) - 1;
                if ((left % 20) == 0) std::clog << "\rTiles left: " << left << "   " << std::flush;
            }
        });
    }

    for (auto& th : threads) th.join();

    auto t2 = curr_time();
    auto ms_int = diff_time<milliseconds>(t1, t2);
    LOG_VERBOSE("Finished ray tracing in {}ms", ms_int.count());

    std::unique_ptr<ImageOutput> out = ImageOutput::create("image.exr");
    DCHECK(out);

    ImageSpec spec(xres, yres, channels, TypeDesc::FLOAT);
    out->open(filename, spec);
    out->write_image(make_cspan(pixels));
    out->close();
}
