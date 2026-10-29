// A command line probe over the primitives of hicDetectStripes, for the
// comparison scripts in cpp/scripts/stripenn_oracle: each subcommand reads
// .npy inputs, runs one primitive of stripes_numpy / stripes_impl and writes
// .npy outputs, so that Python compares them bit for bit with numpy, scipy,
// scikit-image and OpenCV on the same data. Not installed.

#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

#include "../tools/stripes_npy.hpp"
#include "../tools/stripes_numpy.hpp"

namespace np = hicx::stripes::np;
namespace npy = hicx::stripes::npy;

namespace {

np::ImageF load_image(const std::string& path) {
    const npy::Array a = npy::load(path);
    np::ImageF img(static_cast<int>(a.shape[0]), static_cast<int>(a.shape[1]));
    img.v = npy::as<float>(a);
    return img;
}

void save_image(const std::string& path, const np::ImageF& img) {
    npy::save(path, img.v, {static_cast<std::size_t>(img.rows), static_cast<std::size_t>(img.cols)});
}

}  // namespace

int main(int argc, char** argv) {
    if (argc < 2) {
        std::fprintf(stderr, "usage: stripes_probe <command> ...\n");
        return 2;
    }
    const std::string cmd = argv[1];
    try {
        if (cmd == "kernel") {  // SIGMA OUT
            const double sigma = std::atof(argv[2]);
            const int radius = static_cast<int>(4.0 * sigma + 0.5);
            const auto k = np::gaussian_kernel1d(sigma, radius);
            npy::save(argv[3], k, {k.size()});
        } else if (cmd == "gauss") {  // IN SIGMA OUT
            save_image(argv[4], np::gaussian_filter_constant(load_image(argv[2]), std::atof(argv[3])));
        } else if (cmd == "canny") {  // IN SIGMA OUTDIR
            const np::ImageF img = load_image(argv[2]);
            const np::CannyStages st = np::canny_stages(img, std::atof(argv[3]));
            const std::string dir = argv[4];
            save_image(dir + "/smoothed.npy", st.smoothed);
            save_image(dir + "/isobel.npy", st.isobel);
            save_image(dir + "/jsobel.npy", st.jsobel);
            save_image(dir + "/magnitude.npy", st.magnitude);
            save_image(dir + "/suppressed.npy", st.suppressed);
            npy::save(dir + "/edges.npy", st.edges,
                      {static_cast<std::size_t>(img.rows), static_cast<std::size_t>(img.cols)});
        } else if (cmd == "box") {  // IN(f8 2D) K OUT
            const npy::Array a = npy::load(argv[2]);
            const auto plane = npy::as<double>(a);
            const auto out = np::box_filter(plane, static_cast<int>(a.shape[0]),
                                            static_cast<int>(a.shape[1]), std::atoi(argv[3]));
            npy::save(argv[4], out, a.shape);
        } else if (cmd == "gray") {  // IN(f4 n x 3) OUT
            const npy::Array a = npy::load(argv[2]);
            const auto px = npy::as<float>(a);
            std::vector<float> out(a.shape[0]);
            for (std::size_t i = 0; i < out.size(); ++i) {
                out[i] = np::rgb_to_gray(px[3 * i], px[3 * i + 1], px[3 * i + 2]);
            }
            npy::save(argv[3], out, {out.size()});
        } else if (cmd == "blocksums") {  // BIG(f8 2D) QUERIES(i8 n x 4: r0 c0 R C) OUT(f8 n x 2: sum mean)
            const npy::Array big = npy::load(argv[2]);
            const auto data = npy::as<double>(big);
            const npy::Array qa = npy::load(argv[3]);
            const auto q = npy::as<std::int64_t>(qa);
            std::vector<double> out(qa.shape[0] * 2);
            for (std::size_t i = 0; i < qa.shape[0]; ++i) {
                const double* base = data.data() + static_cast<std::size_t>(q[4 * i]) * big.shape[1] +
                                     static_cast<std::size_t>(q[4 * i + 1]);
                out[2 * i] = np::sum_block(base, static_cast<std::size_t>(q[4 * i + 2]),
                                           static_cast<std::size_t>(q[4 * i + 3]), big.shape[1]);
                out[2 * i + 1] = np::mean_block(base, static_cast<std::size_t>(q[4 * i + 2]),
                                                static_cast<std::size_t>(q[4 * i + 3]), big.shape[1]);
            }
            npy::save(argv[4], out, {qa.shape[0], 2});
        } else if (cmd == "rowmeans") {  // BIG QUERIES OUT : np.mean(slice, axis=1), concatenated
            const npy::Array big = npy::load(argv[2]);
            const auto data = npy::as<double>(big);
            const npy::Array qa = npy::load(argv[3]);
            const auto q = npy::as<std::int64_t>(qa);
            std::vector<double> out;
            for (std::size_t i = 0; i < qa.shape[0]; ++i) {
                for (std::int64_t r = 0; r < q[4 * i + 2]; ++r) {
                    const double* row = data.data() + static_cast<std::size_t>(q[4 * i] + r) * big.shape[1] +
                                        static_cast<std::size_t>(q[4 * i + 1]);
                    out.push_back(np::mean_of_row(row, static_cast<std::size_t>(q[4 * i + 3])));
                }
            }
            npy::save(argv[4], out, {out.size()});
        } else if (cmd == "random") {  // SEED N OUT : N random() values
            np::PyRandom rng(std::strtoull(argv[2], nullptr, 10));
            const std::size_t n = std::strtoull(argv[3], nullptr, 10);
            std::vector<double> out(n);
            for (double& v : out) {
                v = rng.random();
            }
            npy::save(argv[4], out, {n});
        } else if (cmd == "median") {  // IN OUT
            const npy::Array a = npy::load(argv[2]);
            const double m = np::median(npy::as<double>(a));
            npy::save(argv[3], &m, {1});
        } else {
            std::fprintf(stderr, "unknown command %s\n", cmd.c_str());
            return 2;
        }
    } catch (const std::exception& e) {
        std::fprintf(stderr, "error: %s\n", e.what());
        return 1;
    }
    return 0;
}
