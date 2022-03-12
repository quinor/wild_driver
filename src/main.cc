#include "external/CLI11.hpp"
#include "external/linalg.h"
#define NANOSVG_IMPLEMENTATION
#include "external/nanosvg.h"

#include <cstdio>

using float2 = linalg::vec<float, 2>;

// 50 steps per mm - as per Wild TA-10's documentation
const float STEPS_PER_MM = 50;


int main(int argc, char** argv)
{
    CLI::App app{"Wild TA-10 plotter SVG driver"};
    app.get_formatter()->column_width(40);

    std::string in_fname;
    int dpi = 96;
    std::string out_fname;
    float scale;
    float2 pre_translate={0,0};
    float2 translate={0,0};
    bool box = false;
    int dry_run = 0;
    bool cut = false;
    int lift_angle = 0;
    int speed = 0;
    float points_per_arch = 50;
    float min_step_size = 0.5;
    bool order_paths;

    // IO options
    app.add_option("input,-i,--input", in_fname, "Input svg file.")
        ->check(CLI::ExistingFile)
        ->required();

    app.add_option("output,-o,--output", out_fname, "Output plotter commands file.")
        ->default_val("out.wild");

    // transforms
    app.add_option<float>("--scale", scale, "Scale of the plot.")
        ->default_val(1)
        ->check(CLI::PositiveNumber);

    app.add_option<float2, std::pair<float, float>>("--pre_translate", pre_translate,
        "Translate before scaling (ie. to center the svg file).")
        ->default_str("0 0");

    app.add_option<float2, std::pair<float, float>>(
        "--translate", translate,
        "Translate after scaling (ie. to move the plot on the work surface).")
        ->default_str("0 0");

    // plotting options
    app.add_flag("--box", box, "Draw the box around the extents of the plot instead. "
        "Forces --dry_run unless otherwise specified.");

    app.add_flag("--dry_run,!--no_dry_run", dry_run, "Dry run (have the tool lifted "
        "when drawing). Off by default.");

    app.add_flag("--cut", dry_run, "Select cutting mode (initialize the cutting head).");

    app.add_option<int>("--lift_angle", lift_angle, "For cutting only. Threshold (in degrees) "
        "after which cutting head will be lifted before changing directions. 0 disables "
        "this feature.")
        ->default_val(20);

    app.add_option<int>("--speed", speed, "Plotting/cutting speed in x*8mm/s.")
        ->check(CLI::Range(1, 37));

    // rasterization options
    app.add_flag("--order_paths", order_paths, "Optimize ordering of rendered paths to save on "
        "free (tool up) travel length. Warning: may slow down for over 5000 paths rendered (O(n^2) complexity).");

    app.add_option<float, int>("--points_per_arch", points_per_arch, "Number of points per arch in rasterization. "
        "Most likely don't touch.")
        ->default_val(50)
        ->check(CLI::PositiveNumber);

    app.add_option<float>("--min_step_size", min_step_size, "Minimal step size for the rasterization "
        "process in mm. Smaller values will increase resolution, but may slow down the plotting "
        "speed due to rs232 transfer rate limits.")
        ->default_val(0.5)
        ->check(CLI::PositiveNumber);

    CLI11_PARSE(app, argc, argv);

    if (box)
        dry_run++;

    dry_run = dry_run > 0;

    auto transform = [&](float2 p) -> float2 {
        return (p+pre_translate)*scale + translate;
    };

    printf("parsing input file %s ...\n", in_fname.c_str());

    // caution. the svg has got x going "right" and y "down". plotters use x "right and y "up", so
    // it has to be compensated for when rasterizing the polyline!
    NSVGimage* svg;
    svg = nsvgParseFromFile(in_fname.c_str(), "mm", dpi);
    float2 svg_size = {svg->width, svg->height};
    // convert points/dots/??? to mm
    svg_size = svg_size * (25.4f / dpi);

    float2 beg = transform({0,0}), end = transform(svg_size);

    printf("input image of size %.2fmm by %.2fmm\n", svg_size.x, svg_size.y);
    printf("transformed extents are:\n    x: (%.2fmm, %.2fmm)\n    y: (%.2fmm, %.2fmm)\n",
        beg.x, end.x, beg.y, end.y);

    if (beg.x < 0 || beg.y < 0 || end.x < 0 || end.y < 0)
    {
        printf("Error! Attempting to draw in area with negative coordinates. "
            "Please translate the plot.");
        return -1;
    }

    std::vector<std::vector<float2>> paths;

    // draw the box in the box mode
    if (box)
    {
        printf("drawing the box outline of the plot only\n");
        paths.emplace_back();
        std::vector<float2>& pln = paths.back();

        pln.emplace_back(beg);
        pln.emplace_back(float2{beg.x, end.y});
        pln.emplace_back(end);
        pln.emplace_back(float2{end.x, beg.y});
        pln.emplace_back(beg);
    }
    // otherwise rasterize the svg
    else
    {
        printf("rasterizing plot...\n");
        for (NSVGshape *shape = svg->shapes; shape != NULL; shape = shape->next)
        {
            for (NSVGpath *path = shape->paths; path != NULL; path = path->next)
            {
                // move to the beginning of the path
                float2 first = {path->pts[0], svg_size.y-path->pts[1]};
                paths.emplace_back();
                std::vector<float2>& pln = paths.back();
                pln.emplace_back(transform(first));

                for (int n = 0; n < path->npts-1; n += 3)
                {
                    float* pt = &path->pts[n*2];
                    float2 a = {pt[0], svg_size.y-pt[1]};
                    float2 b = {pt[2], svg_size.y-pt[3]};
                    float2 c = {pt[4], svg_size.y-pt[5]};
                    float2 d = {pt[6], svg_size.y-pt[7]};
                    float2 x0 = a, x1 = 3*b-3*a, x2 = 3*c - 6*b + 3*a, x3 = d - 3*c + 3*b - a;

                    pln.emplace_back(transform(a));
                    float2 prev = a;
                    for (int i=1; i<points_per_arch; i++)
                    {
                        float p = i/points_per_arch;
                        float2 point = x0 + p*(x1 + p*(x2 + p*x3));
                        if (linalg::length2(point-prev) > min_step_size*min_step_size)
                        {
                            pln.emplace_back(transform(point));
                            prev=point;
                        }
                    }
                    pln.emplace_back(transform(d));
                }
            }
        }
    }

    nsvgDelete(svg);

    if (order_paths)
    {
        printf("reordering paths...\n");
        // intentionally ignore the copy from paths to new_paths and then back
        // it doesn't matter too much performance-wise. TODO: maybe fix.
        std::vector<std::vector<float2>> new_paths;
        std::vector<bool> visited(paths.size(), false);
        visited[0] = true;
        new_paths.push_back(paths[0]);
        float2 last = new_paths[0].back();
        for (size_t i=1; i<paths.size(); i++)
        {
            int best = -1;
            float best_len = 0;
            for (size_t j=0; j<paths.size(); j++)
            {
                if (visited[j])
                    continue;
                float new_len = linalg::length2(last-paths[j][0]);
                if (best == -1 || new_len < best_len)
                {
                    best = j;
                    best_len = new_len;
                }
            }
            last = paths[best].back();
            new_paths.push_back(paths[best]);
            visited[best] = true;
        }
        paths = new_paths;
    }

    // writing the .wild file

    FILE* output = fopen(out_fname.c_str(), "w");
    if (output == nullptr)
    {
        printf("could not open output file: %s\n", out_fname.c_str());
        return -1;
    }

    // .wild file preamble

    // choosing tool
    fprintf(output, ":8%d\r\n", cut ? 1 : 3);
    // circle resolution - default is double (enhanced) resolution
    // probably not strictly needed (this driver doesn't emit circles)
    fprintf(output, ":32\r\n");
    // lifting/lowering time in ms
    fprintf(output, ":525,25\r\n");
    // "Auto pen-lift at angular discontinuities"
    fprintf(output, ":E%d00\r\n", cut ? (lift_angle <= 0 ? -20 : lift_angle) : -20);
    // number of cutting head chosen. Hardcoded to 2.
    fprintf(output, "P2\r\n");
    // move to 0, 0 coordinate first
    fprintf(output, "U0,0\r\n");
    // speed selection. Default of 120mm/s for cutting and 300mm/s (max) for drawing.
    fprintf(output, ":7%d\r\n", speed == 0 ? (cut ? 15 : 37) : speed);

    printf("using the %s\n", cut ? "cutter" : "pen");
    printf("the tool is %s\n", dry_run ? "lifted" : "lowered");
    printf("number of paths rendered: %lu\n", paths.size());

    float2 cur = {0,0};
    float total_up = 0, total_down = 0;
    int points_count = 0;

    for (auto& pln : paths)
    {
        bool first = true;
        for (auto vec : pln)
        {
            bool draw = !first && !dry_run;
            first = false;

            // for stats
            float delta = linalg::length(vec-cur);
            cur = vec;
            if (draw)
                total_down += delta;
            else
                total_up += delta;
            points_count += 1;

            // actual plotter commands
            int x = vec.x*STEPS_PER_MM;
            int y = vec.y*STEPS_PER_MM;
            fprintf(output, "%c%d,%d\r\n", draw ? 'D' : 'U', x, y);
        }
    }

    // supposedly a buffer at the end is needed for the plotter not to stop moving. TODO test
    {
        int x = beg.x*STEPS_PER_MM, y = beg.y*STEPS_PER_MM;
        for (int i=0; i<64; i++)
            fprintf(output, "U%d,%d\r\n", x, y);
    }

    printf("done!\n");
    printf("stats:\n");
    printf("    points count: %d\n", points_count);
    printf("    total path length: %.2fmm\n", total_down);
    printf("    total free travel length: %.2fmm\n", total_up);

    return 0;
}
