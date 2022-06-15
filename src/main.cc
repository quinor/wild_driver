#include "external/CLI11.hpp"
#include "external/linalg.h"
#define NANOSVG_IMPLEMENTATION
#define NANOSVG_ALL_COLOR_KEYWORDS
#include "external/nanosvg.h"

#include <cstdio>
#include <cmath>
#include <random>
#include <unordered_set>
#include <unordered_map>

using float2 = linalg::vec<float, 2>;
using mat2 = linalg::mat<float, 2, 2>;



// my addition to NSVG

typedef union NSVGcolor {
	unsigned int color;
	struct {
		unsigned char r;
		unsigned char g;
		unsigned char b;
		unsigned char a;
	};
} NSVGcolor;


// 50 steps per mm - as per Wild TA-10's documentation
const float STEPS_PER_MM = 50;


int color_distance (uint32_t c1p, uint32_t c2p)
{
    NSVGcolor c1 = {.color=c1p};
    NSVGcolor c2 = {.color=c2p};

    int dr = c1.r - c2.r;
    int dg = c1.g - c2.g;
    int db = c1.b - c2.b;
    return dr*dr + dg*dg + db*db;
}

struct Section {
    float2 first;
    float2 second;
    int dir;
};

bool operator< (const Section& a, const Section& b)
{
    if (a.first.x != b.first.x)
        return a.first.x < b.first.x;

    if (a.second.x != b.second.x)
        return a.second.x < b.second.x;

    if (a.first.y != b.first.y)
        return a.first.y < b.first.y;

    if (a.second.y != b.second.y)
        return a.second.y < b.second.y;

    return a.dir < b.dir;
}

struct Path {
    bool dash = false;
    std::vector<float2> pts;
};



bool hatch_impl (
    std::vector<Path>& paths,
    size_t idx_beg,
    size_t idx_end,
    NSVGshape* shape,
    float hatch_density,
    bool order_paths,
    Path& out
)
{
    std::vector<Section> sections1, sections2;

    static std::mt19937 gen(0);
    std::uniform_real_distribution<> rd(0.f, 1.f);
    float angle = std::acos(-1) * rd(gen); // pi * rand
    mat2 rot = {
        {cosf(angle), -sinf(angle)},
        {sinf(angle), cosf(angle)},
    };
    mat2 unrot = {
        {cosf(-angle), -sinf(-angle)},
        {sinf(-angle), cosf(-angle)},
    };


    for (size_t i=idx_beg; i < idx_end; i++)
    {
        Path& cur = paths[i];
        for (size_t j=1; j<cur.pts.size(); j++)
        {
            float2 p1 = linalg::mul(rot, cur.pts[j-1]);
            float2 p2 = linalg::mul(rot, cur.pts[j]);
            if (p1.y == p2.y)
                continue;
            int dir = 1;
            if (p1.y > p2.y)
            {
                std::swap(p1, p2);
                dir = -1;
            }

            sections1.push_back({p1, p2, dir});
            sections2.push_back({p1, p2, dir});
        }
    }

    if (sections1.size() == 0)
        return false;

    float ymin, ymax;
    ymin = sections1[0].first.y;
    ymax = sections1[0].first.y;
    for (auto& section: sections1)
    {
        ymin = std::min(ymin, section.first.y);
        ymin = std::min(ymin, section.second.y);
        ymax = std::max(ymax, section.first.y);
        ymax = std::max(ymax, section.second.y);
    }


    std::sort(sections1.begin(), sections1.end(), [](
        const Section& a, const Section& b) -> bool
        { return a.first.y < b.first.y; }
    );

    std::sort(sections2.begin(), sections2.end(), [](
        const Section& a, const Section& b) -> bool
        { return a.second.y < b.second.y; }
    );

    std::set<Section> active;

    auto it1 = sections1.begin();
    auto it2 = sections2.begin();

    out.dash = true;

    std::vector<float2> pts;

    for (float y=ymin; y<ymax; y+=hatch_density)
    {
        while (it1 != sections1.end() && it1->first.y <= y)
            active.insert(*(it1++));

        while (it2 != sections2.end() && it2->second.y <= y)
            active.erase(*(it2++));

        std::vector<std::pair<float, int>> edges;
        for (auto& s : active)
            edges.emplace_back(
                (y - s.first.y) / (s.second.y - s.first.y) * (s.second.x - s.first.x) + s.first.x,
                s.dir
            );
        std::sort(edges.begin(), edges.end());
        std::reverse(edges.begin(), edges.end());
        int level = 0;
        float beg;
        for (auto& p : edges)
            if (
                shape->fillRule == NSVG_FILLRULE_EVENODD
                ? (level&1) == 0
                : level == 0
            )
            {
                beg = p.first;
                level += p.second;
            }
            else
            {
                level += p.second;
                if (
                    shape->fillRule == NSVG_FILLRULE_EVENODD
                    ? (level&1) == 0
                    : level == 0
                )
                {
                    pts.push_back(linalg::mul(unrot, {beg, y}));
                    pts.push_back(linalg::mul(unrot, {p.first, y}));
                }
            }
    }

    if (pts.size() == 0)
        return false;

    if (!order_paths)
    {
        out.pts = pts;
        return true;
    }

    std::vector<bool> visited(pts.size()/2, false);
    visited[0] = true;
    out.pts.push_back(pts[0]);
    out.pts.push_back(pts[1]);

    float2 last = out.pts.back();
    for (size_t i=1; i<pts.size()/2; i++)
    {
        int best = -1;
        float best_len = 0;
        for (size_t j=0; j<pts.size()/2; j++)
        {
            if (visited[j])
                continue;
            float new_len = linalg::length2(last-pts[2*j]);
            if (best == -1 || new_len < best_len)
            {
                best = j;
                best_len = new_len;
            }
        }
        last = pts[2*best+1];
        out.pts.push_back(pts[2*best]);
        out.pts.push_back(pts[2*best+1]);
        visited[best] = true;
    }

    return true;
}


int main(int argc, char** argv)
{
    CLI::App app{"Wild TA-10 plotter SVG driver"};
    app.get_formatter()->column_width(40);

    std::string in_fname;
    int dpi = 96;
    std::string out_fname;
    std::string vis_fname;
    float scale = 1;
    float2 pre_translate={0,0};
    float2 translate={0,0};
    std::unordered_set<std::string> color_key_set;
    bool colors_only = false;
    bool box = false;
    int dry_run = 0;
    bool cut = false;
    int lift_angle = 0;
    int speed = 0;
    float min_step_size = 0.5;
    bool order_paths = true;
    bool no_outline = false;
    bool hatch = false;
    float hatch_density = 0;
    bool visualize = true;
    bool json = false;

    // IO options
    app.add_option("input,-i,--input", in_fname, "Input svg file.")
        ->check(CLI::ExistingFile)
        ->required();

    app.add_option("output,-o,--output", out_fname, "Output plotter commands file.")
        ->default_val("out.wild");

    app.add_option("--vis", vis_fname, "Visualization svg file.")
        ->default_val("vis.svg");

    app.add_flag("--json", json, "Instead of regular stdout, output machine-readable json");

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
    app.add_option<std::unordered_set<std::string>>("--color_key", color_key_set,
        "Draw only specified colors.");

    app.add_flag("--colors_only", colors_only, "Only extract colors from the svg.");

    app.add_flag("--box", box, "Draw the box around the extents of the plot instead. "
        "Forces --dry_run unless otherwise specified.");

    app.add_flag("--dry_run,!--no_dry_run", dry_run, "Dry run (have the tool lifted "
        "when drawing). Off by default.");

    app.add_flag("--cut", cut, "Select cutting mode (initialize the cutting head).");

    app.add_option<int>("--lift_angle", lift_angle, "For cutting only. Threshold (in degrees) "
        "after which cutting head will be lifted before changing directions. 0 disables "
        "this feature.")
        ->default_val(20);

    app.add_option<int>("--speed", speed, "Plotting/cutting speed in x*8mm/s.")
        ->check(CLI::Range(1, 37));

    // rasterization options
    app.add_flag("!--no_order_paths", order_paths, "Optimize ordering of rendered paths to save on "
        "free (tool up) travel length. Warning: may slow down for over 5000 paths rendered (O(n^2) complexity).");

    app.add_option<float>("--min_step_size", min_step_size, "Minimal step size for the rasterization "
        "process in mm. Smaller values will increase resolution, but may slow down the plotting "
        "speed due to rs232 transfer rate limits.")
        ->default_val(0.5)
        ->check(CLI::PositiveNumber);

    app.add_flag("--no_outline", no_outline, "Don't trace outlines of shapes");

    app.add_flag("--hatch", hatch, "Hatch the inside of (filled in) svg shapes");

    app.add_option<float>("--hatch_density", hatch_density, "Distance (in mm) between two consecutive hatch lines.")
        ->default_val(2)
        ->check(CLI::PositiveNumber);

    // misc options
    app.add_flag("!--no_visualize", visualize, "Visualize the path as an svg file.");


    CLI11_PARSE(app, argc, argv);

    FILE* text_out = stdout;
    FILE* json_out = fopen("/dev/null", "w");

    bool color_key = color_key_set.size() != 0;

    if (json)
        std::swap(text_out, json_out);

    if (box)
        dry_run++;

    dry_run = dry_run > 0;

    auto transform = [&](float2 p) -> float2 {
        return (p+pre_translate)*scale + translate;
    };

    fprintf(text_out, "parsing input file %s ...\n", in_fname.c_str());

    fprintf(json_out, "{\n");

    // caution. the svg has got x going "right" and y "down". plotters use x "right and y "up", so
    // it has to be compensated for when rasterizing the polyline!
    NSVGimage* svg;
    svg = nsvgParseFromFile(in_fname.c_str(), "mm", dpi);
    float2 svg_size = {svg->width, svg->height};
    // convert points/dots/??? to mm
    svg_size = svg_size * (25.4f / dpi);

    float2 beg = transform({0,0}), end = transform(svg_size);

    fprintf(text_out, "input image of size %.2fmm by %.2fmm\n", svg_size.x, svg_size.y);
    fprintf(text_out, "transformed extents are:\n    x: (%.2fmm, %.2fmm)\n    y: (%.2fmm, %.2fmm)\n",
        beg.x, end.x, beg.y, end.y);

    fprintf(json_out, "    \"input_size\": { \"x\": %.2f, \"y\": %.2f },\n", svg_size.x, svg_size.y);
    fprintf(json_out,
        "    \"transformed_size\": { "
        "\"begin\": { \"x\": %.2f, \"y\": %.2f }, "
        "\"end\": { \"x\": %.2f, \"y\": %.2f } },\n",
        beg.x, beg.y, end.x, end.y
    );


    if (beg.x < 0 || beg.y < 0 || end.x < 0 || end.y < 0)
    {
        fprintf(text_out, "Error! Attempting to draw in area with negative coordinates. "
            "Please translate the plot.\n");
        fprintf(json_out, "    \"error\": \"Attempting to draw in area with negative coordinates!\"\n}\n");
        return -1;
    }

    std::unordered_set<uint32_t> all_colors;
    std::unordered_map<uint32_t, NSVGNamedColor> color_names;
    std::unordered_map<uint32_t, bool> color_enabled;

    // gather all colors
    for (NSVGshape *shape = svg->shapes; shape != NULL; shape = shape->next)
    {
        if (!(shape->flags & NSVG_FLAGS_VISIBLE))
            continue;
        if (shape->stroke.type == NSVG_PAINT_COLOR)
            all_colors.insert(shape->stroke.color);
        if (shape->fill.type == NSVG_PAINT_COLOR)
            all_colors.insert(shape->fill.color);
    }

    fprintf(text_out, "color names detected from the image:\n");
    fprintf(json_out, "    \"colors\": [");
    size_t ncolors = sizeof(nsvg__colors)/sizeof(NSVGNamedColor);
    bool first = true;
    for (auto color : all_colors)
    {
        NSVGcolor unpacked = {.color=color};

        NSVGNamedColor best = nsvg__colors[0];
        for (size_t i=0; i<ncolors; i++)
        {
            NSVGNamedColor cur = nsvg__colors[i];
            if (color_distance(color, best.color) > color_distance(color, cur.color))
                best = cur;
        }
        char hexname[8];
        sprintf(hexname, "#%02x%02x%02x", unpacked.r, unpacked.g, unpacked.b);

        color_names[color] = best;
        color_enabled[color] = (color_key_set.count(best.name) + color_key_set.count(hexname) > 0);
        fprintf(text_out, "    rgb(%3d, %3d, %3d) #%02x%02x%02x %16s %s\n",
            unpacked.r, unpacked.g, unpacked.b,
            unpacked.r, unpacked.g, unpacked.b,
            best.name, (color_enabled[color] || !color_key) ? "(enabled)" : ""
        );
        fprintf(json_out, "%s\n", first ? "" : ",");
        first = false;
        fprintf(json_out,
            "        { \"rgb\": { \"r\": %d, \"g\": %d, \"b\": %d }, "
            "\"hex\": \"#%02x%02x%02x\", \"name\": \"%s\", \"enabled\": %s }",
            unpacked.r, unpacked.g, unpacked.b,
            unpacked.r, unpacked.g, unpacked.b,
            best.name, (color_enabled[color] || !color_key) ? "true" : "false"
        );
    }

    if (colors_only)
    {
        fprintf(json_out, "\n    ]\n");
        fprintf(json_out, "}\n");
        return 0;
    }
    else
        fprintf(json_out, "\n    ],\n");

    std::vector<Path> paths;

    // draw the box in the box mode
    if (box)
    {
        fprintf(text_out, "drawing the box outline of the plot only\n");
        fprintf(json_out, "    \"box\": true,\n");
        paths.emplace_back();
        std::vector<float2>& pln = paths.back().pts;

        pln.emplace_back(beg);
        pln.emplace_back(float2{beg.x, end.y});
        pln.emplace_back(end);
        pln.emplace_back(float2{end.x, beg.y});
        pln.emplace_back(beg);
    }
    // otherwise rasterize the svg
    else
    {
        fprintf(text_out, "rasterizing plot%s...\n", hatch ? " and hatching" : "");
        fprintf(json_out, "    \"hatch\": %s,\n", hatch ? "true" : "false");

        std::function<void(std::vector<float2>&, float2, float2, float2, float2, int)> rasterize_bezier =
        [&](std::vector<float2>& pln, float2 a, float2 b, float2 c, float2 d, int n) -> void
        {
            float2 p = (a+b)*0.5f, q=(b+c)*0.5f, r=(c+d)*0.5f;
            float2 k = (p+q)*0.5f, l = (q+r)*0.5f;
            float2 x = (k+l)*0.5f;
            if (
                fabs(a.x + c.x - b.x - b.x) +
                fabs(a.y + c.y - b.y - b.y) +
                fabs(b.x + d.x - c.x - c.x) +
                fabs(b.y + d.y - c.y - c.y)
                < min_step_size * 0.5
                || n > 10
            )
            {
                pln.emplace_back(transform(x));
                return;
            }
            rasterize_bezier(pln, a, p, k, x, n+1);
            rasterize_bezier(pln, x, l, r, d, n+1);
        };

        for (NSVGshape *shape = svg->shapes; shape != NULL; shape = shape->next)
        {
            if (!(shape->flags & NSVG_FLAGS_VISIBLE))
                continue;

            size_t idx_beg = paths.size();

            unsigned int outline_color =
                shape->stroke.type == NSVG_PAINT_COLOR
                ? shape->stroke.color
                : shape->fill.color;

            if (!no_outline && (!color_key || color_enabled[outline_color]))
            {
                for (NSVGpath *path = shape->paths; path != NULL; path = path->next)
                {
                    paths.emplace_back();
                    std::vector<float2>& pln = paths.back().pts;

                    for (int n = 0; n < path->npts-1; n += 3)
                    {
                        float* pt = &path->pts[n*2];
                        float2 a = {pt[0], svg_size.y-pt[1]};
                        float2 b = {pt[2], svg_size.y-pt[3]};
                        float2 c = {pt[4], svg_size.y-pt[5]};
                        float2 d = {pt[6], svg_size.y-pt[7]};

                        if (n == 0)
                            pln.emplace_back(transform(a));
                        rasterize_bezier(pln, a, b, c, d, 0);
                        pln.emplace_back(transform(d));
                    }
                }
            }

            unsigned int fill_color = shape->fill.color;
            if (shape->fill.type != NSVG_PAINT_COLOR || !hatch || (color_key && !color_enabled[fill_color]))
                continue;

            size_t idx_end = paths.size();
            paths.emplace_back();
            if (!hatch_impl(paths, idx_beg, idx_end, shape, hatch_density, order_paths, paths.back()))
                paths.pop_back();
        }
    }

    nsvgDelete(svg);

    if (paths.size() == 0)
    {
        fprintf(text_out, "Error! Zero paths rendered!\n");
        fprintf(json_out, "    \"error\": \"Zero paths rendered!\"\n}\n");
        return -1;
    }

    if (order_paths)
    {
        fprintf(text_out, "reordering paths...\n");
        // intentionally ignore the copy from paths to new_paths and then back
        // it doesn't matter too much performance-wise. TODO: maybe fix.
        std::vector<Path> new_paths;
        std::vector<bool> visited(paths.size(), false);
        visited[0] = true;
        new_paths.push_back(paths[0]);
        float2 last = new_paths[0].pts.back();
        for (size_t i=1; i<paths.size(); i++)
        {
            int best = -1;
            float best_len = 0;
            for (size_t j=0; j<paths.size(); j++)
            {
                if (visited[j])
                    continue;
                float new_len = linalg::length2(last-paths[j].pts[0]);
                if (best == -1 || new_len < best_len)
                {
                    best = j;
                    best_len = new_len;
                }
            }
            last = paths[best].pts.back();
            new_paths.push_back(paths[best]);
            visited[best] = true;
        }
        paths = new_paths;
    }

    if (visualize)
    {
        fprintf(text_out, "saving visualization as %s\n", vis_fname.c_str());
        FILE* vis = fopen(vis_fname.c_str(), "w");

        float2 size = (end-beg);
        float yyy = end.y + beg.y;

        fprintf(vis, "<?xml version=\"1.0\" standalone=\"no\"?>\n");
        fprintf(vis, "<!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 1.1//EN\"\n");
        fprintf(vis, "\"http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd\">\n");
        fprintf(vis, "<svg width=\"%fmm\" height=\"%fmm\"\n", size.x, size.y);
        fprintf(vis, "xmlns=\"http://www.w3.org/2000/svg\" version=\"1.1\"\n");
        fprintf(vis, "viewBox=\"%f %f %f %f\">\n", beg.x, beg.y, end.x, end.y);

        for (size_t i=0; i<paths.size(); i++)
        {
            if (i != 0)
            {
                float2 p1 = paths[i-1].pts.back();
                float2 p2 = paths[i].pts.front();
                fprintf(vis, "<polyline points=\"");
                fprintf(vis, "%f,%f ", p1.x, yyy-p1.y);
                fprintf(vis, "%f,%f ", p2.x, yyy-p2.y);
                fprintf(vis, "\" fill=\"none\" stroke=\"red\" stroke-width=\"0.1\"/>\n");
            }
            Path& cur = paths[i];
            if (cur.dash)
            {
                for (size_t j=1; j<cur.pts.size(); j++)
                {
                    float2 p1 = cur.pts[j-1];
                    float2 p2 = cur.pts[j];
                    fprintf(vis, "<polyline points=\"");
                    fprintf(vis, "%f,%f ", p1.x, yyy-p1.y);
                    fprintf(vis, "%f,%f ", p2.x, yyy-p2.y);
                    fprintf(vis, "\" fill=\"none\" stroke=\"%s\" stroke-width=\"0.1\"/>\n", (j&1) ? "black" : "red");
                }
            }
            else
            {
                fprintf(vis, "<polyline points=\"");
                for (auto p : cur.pts)
                    fprintf(vis, "%f,%f ", p.x, yyy-p.y);
                fprintf(vis, "\" fill=\"none\" stroke=\"black\" stroke-width=\"0.1\"/>\n");
            }
        }

        fprintf(vis, "</svg>\n");
        fclose(vis);
    }

    // writing the .wild file

    FILE* output = fopen(out_fname.c_str(), "w");
    if (output == nullptr)
    {
        fprintf(text_out, "could not open output file: %s\n", out_fname.c_str());
        return -1;
    }

    // .wild file preamble

    // choosing tool
    fprintf(output, ":8%d\r", cut ? 1 : 3);
    // circle resolution - default is double (enhanced) resolution
    // probably not strictly needed (this driver doesn't emit circles)
    fprintf(output, ":32\r");
    // lifting/lowering time in ms
    fprintf(output, ":525,25\r");
    // "Auto pen-lift at angular discontinuities"
    fprintf(output, ":E%d00\r", cut ? (lift_angle <= 0 ? -20 : lift_angle) : -20);
    // index of tool head chosen. Hardcoded to 2.
    fprintf(output, "P2\r");
    // move to 0, 0 coordinate first
    fprintf(output, "U0,0\r");
    // speed selection. Default of 120mm/s for cutting and 300mm/s (max) for drawing.
    fprintf(output, ":7%d\r", speed == 0 ? (cut ? 15 : 37) : speed);

    fprintf(text_out, "using the %s\n", cut ? "cutter" : "pen");
    fprintf(text_out, "the tool is %s\n", dry_run ? "lifted" : "lowered");

    fprintf(json_out, "    \"cut\": %s,\n", cut ? "true" : "false");
    fprintf(json_out, "    \"dry_run\": %s,\n", dry_run ? "true" : "false");

    float2 cur = {0,0};
    int lastx = 0, lasty = 0;
    float total_up = 0, total_down = 0;
    int points_count = 0;

    for (auto& pln : paths)
    {
        bool first = true;
        for (size_t i=0; i<pln.pts.size(); i++)
        {
            float2 vec = pln.pts[i];
            bool draw = !first && !dry_run && (pln.dash <= (i&1));
            first = false;

            if (!first && linalg::length2(cur-vec) < min_step_size*min_step_size)
                continue;

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
            int dx = x-lastx;
            int dy = y-lasty;
            if (dx > 8191 || dx < -8192 || dy > 8191 || dy < -8192)
                fprintf(output, "%c%d,%d\r", draw ? 'D' : 'U', x, y);
            else
            {
                char xu = (dx>>7)&0x7f, xl = dx&0x7f, yu = (dy>>7)&0x7f, yl = dy&0x7f;
                fprintf(output, "%c%c%c%c%c\r", draw ? 'S' : 'T', xu|0x80, xl|0x80, yu|0x80, yl|0x80);
            }

            lastx = x;
            lasty = y;
        }
    }

    // A buffer at the end is needed for the plotter not to stop moving. Weird.
    {
        int x = beg.x*STEPS_PER_MM, y = beg.y*STEPS_PER_MM;
        for (int i=0; i<1024; i++)
            fprintf(output, "U%d,%d\r", x, y);
    }

    fclose(output);

    fprintf(json_out, "    \"stats\": {\n");
    fprintf(json_out, "        \"paths\": %lu,\n", paths.size());
    fprintf(json_out, "        \"points\": %d,\n", points_count);
    fprintf(json_out, "        \"length_down\": %.2f,\n", total_down);
    fprintf(json_out, "        \"length_up\": %.2f\n", total_up);
    fprintf(json_out, "    }\n");
    fprintf(json_out, "}\n");

    fprintf(text_out, "done!\n");
    fprintf(text_out, "stats:\n");
    fprintf(text_out, "    paths count: %lu\n", paths.size());
    fprintf(text_out, "    points count: %d\n", points_count);
    fprintf(text_out, "    total path length: %.2fmm\n", total_down);
    fprintf(text_out, "    total free travel length: %.2fmm\n", total_up);

    return 0;
}
