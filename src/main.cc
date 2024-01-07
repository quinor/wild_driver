#include "external/CLI11.hpp"
#include "svg.hh"
#include "math.hh"
#include "json_output.hh"

#include <cstdio>

// 50 steps per mm - as per Wild TA-10's documentation
const float STEPS_PER_MM = 50;


int main(int argc, char** argv)
{
    CLI::App app{"Wild TA-10 plotter SVG driver"};
    app.get_formatter()->column_width(40);

    std::string in_fname;
    std::string out_fname;
    std::string vis_fname;
    float scale = 1;
    float2 translate={0,0};
    std::unordered_set<std::string> color_key_set;
    bool colors_only = false;
    bool shrink_to_size = false;
    bool mirror = false;
    bool box = false;
    int dry_run = 0;
    bool cut = false;
    int lift_angle = 0;
    int speed = 0;
    float min_step_size = 0.5;
    bool no_outline = false;
    bool hatch = false;
    float hatch_density = 0;
    bool visualize = true;

    // IO options
    app.add_option("input,-i,--input", in_fname, "Input svg file.")
        ->check(CLI::ExistingFile)
        ->required();

    app.add_option("output,-o,--output", out_fname, "Output plotter commands file.")
        ->default_val("out.wild");

    app.add_option("--vis", vis_fname, "Visualization svg file.")
        ->default_val("vis.svg");

    // transforms
    app.add_option<float>("--scale", scale, "Scale of the plot.")
        ->default_val(1)
        ->check(CLI::PositiveNumber);

    app.add_option<float2, std::pair<float, float>>(
        "--translate", translate,
        "Translate after scaling (ie. to move the plot on the work surface).")
        ->default_str("0 0");

    app.add_flag("--shrink_to_size", shrink_to_size, "Calculate the bounding box based on what's in "
        "the svg and not the canvas size.");

    app.add_flag("--mirror", mirror, "Mirror the image (useful for cutting).");

    // plotting options
    app.add_option<std::unordered_set<std::string>>("--color_key", color_key_set,
        "Draw only specified colors.");

    app.add_flag("--colors_only", colors_only, "Only extract colors from the svg.");

    app.add_flag("--box", box, "Draw the box around the extents of the plot instead. "
        "Forces --dry_run.");

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

    if (box)
        dry_run++;

    dry_run = dry_run > 0;

    // caution. the svg has got x going "right" and y "down". plotters use x "right" and y "up", so
    // it has to be compensated for when rasterizing the polyline!

    Image svg = Image(
        in_fname.c_str(),
        shrink_to_size,
        colors_only,
        scale,
        mirror,
        min_step_size,
        color_key_set,
        hatch,
        hatch_density
    );

    // colors data already dumped during parsing/rasterization so we can just exit
    if (colors_only)
        return 0;

    if (visualize)
    {
        json_log("Saving visualization to file.", vis_fname);
        json_output("vis_file", vis_fname);

        FILE* vis = fopen(vis_fname.c_str(), "w");

        float2 size = svg.end-svg.beg;
        fprintf(vis, "<?xml version=\"1.0\" standalone=\"no\"?>\n");
        fprintf(vis, "<!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 1.1//EN\"\n");
        fprintf(vis, "\"http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd\">\n");
        fprintf(vis, "<svg xmlns=\"http://www.w3.org/2000/svg\" version=\"1.1\"\n");
        fprintf(vis, "viewBox=\"%f %f %f %f\">\n", svg.beg.x, svg.beg.y, size.x, size.y);
        fprintf(
            vis,
            "<rect x=\"%f\" y=\"%f\" width=\"100%%\" height=\"100%%\" fill=\"#def2ff\" "
            "stroke=\"black\" stroke-dasharray=\"1\" stroke-width=\"0.1\"/>",
            svg.beg.x,
            svg.beg.y
        );

        float2 last = svg.beg;
        for (auto& shape : svg.shapes)
        {
            if (shape.outline.size())
            {
                for (auto& path : shape.outline)
                {
                    float2 cur = path.front();
                    fprintf(vis, "<polyline points=\"");
                    fprintf(vis, "%f,%f ", last.x, last.y);
                    fprintf(vis, "%f,%f ", cur.x, cur.y);
                    fprintf(vis, "\" fill=\"none\" stroke=\"red\" stroke-dasharray=\"1\" stroke-width=\"0.1\"/>\n");

                    fprintf(vis, "<polyline points=\"");
                    for (auto p : path)
                        fprintf(vis, "%f,%f ", p.x, p.y);
                    fprintf(vis, "\" fill=\"none\" stroke=\"%s\" stroke-width=\"0.2\"/>\n",
                        svg.colors[shape.outline_color].hexname.c_str());

                    last = path.back();
                }
            }
            if (shape.hatching.size())
            {
                for (auto segment : shape.hatching)
                {
                    fprintf(vis, "<polyline points=\"");
                    fprintf(vis, "%f,%f ", last.x, last.y);
                    fprintf(vis, "%f,%f ", segment.first.x, segment.first.y);
                    fprintf(vis, "\" fill=\"none\" stroke=\"red\" stroke-dasharray=\"1\" stroke-width=\"0.1\"/>\n");

                    fprintf(vis, "<polyline points=\"");
                    fprintf(vis, "%f,%f ", segment.first.x, segment.first.y);
                    fprintf(vis, "%f,%f ", segment.second.x, segment.second.y);
                    fprintf(vis, "\" fill=\"none\" stroke=\"%s\" stroke-width=\"0.2\"/>\n",
                        svg.colors[shape.hatching_color].hexname.c_str());

                    last = segment.second;
                }
            }
        }
        fprintf(vis, "</svg>\n");
        fclose(vis);
    }

    // points to plot - the boolean flag specifies pen up (0) / down (1)
    std::vector<std::pair<float2, bool>> plt;

    // draw the box in the box mode
    if (box)
    {
        json_log("Tracing the bounding box.");
        plt.emplace_back(float2{svg.beg.x, svg.end.y}, 0);
        plt.emplace_back(svg.end, 0);
        plt.emplace_back(float2{svg.end.x, svg.beg.y}, 0);
        plt.emplace_back(svg.beg, 0);
        plt.emplace_back(float2{svg.beg.x, svg.end.y}, 0);
    }
    else
    {
        json_log("Tracing the polyline.");
        // draw the polyline otherwise
        for (auto& shape : svg.shapes)
        {
            if (shape.outline.size())
            {
                for (auto& path : shape.outline)
                {
                    bool init = 0;
                    for (auto p : path)
                    {
                        plt.emplace_back(p, init);
                        init = 1;
                    }
                }
            }
            if (shape.hatching.size())
            {
                for (auto segment : shape.hatching)
                {
                    plt.emplace_back(segment.first, 0);
                    plt.emplace_back(segment.second, 1);
                }
            }
        }
    }

    json_output("settings", {
        {"box", box},
        {"dry_run", dry_run},
        {"cut", cut}
    });

    // apply the plotter coordinate system:
    // * make svg.beg.x and svg.end.y be the origin
    // * flip the y axis
    // * apply the translate
    json_output("translate", translate);
    for (auto& p : plt)
    {
        p.first.x = p.first.x - svg.beg.x + translate.x;
        p.first.y = svg.end.y - p.first.y + translate.y;
    }

    // writing the .wild file

    json_log("Saving the output to file.", out_fname);
    json_output("out_file", out_fname);

    FILE* output = fopen(out_fname.c_str(), "w");
    if (output == nullptr)
    {
        json_fail("Could not open output file.", {{"filename", out_fname}});
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
    fprintf(output, ":7%d,37\r", speed == 0 ? (cut ? 15 : 37) : speed);

    float2 cur = {0,0};
    int lastx = 0, lasty = 0;
    float total_up = 0, total_down = 0;

    for (auto& p: plt)
    {
        float2 vec = p.first;
        bool draw = p.second && (!dry_run);

        // for stats
        float delta = linalg::length(vec-cur);
        cur = vec;
        if (draw)
            total_down += delta;
        else
            total_up += delta;

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

    // A buffer at the end is needed for the plotter not to stop moving. Weird.
    {
        for (int i=0; i<1024; i++)
            fprintf(output, "U0,0\r");
    }

    fclose(output);

    json_output("stats", {
        {"points", plt.size()},
        {"length_down", total_down},
        {"length_up", total_up},
    });

    json_log("done!");

    return 0;
}
