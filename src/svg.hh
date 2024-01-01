#pragma once

#include <vector>
#include <unordered_set>
#include <unordered_map>
#include <string>

#include "external/nanosvg.h"
#include "math.hh"


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

std::pair<uint32_t, std::string> best_color_match(uint32_t color);


// parser/rasterizer implementation

struct Shape
{
    NSVGshape *shape;

    // each element is a polyline
    std::vector<std::vector<float2>> outline;
    uint32_t outline_color;

    // each element is a single hatchline
    std::vector<std::pair<float2, float2>> hatching;
    uint32_t hatching_color;

    void hatch(float hatch_interval, bool order_paths);
};

struct Color
{
    NSVGcolor color;
    std::string hexname;
    std::string bestname;
    bool enabled;
};

static inline void to_json(json& j, const Color& v)
{
    j = json{
        {"rgb", {
            {"r", v.color.r},
            {"g", v.color.g},
            {"b", v.color.b}
        }},
        {"hex", v.hexname},
        {"name", v.bestname},
        {"enabled", v.enabled}
    };
}


struct Image
{
    NSVGimage* svg;

    // extents
    float2 beg;
    float2 end;

    std::vector<Shape> shapes;

    std::unordered_map<uint32_t, Color> colors;

    std::unordered_set<std::string> color_key_set;

    explicit Image(
        const char* fname,
        bool shrink_to_size,
        bool colors_only,
        float scale,
        bool mirror,
        float min_step_size,
        std::unordered_set<std::string> color_key_set,
        bool hatch,
        float hatch_interval
    );

    ~Image();

private:
    void add_color(uint32_t c);
};
