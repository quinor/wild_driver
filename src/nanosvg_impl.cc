// enable nanosvg implementation in this source file
#define NANOSVG_IMPLEMENTATION
#define NANOSVG_ALL_COLOR_KEYWORDS
#include "svg.hh"


namespace
{

int color_distance (uint32_t c1p, uint32_t c2p)
{
    NSVGcolor c1 = {.color=c1p};
    NSVGcolor c2 = {.color=c2p};

    int dr = c1.r - c2.r;
    int dg = c1.g - c2.g;
    int db = c1.b - c2.b;
    return dr*dr + dg*dg + db*db;
}

}

std::pair<uint32_t, std::string> best_color_match(uint32_t color)
{
    size_t ncolors = sizeof(nsvg__colors)/sizeof(NSVGNamedColor);

    NSVGNamedColor best = nsvg__colors[0];

    for (size_t i=0; i<ncolors; i++)
    {
        NSVGNamedColor cur = nsvg__colors[i];
        if (color_distance(color, best.color) > color_distance(color, cur.color))
            best = cur;
    }
    return {best.color, best.name};
}
