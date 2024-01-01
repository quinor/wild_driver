#include "svg.hh"
#include "json_output.hh"
#include "external/linalg.h"

#include <limits>
#include <random>
#include <algorithm>
#include <set>


namespace
{

const float DPI = 96;


uint32_t remove_a(uint32_t x)
{
    NSVGcolor c;
    c.color = x;
    c.a = 0;
    return c.color;
}

void rasterize_bezier (
    std::vector<float2>& pln,
    float2 a,
    float2 b,
    float2 c,
    float2 d,
    float min_step_size,
    int depth
)
{
    float2 p = (a+b)*0.5f, q=(b+c)*0.5f, r=(c+d)*0.5f;
    float2 k = (p+q)*0.5f, l = (q+r)*0.5f;
    float2 x = (k+l)*0.5f;
    if (
        fabs(a.x + c.x - b.x - b.x) +
        fabs(a.y + c.y - b.y - b.y) +
        fabs(b.x + d.x - c.x - c.x) +
        fabs(b.y + d.y - c.y - c.y)
        < min_step_size * 0.25
        || depth > 10
    )
    {
        // assumes pln is not empty which is supposed to be true - the beginning of a polyline is already in
        auto prev = pln.back();
        auto next = x;

        if (linalg::length2(next-prev) > min_step_size * min_step_size)
            pln.emplace_back(x);
        return;
    }
    rasterize_bezier(pln, a, p, k, x, min_step_size, depth+1);
    rasterize_bezier(pln, x, l, r, d, min_step_size, depth+1);
};


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


}


void Shape::hatch (
    float hatch_interval,
    bool order_paths
)
{
    std::vector<Section> sections1, sections2;

    // hatch at a (random) angle.
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

    // generate all sections and sort them
    for (auto& cur : outline)
    {
        for (size_t j=1; j<cur.size(); j++)
        {
            float2 p1 = linalg::mul(rot, cur[j-1]);
            float2 p2 = linalg::mul(rot, cur[j]);
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
        return;

    // min/max y range
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

    // sort sections by when they spawn on the Y axis...
    std::sort(sections1.begin(), sections1.end(), [](
        const Section& a, const Section& b) -> bool
        { return a.first.y < b.first.y; }
    );

    // and when they despawn
    std::sort(sections2.begin(), sections2.end(), [](
        const Section& a, const Section& b) -> bool
        { return a.second.y < b.second.y; }
    );

    // all the currently existing sections
    std::set<Section> active;

    auto it1 = sections1.begin();
    auto it2 = sections2.begin();

    // attempt a hatchline from ymin to ymax every hatch_interval
    for (float y=ymin; y<ymax; y+=hatch_interval)
    {
        // bookeeping of sections1 and sections2
        while (it1 != sections1.end() && it1->first.y <= y)
            active.insert(*(it1++));

        while (it2 != sections2.end() && it2->second.y <= y)
            active.erase(*(it2++));

        // generate all intersection points with the current Y level
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

        // create the hatchlines
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
                    hatching.emplace_back(
                        linalg::mul(unrot, {beg, y}),
                        linalg::mul(unrot, {p.first, y})
                    );
                }
            }
    }

    // no hatchlines created
    if (hatching.size() == 0)
        return;

    if (!order_paths)
        return;

    // ordering by nearest neighbor
    for (size_t i=1; i<hatching.size(); i++)
    {
        float2 last = hatching[i-1].second;
        for (size_t j=i+1; j<hatching.size(); j++)
            if (
                linalg::length2(last-hatching[i].first) >
                linalg::length2(last-hatching[j].first)
            )
                std::swap(hatching[i], hatching[j]);
    }
}


void Image::add_color(uint32_t c)
{
    if (colors.find(c) != colors.end())
    {
        return;
    }

    Color cur;
    cur.color.color = c;

    char hexname[8];
    sprintf(hexname, "#%02x%02x%02x", cur.color.r, cur.color.g, cur.color.b);
    cur.hexname = hexname;

    auto best = best_color_match(c);
    cur.bestname = best.second;

    cur.enabled = true;

    if (
        color_key_set.size() &&
        color_key_set.count(cur.hexname) == 0 &&
        color_key_set.count(cur.bestname) == 0
    )
        cur.enabled = false;

    colors[c] = cur;
}


Image::Image(
    const char* fname,
    bool shrink_to_size,
    bool colors_only,
    float scale,
    bool mirror,
    float min_step_size,
    std::unordered_set<std::string> m_color_key_set,
    bool hatch,
    float hatch_interval
)
: color_key_set(m_color_key_set)
{
    // rasterization is done before scaling
    min_step_size /= scale;
    hatch_interval /= scale;

    bool order_paths = true;

    json_log("Parsing input file", fname);

    svg = nsvgParseFromFile(fname, "mm", DPI);

    beg = {0, 0};
    end = {svg->width, svg->height};
    end *= 25.4f / DPI;

    json_log("Parsed svg of size (mm)", end);
    json_output("input_size", {{"begin", beg}, {"end", end}, {"size", end-beg}});

    // add the default color - black
    add_color(0);

    for (NSVGshape *shape = svg->shapes; shape != NULL; shape = shape->next)
    {
        if (!(shape->flags & NSVG_FLAGS_VISIBLE))
            continue;

        shapes.emplace_back();
        Shape& cur = shapes.back();
        cur.shape = shape;

        // the default (fallback) color is black
        cur.outline_color = cur.hatching_color = 0;

        if (shape->stroke.type == NSVG_PAINT_COLOR)
        {
            cur.outline_color = remove_a(shape->stroke.color);
            add_color(cur.outline_color);
        }

        if (shape->fill.type == NSVG_PAINT_COLOR)
        {
            cur.hatching_color = remove_a(shape->fill.color);
            add_color(cur.hatching_color);
        }

        // fall back to fill color if outline is not colored
        cur.outline_color =
            shape->stroke.type == NSVG_PAINT_COLOR
            ? cur.outline_color
            : cur.hatching_color;

        if (colors_only)
            continue;

        if (colors[cur.outline_color].enabled)
        {
            for (NSVGpath *path = shape->paths; path != NULL; path = path->next)
            {
                cur.outline.emplace_back();
                std::vector<float2>& pln = cur.outline.back();

                for (int n = 0; n < path->npts-1; n += 3)
                {
                    float* pt = &path->pts[n*2];
                    float2 a = {pt[0], pt[1]};
                    float2 b = {pt[2], pt[3]};
                    float2 c = {pt[4], pt[5]};
                    float2 d = {pt[6], pt[7]};

                    if (n == 0)
                        pln.emplace_back(a);
                    rasterize_bezier(pln, a, b, c, d, min_step_size, 0);
                    pln.emplace_back(d);
                }
            }

            if (!order_paths)
                continue;

            // ordering by nearest neighbor
            for (size_t i=1; i<cur.outline.size(); i++)
            {
                float2 last = cur.outline[i-1].back();
                for (size_t j=i+1; j<cur.outline.size(); j++)
                    if (
                        linalg::length2(last-cur.outline[i].front()) >
                        linalg::length2(last-cur.outline[j].front())
                    )
                        std::swap(cur.outline[i], cur.outline[j]);
            }
        }

        if (colors[cur.hatching_color].enabled && shape->fill.type == NSVG_PAINT_COLOR && hatch)
        {
            cur.hatch(hatch_interval, order_paths);
        }
    }

    json_log("Parsed n shapes", shapes.size());
    json_log("Parsed n colors", colors.size());

    {
        json color_list = json::array();
        for (auto& e : colors)
            color_list.push_back(e.second);
        json_output("colors", color_list);
    }

    if (colors_only)
        return;

    beg = beg*scale;
    end = end*scale;

    auto transform = [beg=beg, end=end, mirror, scale](float2 val) -> float2
    {
        float2 ret = val * scale;
        if (mirror)
            ret.x = beg.x + end.x - ret.x;
        return ret;
    };

    bool any_exist = false;

    constexpr float inf = std::numeric_limits<double>::infinity();
    float2 low={inf, inf}, hi={-inf, -inf};
    for (auto& shape : shapes)
    {
        for (auto& pts : shape.outline)
            for (float2& v : pts)
            {
                any_exist = true;
                v = transform(v);
                low = linalg::min(low, v);
                hi = linalg::max(hi, v);
            }
        for (auto& p: shape.hatching)
        {
            p.first = transform(p.first);
            p.second = transform(p.second);
        }
    }

    if (!any_exist)
    {
        json_fail("No paths rendered. Do you have non-path objects (ie. text)?", {});
    }

    if (shrink_to_size)
    {
        beg = low;
        end = hi;
    }
    if (
        low.x < beg.x || low.y < beg.y ||
        hi.x > end.x || hi.y > end.y
    )
    {
        json_fail("Points outside bounding box. Maybe use `--shrink_to_size`?", {
            {"bounding_box", {
                {"begin", beg},
                {"end", end}
            }},
            {"min_max_coords", {
                {"begin", beg},
                {"end", end}
            }},
        });
    }

    json_log("Transformed svg extents", beg, end);
    json_output("transformed_size", {{"begin", beg}, {"end", end}, {"size", end-beg}});
}

Image::~Image()
{
    nsvgDelete(svg);
}
