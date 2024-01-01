#pragma once
#include "external/linalg.h"
#include "json_output.hh"

using float2 = linalg::vec<float, 2>;
using float3 = linalg::vec<float, 3>;
using mat2 = linalg::mat<float, 2, 2>;
using mat3 = linalg::mat<float, 3, 3>;


namespace linalg
{

static inline void to_json(json& j, const float2& v)
{
    j = json{
        {"x", v.x},
        {"y", v.y}
    };
}

}
