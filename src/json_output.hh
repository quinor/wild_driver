#pragma once

#include "external/json.hpp"
#include <iostream>
#include <cstdlib>

using json = nlohmann::ordered_json;


void json_fail(const std::string& text, const json& detail);
void json_log_raw(const json& message);
void json_output(const std::string& key, const json& message);

template<typename... Types>
void json_log(Types... args)
{
    json_log_raw(json::array({args...}));
}
