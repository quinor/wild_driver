#include "json_output.hh"
#include <functional>

namespace
{

json& json_out()
{
    static std::unique_ptr<json, std::function<void(json*)>> out(
        new json({
            {"success", true},
            {"log", json::array()},
            {"output", json::object()}
        }),
        [](json* j) -> void
        {
            std::cout << std::setw(2) << *j << std::endl;
        }
    );
    return *out;
}

}


void json_fail(const std::string& text, const json& detail)
{
    json_out()["success"] = false;
    json_out()["error"] = json({
        {"text", text},
        {"detail", detail}
    });
    ::exit(-1);
}

void json_log_raw(const json& message)
{
    json_out()["log"].push_back(message);
}

void json_output(const std::string& key, const json& message)
{
    if (json_out()["output"].contains(key))
        json_fail("`json_output`: the key already exists!", {
            {"key", key},
            {"message", message}
        });

    json_out()["output"][key] = message;
}
