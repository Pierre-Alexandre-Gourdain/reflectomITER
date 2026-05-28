#include "io/IWrapActor.H"

#include <cstdlib>
#include <fstream>
#include <sstream>
#include <stdexcept>

namespace reflectomiter::io::iwrap
{

namespace
{
int g_step = 0;
std::string g_state = "initialized";

std::string getenv_or_empty(const char* name)
{
    const char* value = std::getenv(name);
    return value ? std::string(value) : std::string();
}

bool getenv_is_true(const char* name)
{
    const std::string value = getenv_or_empty(name);
    return value == "1" || value == "ON" || value == "on" ||
           value == "TRUE" || value == "true" || value == "YES" ||
           value == "yes";
}

std::string shell_quote(const std::string& s)
{
    std::string out = "'";
    for (char c : s) {
        if (c == '\'') {
            out += "'\\''";
        } else {
            out += c;
        }
    }
    out += "'";
    return out;
}

} // namespace

RunOptions options_from_environment()
{
    RunOptions options;

    if (!getenv_or_empty("REFLECTOMITER_EXECUTABLE").empty()) {
        options.executable = getenv_or_empty("REFLECTOMITER_EXECUTABLE");
    }

    if (!getenv_or_empty("REFLECTOMITER_INPUT").empty()) {
        options.input_file = getenv_or_empty("REFLECTOMITER_INPUT");
    }

    if (!getenv_or_empty("REFLECTOMITER_LAUNCHER").empty()) {
        options.launcher = getenv_or_empty("REFLECTOMITER_LAUNCHER");
    }

    if (!getenv_or_empty("REFLECTOMITER_LAUNCHER_ARGS").empty()) {
        options.launcher_arguments = getenv_or_empty("REFLECTOMITER_LAUNCHER_ARGS");
    }

    if (!getenv_or_empty("REFLECTOMITER_WORKDIR").empty()) {
        options.working_directory = getenv_or_empty("REFLECTOMITER_WORKDIR");
    }

    if (!getenv_or_empty("REFLECTOMITER_IWRAP_LOG").empty()) {
        options.log_file = getenv_or_empty("REFLECTOMITER_IWRAP_LOG");
    }

    options.dry_run = getenv_is_true("REFLECTOMITER_DRY_RUN");

    return options;
}

std::string build_launch_command(const RunOptions& options)
{
    std::ostringstream cmd;

    cmd << "cd " << shell_quote(options.working_directory) << " && ";

    if (!options.launcher.empty()) {
        cmd << options.launcher << " ";
    }

    if (!options.launcher_arguments.empty()) {
        cmd << options.launcher_arguments << " ";
    }

    cmd << shell_quote(options.executable) << " "
        << shell_quote(options.input_file)
        << " > " << shell_quote(options.log_file)
        << " 2>&1";

    return cmd.str();
}

int launch_reflectomiter(const RunOptions& options,
                         std::string& status_message)
{
    const std::string command = build_launch_command(options);

    {
        std::ofstream log(options.log_file, std::ios::app);
        log << "[iWrap adapter] command: " << command << "\n";
        log << "[iWrap adapter] dry_run: " << (options.dry_run ? "true" : "false") << "\n";
    }

    if (options.dry_run) {
        status_message = "DRY RUN: " + command;
        return 0;
    }

    const int ret = std::system(command.c_str());

    if (ret == 0) {
        status_message = "reflectomITER completed successfully";
    } else {
        status_message = "reflectomITER failed with system return code " + std::to_string(ret);
    }

    return ret;
}

void init_code(int& status_code, std::string& status_message)
{
    g_step = 0;
    g_state = "initialized";

    status_code = 0;
    status_message = "reflectomITER iWrap adapter initialized";
}

void clean_up(int& status_code, std::string& status_message)
{
    g_state = "finalized";

    status_code = 0;
    status_message = "reflectomITER iWrap adapter finalized";
}

void get_code_state(std::string& state_out,
                    int& status_code,
                    std::string& status_message)
{
    state_out = g_state + ";step=" + std::to_string(g_step);

    status_code = 0;
    status_message = "state returned";
}

void restore_code_state(const std::string& state,
                        int& status_code,
                        std::string& status_message)
{
    g_state = state;

    status_code = 0;
    status_message = "state restored";
}

void get_timestamp_cpp(double& timestamp_out,
                       int& status_code,
                       std::string& status_message)
{
    timestamp_out = static_cast<double>(g_step);

    status_code = 0;
    status_message = "timestamp returned";
}

void code_step(int& status_code, std::string& status_message)
{
    const RunOptions options = options_from_environment();

    status_code = launch_reflectomiter(options, status_message);

    if (status_code == 0) {
        ++g_step;
        g_state = "completed";
    } else {
        g_state = "failed";
    }
}

} // namespace reflectomiter::io::iwrap
