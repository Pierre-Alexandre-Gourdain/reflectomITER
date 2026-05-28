#ifndef REFLECTOMITER_IO_IWRAP_ACTOR_H
#define REFLECTOMITER_IO_IWRAP_ACTOR_H

#include <string>

namespace reflectomiter::io::iwrap
{

struct RunOptions
{
    std::string executable{"./reflectomITER"};
    std::string input_file{"inputs.iwrap_reflectomITER"};
    std::string launcher{""};             // e.g. "mpirun", "srun"
    std::string launcher_arguments{""};   // e.g. "-np 2"
    std::string working_directory{"."};
    std::string log_file{"reflectomITER_iwrap.log"};
    bool dry_run{false};
};

RunOptions options_from_environment();

std::string build_launch_command(const RunOptions& options);

int launch_reflectomiter(const RunOptions& options,
                         std::string& status_message);

void init_code(int& status_code, std::string& status_message);

void clean_up(int& status_code, std::string& status_message);

void get_code_state(std::string& state_out,
                    int& status_code,
                    std::string& status_message);

void restore_code_state(const std::string& state,
                        int& status_code,
                        std::string& status_message);

void get_timestamp_cpp(double& timestamp_out,
                       int& status_code,
                       std::string& status_message);

void code_step(int& status_code, std::string& status_message);

} // namespace reflectomiter::io::iwrap

// iWrap's generated C++ wrapper calls these functions in the global namespace.
// Keep the real implementation namespaced, but expose global aliases here.
using reflectomiter::io::iwrap::init_code;
using reflectomiter::io::iwrap::clean_up;
using reflectomiter::io::iwrap::code_step;
using reflectomiter::io::iwrap::get_code_state;
using reflectomiter::io::iwrap::restore_code_state;
using reflectomiter::io::iwrap::get_timestamp_cpp;

#endif // REFLECTOMITER_IO_IWRAP_ACTOR_H
