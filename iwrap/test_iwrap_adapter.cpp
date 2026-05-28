#include "io/IWrapActor.H"

#include <iostream>
#include <string>

int main()
{
    int status = 0;
    std::string message;

    reflectomiter::io::iwrap::init_code(status, message);
    std::cout << "init: " << status << " | " << message << "\n";

    reflectomiter::io::iwrap::code_step(status, message);
    std::cout << "step: " << status << " | " << message << "\n";

    std::string state;
    reflectomiter::io::iwrap::get_code_state(state, status, message);
    std::cout << "state: " << state << "\n";

    reflectomiter::io::iwrap::clean_up(status, message);
    std::cout << "cleanup: " << status << " | " << message << "\n";

    return status;
}
