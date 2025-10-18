#include <pybind11/pybind11.h>

#include "logging.hpp"

namespace py = pybind11;

// clang-format off
namespace lahuta::bindings {
void bind_logger(py::module &m) {
  py::class_<Logger> logger(m, "Logger");

  py::enum_<Logger::LogLevel>(logger, "LogLevel")
    .value("Trace",    Logger::LogLevel::Trace)
    .value("Debug",    Logger::LogLevel::Debug)
    .value("Info",     Logger::LogLevel::Info)
    .value("Warn",     Logger::LogLevel::Warn)
    .value("Error",    Logger::LogLevel::Error)
    .value("Critical", Logger::LogLevel::Critical)
    .value("Off",      Logger::LogLevel::Off);

  py::enum_<Logger::FormatStyle>(logger, "FormatStyle")
    .value("Simple",   Logger::FormatStyle::Simple)
    .value("Detailed", Logger::FormatStyle::Detailed);

  logger
    .def_static("get_instance", &Logger::get_instance, py::return_value_policy::reference)
    .def("set_format",    &Logger::set_format)
    .def("set_log_level", &Logger::set_log_level)
    .def("log", [](Logger &logger, Logger::LogLevel level, const std::string &msg) {
      return logger.log(level, msg);
    });
}
} // namespace lahuta::bindings
