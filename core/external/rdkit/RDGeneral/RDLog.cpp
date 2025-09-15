//
// Copyright (C)  2005-2022 Greg Landrum and other RDKit contributors
//
//  @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//
#include "RDLog.h"
#include <iomanip>
#include <string>
#include <ctime>
#include <iostream>
#include <algorithm>
#include <sstream>

RDLogger rdAppLog = nullptr;
RDLogger rdDebugLog = nullptr;
RDLogger rdInfoLog = nullptr;
RDLogger rdErrorLog = nullptr;
RDLogger rdWarningLog = nullptr;
RDLogger rdStatusLog = nullptr;
namespace RDLog {

namespace {
const std::vector<RDLogger *> allLogs = {&rdAppLog,     &rdDebugLog,
                                         &rdInfoLog,    &rdErrorLog,
                                         &rdWarningLog, &rdStatusLog};
}

LogStateSetter::LogStateSetter() {
  for (auto i = 0u; i < allLogs.size(); ++i) {
    if (*allLogs[i] && (*allLogs[i])->df_enabled) {
      d_origState |= 1 << i;
      (*allLogs[i])->df_enabled = false;
    }
  }
}

LogStateSetter::LogStateSetter(RDLoggerList toEnable) : LogStateSetter() {
  for (auto i = 0u; i < allLogs.size(); ++i) {
    if (*allLogs[i] && std::find(toEnable.begin(), toEnable.end(),
                                 *allLogs[i]) != toEnable.end()) {
      d_origState ^= 1 << i;
      (*allLogs[i])->df_enabled = true;
    }
  }
}

LogStateSetter::~LogStateSetter() {
  for (auto i = 0u; i < allLogs.size(); ++i) {
    if (*allLogs[i]) {
      (*allLogs[i])->df_enabled ^= d_origState >> i & 1;
    }
  }
}
}  // namespace RDLog
namespace boost {
namespace logging {

void enable_logs(const char *arg) { enable_logs(std::string(arg)); };
void enable_logs(const std::string &arg) {
  // Yes... this is extremely crude
  if (arg == "rdApp.debug" || arg == "rdApp.*") {
    if (rdDebugLog) {
      rdDebugLog->df_enabled = true;
    }
  }
  if (arg == "rdApp.info" || arg == "rdApp.*") {
    if (rdInfoLog) {
      rdInfoLog->df_enabled = true;
    }
  }
  if (arg == "rdApp.warning" || arg == "rdApp.*") {
    if (rdWarningLog) {
      rdWarningLog->df_enabled = true;
    }
  }
  if (arg == "rdApp.error" || arg == "rdApp.*") {
    if (rdErrorLog) {
      rdErrorLog->df_enabled = true;
    }
  }
};
void disable_logs(const char *arg) { disable_logs(std::string(arg)); };
void disable_logs(const std::string &arg) {
  // Yes... this is extremely crude
  if (arg == "rdApp.debug" || arg == "rdApp.*") {
    if (rdDebugLog) {
      rdDebugLog->df_enabled = false;
    }
  }
  if (arg == "rdApp.info" || arg == "rdApp.*") {
    if (rdInfoLog) {
      rdInfoLog->df_enabled = false;
    }
  }
  if (arg == "rdApp.warning" || arg == "rdApp.*") {
    if (rdWarningLog) {
      rdWarningLog->df_enabled = false;
    }
  }
  if (arg == "rdApp.error" || arg == "rdApp.*") {
    if (rdErrorLog) {
      rdErrorLog->df_enabled = false;
    }
  }
};

bool is_log_enabled(RDLogger log) { return log && log->df_enabled; }

void get_log_status(std::ostream &ss, const std::string &name, RDLogger log) {
  ss << name << ":";
  if (log) {
    if (log->df_enabled) {
      ss << "enabled";
    } else {
      ss << "disabled";
    }
  } else {
    ss << "unitialized";
  }
}

std::string log_status() {
  std::stringstream ss;
  get_log_status(ss, "rdApp.debug", rdDebugLog);
  ss << std::endl;
  get_log_status(ss, "rdApp.info", rdInfoLog);
  ss << std::endl;
  get_log_status(ss, "rdApp.warning", rdWarningLog);
  ss << std::endl;
  get_log_status(ss, "rdApp.error", rdErrorLog);
  return ss.str();
}

}  // namespace logging
}  // namespace boost

namespace RDLog {
void InitLogs() {
  rdDebugLog = std::make_shared<boost::logging::rdLogger>(&std::cerr);
  rdDebugLog->df_enabled = false;
  rdInfoLog = std::make_shared<boost::logging::rdLogger>(&std::cout);
  rdInfoLog->df_enabled = false;
  rdWarningLog = std::make_shared<boost::logging::rdLogger>(&std::cerr);
  rdErrorLog = std::make_shared<boost::logging::rdLogger>(&std::cerr);
}

std::ostream &toStream(std::ostream &logstrm) {
  char buffer[16];
  time_t t = time(nullptr);
// localtime() is thread safe on windows, but not on *nix
#ifdef WIN32
  strftime(buffer, 16, "[%T] ", localtime(&t));
#else
  struct tm buf;
  strftime(buffer, 16, "[%T] ", localtime_r(&t, &buf));
#endif
  return logstrm << buffer;
}
}  // namespace RDLog
