//
// Copyright (C)  2005-2022 Greg Landrum and other RDKit contributors
//
//  @@ All Rights Reserved @@
//  This file is part of the RDKit.
//  The contents are covered by the terms of the BSD license
//  which is included in the file license.txt, found at the root
//  of the RDKit source tree.
//

#ifndef RDLOG_H_29JUNE2005
#define RDLOG_H_29JUNE2005

#include "SimpleTee.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <cstdint>

namespace boost {
namespace logging {

class rdLogger {
 public:
  std::ostream *dp_dest;
  bool df_owner;
  bool df_enabled;

  std::ofstream *dp_teeHelperStream;
  TeeStream *teestream;

  rdLogger(std::ostream *dest, bool owner = false)
      : dp_dest(dest),
        df_owner(owner),
        df_enabled(true),
        dp_teeHelperStream(nullptr),
        teestream(nullptr) {}

  //! Sets a stream to tee the output to.
  void SetTee(std::ostream &stream) {
    if (dp_dest) {
      ClearTee();
      teestream = new TeeStream(*dp_dest, stream);
    }
  }

  //! Sets a filename to tee the output to.
  void SetTee(const char *filename) {
    if (dp_dest) {
      auto s = new std::ofstream(filename);
      if (!s->is_open()) {
        delete s;
        return;
      }
      SetTee(*s);
      dp_teeHelperStream = s;
    }
  }

  //! Sets a filename to tee the output to.
  void SetTee(const std::string &filename) { return SetTee(filename.c_str()); }

  //! Remove our tee if it's set.
  void ClearTee() {
    if (dp_dest) {
      delete teestream;
      teestream = nullptr;
      if (dp_teeHelperStream) {
        dp_teeHelperStream->close();
        delete dp_teeHelperStream;
        dp_teeHelperStream = nullptr;
      }
    }
  }
  ~rdLogger() {
    if (dp_dest) {
      dp_dest->flush();
      ClearTee();
      if (df_owner) {
        delete dp_dest;
      }
      dp_dest = nullptr;
    }
  }

 private:
  // disable copy ctor and assignment
  rdLogger(const rdLogger &);
  rdLogger &operator=(const rdLogger &);
};
void enable_logs(const char *arg);
void enable_logs(const std::string &arg);
void disable_logs(const char *arg);
void disable_logs(const std::string &arg);
std::string log_status();
}  // namespace logging
}  // namespace boost
namespace RDLog {
std::ostream &toStream(std::ostream &);
}
#define BOOST_LOG(__arg__)                                      \
  if ((__arg__) && (__arg__->dp_dest) && (__arg__->df_enabled)) \
  RDLog::toStream((__arg__->teestream) ? *(__arg__->teestream)  \
                                       : *(__arg__->dp_dest))

using RDLogger = std::shared_ptr<boost::logging::rdLogger>;

extern RDLogger rdAppLog;
extern RDLogger rdDebugLog;
extern RDLogger rdInfoLog;
extern RDLogger rdErrorLog;
extern RDLogger rdWarningLog;
extern RDLogger rdStatusLog;

namespace RDLog {
void InitLogs();

using RDLoggerList = std::vector<RDLogger>;
class LogStateSetter {
 public:
  //! enables only the logs in the list, the current state will be restored when
  //! this object is destroyed
  LogStateSetter(RDLoggerList toEnable);
  //! disables all logs, the current state will be restored when this object is
  //! destroyed
  LogStateSetter();
  ~LogStateSetter();

  LogStateSetter(const LogStateSetter&) = delete;
  LogStateSetter& operator=(const LogStateSetter&) = delete;

 private:
  std::uint64_t d_origState = 0;
};

inline void deprecationWarning(const std::string& message) {
  BOOST_LOG(rdWarningLog) << "DEPRECATION WARNING: " << message << std::endl;
}

}  // namespace RDLog
#endif
