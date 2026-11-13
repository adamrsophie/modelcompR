// Stub of HiGHS io/HighsIO.h for the feasibility jump harness: logging is a no-op.
#ifndef FJ_STUB_HIGHSIO_H_
#define FJ_STUB_HIGHSIO_H_
#include "lp_data/HConst.h"
struct HighsLogOptions {};
inline void highsLogDev(const HighsLogOptions&, const HighsLogType, const char*, ...) {}
#endif
