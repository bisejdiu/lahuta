#ifndef PYTHONCAPI_COMPAT_H
#define PYTHONCAPI_COMPAT_H

#include <Python.h>

#ifdef __cplusplus
extern "C" {
#endif

// See: https://github.com/python/pythoncapi-compat/blob/90c06a4cae557bdbfa4f231a781d2b5c1a8f6d1c/pythoncapi_compat.h#L867
#if (PY_VERSION_HEX >= 0x030A0000) && (PY_VERSION_HEX < 0x030D0000)
extern int _Py_IsFinalizing(void);
static inline int Py_IsFinalizing(void) {
    return _Py_IsFinalizing();
}
#endif

#ifdef __cplusplus
}
#endif

#endif // PYTHONCAPI_COMPAT_H
