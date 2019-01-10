#ifndef MKE2_INCLUDE_DEBUG_H_
#define MKE2_INCLUDE_DEBUG_H_

#include <cstdio>
#include <exception>
#include <iostream>

#ifdef NDEBUG
#define DEBUG_FLAG false
#define DBG if (false)
#else
#define DEBUG_FLAG true
#define DBG if (true)
#endif

#define cdbg                                                                  \
    if (!DEBUG_FLAG) {                                                        \
    } else                                                                    \
        std::cout

#define throw_fmt(fmt_string, ...)                                            \
    do {                                                                      \
        char what[BUFSIZ];                                                    \
        snprintf(what, BUFSIZ, "(%s:%d)" fmt_string, __FILE__, __LINE__,      \
                 __VA_ARGS__);                                                \
        throw std::runtime_error(what);                                       \
    } while (0)

#define check_if(condition, fmt_string, ...)                                  \
    do {                                                                      \
        if (!(condition)) {                                                   \
            throw_fmt(fmt_string, __VA_ARGS__);                               \
        }                                                                     \
    } while (0)

#define debug_fmt(fmt_string, ...)                                            \
    do {                                                                      \
        if (DEBUG_FLAG) {                                                     \
            fprintf(stdout, fmt_string, __VA_ARGS__);                         \
        }                                                                     \
    } while (0)

#endif // MKE2_INCLUDE_DEBUG_H_
