#ifndef MKE2_INCLUDE_DEBUG_HPP_
#define MKE2_INCLUDE_DEBUG_HPP_

#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <stdexcept>

#ifdef NDEBUG
#define DEBUG_FLAG false
#define DBG if (false)
#else
#define DEBUG_FLAG true
#define DBG if (true)
#endif

#define cerrd                                                                 \
    if (!DEBUG_FLAG) {                                                        \
    } else                                                                    \
        std::cerr

#define coutd                                                                 \
    if (!DEBUG_FLAG) {                                                        \
    } else                                                                    \
        std::cerr

#define exit_fmt(fmt_string, ...)                                             \
    do {                                                                      \
        fprintf_s(stderr, "CRITICAL: (%s:%d) " fmt_string, __FILE__,          \
                  __LINE__, ##__VA_ARGS__);                                   \
        std::exit(1);                                                         \
    } while (0)

#define check_if(condition, fmt_string, ...)                                  \
    do {                                                                      \
        if (!(condition)) {                                                   \
            exit_fmt(fmt_string, ##__VA_ARGS__);                              \
        }                                                                     \
    } while (0)

#define debug_fmt(fmt_string, ...)                                            \
    do {                                                                      \
        if (DEBUG_FLAG) {                                                     \
            fprintf_s(stderr, fmt_string, ##__VA_ARGS__);                     \
        }                                                                     \
    } while (0)

#endif // MKE2_INCLUDE_DEBUG_HPP_
