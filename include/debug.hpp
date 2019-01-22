#ifndef MKE2_INCLUDE_DEBUG_HPP_
#define MKE2_INCLUDE_DEBUG_HPP_

#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <stdexcept>

#ifdef NDEBUG

#define NOEXCEPTD noexcept
static constexpr auto debug_flag = false;

#else // NDEBUG

#define NOEXCEPTD
static constexpr auto debug_flag = true;

#endif // NDEBUG

#define DBG if (debug_flag)

#define cerrd                                                                 \
    if (!debug_flag) {                                                        \
    } else                                                                    \
        std::cerr

#define coutd                                                                 \
    if (!debug_flag) {                                                        \
    } else                                                                    \
        std::cout

#define throw_fmt(fmt_string, ...)                                            \
    do {                                                                      \
        char what[BUFSIZ];                                                    \
        snprintf(what, BUFSIZ, "(%s:%d) " fmt_string, __FILE__, __LINE__,     \
                 ##__VA_ARGS__);                                              \
        throw std::runtime_error(what);                                       \
    } while (0)

#define check_if(condition, fmt_string, ...)                                  \
    do {                                                                      \
        if (!(condition)) {                                                   \
            throw_fmt(fmt_string, ##__VA_ARGS__);                             \
        }                                                                     \
    } while (0)

#ifdef NDEBUG

#define debug_fmt(fmt_string, ...)
#define check_ifd(condition, fmt_string, ...)

#else // NDEBUG

#define debug_fmt(fmt_string, ...)                                            \
    do {                                                                      \
        fprintf_s(stderr, fmt_string, ##__VA_ARGS__);                         \
    } while (0)

#define check_ifd(condition, fmt_string, ...)                                 \
    check_if(condition, fmt_string, ##__VA_ARGS__)

#endif // NDEBUG

#endif // MKE2_INCLUDE_DEBUG_HPP_
