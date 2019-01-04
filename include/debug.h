#ifndef MKE2_SRC_UTIL_DEBUG_H_
#define MKE2_SRC_UTIL_DEBUG_H_

#include <exception>
#include <iostream>

#include <fmt/ostream.h>

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

#define throw_fmt(exception, fmt_string, ...)                                 \
    do {                                                                      \
        throw exception(fmt::format("({}:{}) " fmt_string, __FILE__,          \
                                    __LINE__, __VA_ARGS__));                  \
    } while (0)

#define check_if(condition, fmt_string, ...)                                  \
    do {                                                                      \
        if (!(condition)) {                                                   \
            throw_fmt(std::runtime_error, fmt_string, __VA_ARGS__);           \
        }                                                                     \
    } while (0)

#define debug_fmt(fmt_string, ...)                                            \
    do {                                                                      \
        DBG                                                                   \
        {                                                                     \
            fmt::print(std::cout, fmt_string "\n", __VA_ARGS__);              \
        }                                                                     \
    } while (0)

#endif // MKE2_SRC_UTIL_DEBUG_H_
