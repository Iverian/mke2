enable_testing()


file(GLOB TEST_SOURCES "src/test_*.cpp")

###################################################
# Make all tests in single executable

add_executable(all_tests
    "src/all_tests.cpp"
    ${TEST_SOURCES}
)

target_include_directories(all_tests
    PUBLIC
        "../src"
)

target_link_libraries(all_tests
    PUBLIC
        ${LIBPROXY}
        CONAN_PKG::gtest
)

add_test(all_tests all_tests)
###################################################
