# function add_star_test(name)
#
# Adds a new executable test called <name>_test.
# Assumes thes source code is in a file called <name>_test.
function (add_star_test name)
  set(TEST_NAME ${name}_test)
  add_executable(${TEST_NAME}
    ${TEST_NAME}.cpp
    )
  target_link_libraries(${TEST_NAME}
    PRIVATE
    star::star++
    gtest::gtest
    )
  if (NOT APPLE AND NOT WIN32)
    target_link_libraries(${TEST_NAME} PUBLIC m)
  endif()
  gtest_discover_tests(${TEST_NAME})
endfunction()

add_star_test(vector3)
add_star_test(matrix3)
add_star_test(matrix)
add_star_test(quaternion)
add_star_test(quaternion_class)
add_star_test(matrix_class)
add_star_test(rotmat_class)

add_executable(vector3 vector3_main.c)
target_link_libraries(vector3 PRIVATE star::star)