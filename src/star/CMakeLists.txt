add_library(
  star

  star.h
  typedefs.h

  vector3.c
  vector3.h

  vector4.c
  vector4.h

  matrix3.c
  matrix3.h

  matrix4.c
  matrix4.h

  quaternion.c
  quaternion.h

  matrix4.c matrix4.h

  matrix43.c matrix43.h
)
target_compile_definitions(star PUBLIC STAR_FLOAT=${STAR_FLOAT})
target_include_directories(star PUBLIC ${PROJECT_SOURCE_DIR}/src)

add_library(star::star ALIAS star)

add_library(
  star++

  star.hpp

  Vec3.cpp
  Vec3.hpp

  Vec4.cpp
  Vec4.hpp

  Quaternion.cpp
  Quaternion.hpp

  Mat3.cpp
  Mat3.hpp

  Transpose.cpp
  Transpose.hpp

  matrix_multiplication.cpp
  matrix_multiplication.hpp
  Mat4.cpp Mat4.hpp Mat43.cpp Mat43.hpp RotMat.hpp)
target_link_libraries(star++ PUBLIC star::star)
target_include_directories(star++ PUBLIC ${PROJECT_SOURCE_DIR}/src)
add_library(star::star++ ALIAS star++)
