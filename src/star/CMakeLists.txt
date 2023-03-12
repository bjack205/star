add_library(
  star

  star.h
  typedefs.h

  vector3.c
  vector3.h

  quaternion.c
  quaternion.h
)
target_compile_definitions(star PUBLIC STAR_FLOAT=${STAR_FLOAT})
target_include_directories(star PUBLIC ${PROJECT_SOURCE_DIR}/src)

add_library(star::star ALIAS star)

add_library(
  star++

  star.hpp

  Vec3.cpp
  Vec3.hpp
)
target_link_libraries(star++ PRIVATE star::star)
target_include_directories(star++ PUBLIC ${PROJECT_SOURCE_DIR}/src)
add_library(star::star++ ALIAS star++)
