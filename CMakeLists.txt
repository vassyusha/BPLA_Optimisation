cmake_minimum_required(VERSION 3.10)

project(BPLA_Optimisation VERSION 1.0 LANGUAGES CXX)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build" FORCE)
endif()

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Поиск OpenMP
find_package(OpenMP REQUIRED)

# Явно перечисляем файлы. Если добавите новый .cpp — просто допишите его сюда.
set(SOURCES 
    main.cpp
)

# Если у вас есть папка include, добавляем её
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)

add_executable(${PROJECT_NAME} ${SOURCES})

if(OpenMP_CXX_FOUND)
    target_link_libraries(${PROJECT_NAME} PUBLIC OpenMP::OpenMP_CXX)
endif()