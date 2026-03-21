set(CMAKE_C_COMPILER /opt/homebrew/bin/gcc)
set(CMAKE_CXX_COMPILER /opt/homebrew/bin/g++)

cmake_minimum_required(VERSION 3.18)

project(nicole)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# 再帰的に.cppを取得
file(GLOB_RECURSE SRC_FILES
    src/*.cpp
)

# ライブラリ作成
add_library(nicole ${SRC_FILES})

# includeディレクトリ設定
target_include_directories(nicole
    PUBLIC
        ${PROJECT_SOURCE_DIR}/include
)

# 警告（任意・おすすめ）
if (CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
    target_compile_options(nicole PRIVATE -Wall -Wextra -Wpedantic)
endif()
