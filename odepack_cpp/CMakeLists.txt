# ソースファイル
file(GLOB_RECURSE SRC_FILES src/*.cpp)

# 静的ライブラリ作成
add_library(odepack_cpp STATIC ${SRC_FILES})

# includeディレクトリ設定
target_include_directories(odepack_cpp
    PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include/
)

# 警告
# if (CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
#     target_compile_options(odepack_cpp PRIVATE -Wall -Wextra -Wpedantic)
# endif()
