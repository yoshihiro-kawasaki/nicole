# ソースファイル
file(GLOB_RECURSE SRC_FILES src/*.cpp)

# ライブラリ作成
add_library(nicole STATIC ${SRC_FILES})

# includeディレクトリ設定
target_include_directories(nicole
    PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include/
)

target_link_libraries(nicole PUBLIC odepack_cpp)

# 警告
# if (CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
#     target_compile_options(nicole_lib PRIVATE -Wall -Wextra -Wpedantic)
# endif()
