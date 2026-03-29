#!/usr/bin/env bash

cpp_file="${1}"

if [ ! -e "${cpp_file}" ]; then
    echo "ファイルが指定されていません: ${cpp_file}"
    exit 1
fi

g++ -O3 "${cpp_file}" \
    -I ../../nicole/include \
    -L ../../build/nicole -lnicole \
    -I ../../odepack_cpp/include \
    -L ../../build/odepack_cpp -lodepack_cpp \
    -o run
