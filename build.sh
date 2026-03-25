#!/usr/bin/env bash

# エラーが発生したら即座に終了
set -e

# 古いビルドディレクトリの削除
echo "Cleaning build directory..."
rm -rf build

# CMakeの構成
# -S . : ソースは現在のディレクトリ
# -B build : ビルド用の中間ファイルはbuildディレクトリへ
# -D... : リリースモード（-O3など）を指定
echo "Running CMake..."
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release

# ビルド
# --build build : build ディレクトリ内の設定に基づいてコンパイルを実行
echo "Building..."
cmake --build build

echo "Done!"
