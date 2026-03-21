#!/usr/bin/env bash

set -e

# venv有効化
# source venv/bin/activate

echo "Cleaning build directory..."
rm -rf build

echo "Creating build directory..."
mkdir build
cd build

echo "Running CMake..."
cmake ../nicole

echo "Building..."
cmake --build .

echo "Done!"
