#
# odepack_cppの取得（更新）
# ./scripts/clone_odepack_cpp.sh
#

#!/usr/bin/env bash

set -e

if [ -d odepack_cpp ]; then
    rm -r odepack_cpp
fi

mkdir -p tmp
cd tmp
git clone https://github.com/yoshihiro-kawasaki/odepack_cpp.git
mv odepack_cpp/odepack_cpp ../
cd ..
rm -rf tmp
