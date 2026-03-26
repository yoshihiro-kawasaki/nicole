Calculation of chemical reaction and non-ideal MHD resistivity for astropyscis

odepack_cppの取得（更新）
```bash
./scripts/clone_odepack_cpp.sh
```

```bash
./build.sh
g++ -O3 test.cpp -I ./nicole/include -L ./build/nicole -lnicole -I ./odepack_cpp/include -L ./build/odepack_cpp -lodepack_cpp
```