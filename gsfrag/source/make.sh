#! /bin/sh
g++ -I/usr/include/x86_64-linux-gnu/c++/4.8 -static -O3 global.cpp graph.cpp strsrcSDF.cpp struct.cpp gsfrag.cpp -o gsfrag
g++  -I/usr/include/x86_64-linux-gnu/c++/4.8 -static -O3 global.cpp graph.cpp strsrcSDF.cpp struct.cpp gsfragl.cpp -o gsfragl
