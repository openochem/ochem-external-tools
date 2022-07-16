#! /bin/sh
g++ -static -O3 global.cpp graph.cpp strsrcSDF.cpp struct.cpp gsfrag.cpp -o gsfrag
g++ -static -O3 global.cpp graph.cpp strsrcSDF.cpp struct.cpp gsfragl.cpp -o gsfragl
