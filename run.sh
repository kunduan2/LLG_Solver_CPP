#!/bin/bash
clear

mkdir -p build

g++ \
    apps/main.cpp \
    src/*.cpp \
    -I include/ \
    -o build/a.out \
&& ./build/a.out

