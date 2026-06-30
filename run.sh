#!/bin/bash
clear

mkdir -p build  # create the build directory (-p: no error if it already exists)

g++ \
    apps/run_hysteresis.cpp \
    src/*.cpp \
    -I include/ \
    -o build/a.out \
&& ./build/a.out





# g++ \
#     apps/main.cpp \   # main source file with the entry point (main())
#     src/*.cpp \                 # all other source files in src/
#     -I include/ \               # write the compiled executable to build/a.out
#     -o build/a.out \            # write the compiled executable to build/a.out
# && ./build/a.out                # if compilation succeeds, run the executable



# -I include/ tells the compiler to add the include/ directory to the list of
# places it searches for header files. So when your source files do
# #include "something.hpp", the compiler will look inside include/ to find
# something.hpp. Without it, the compiler only searches the default system
# paths and the current directory, so it wouldn't find your project's headers
# in include/.