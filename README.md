mathResearch
==================

A repository for code for math research. As of right now, it holds code for CUDA operations on matrices done in python, as well as code for optimized demazure operators in Sage.

Run in the demazure directory:
    make

This will generate a .so, .c, and .pyx file so that opGroup.spyx can reference the code from modPolyCy.spyx. To remove these files, run in the demazure directory:
    make clean
