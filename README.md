# PUFEM-for-advection-diffusion
PUFEM implementation for advection diffusion equation (beta)s

Prerequisites
=========================

1. gcc (GNU or Clang)
2. Eigen
3. python

Procedures: Ubuntu & Mac OS X (Yosemite / El Capitan / Sierra)
=========================

1. Please install gcc from GNU or Clang.
2. Please specify location of Eigen in makefile.
3. then type `make` and excute the program.
4. If you want to visualize the error or result, you can compile and run `plot.py` to get pictures.
5. (IMPORTANT) If you want to test code with different parameters, please go to `quadeval.cpp` file.

Cases
=========================

1. Two folders **CIRC_qr24_AD0** **CIRC_qr48_AD0** are independent. **CIRC_qr24_AD0** is for quadrature rule of order 24. **CIRC_qr48_AD0** is for quadrature rule of order 48. For tests in both cases, please do `make` and compilation in corresponding folders.

Remarks
=========================
1. This code is not quite stable now. That means it may not work in extreme cases such as very few patches or very large matrices. Some inappropriate choices of parameters might lead to errors.
2. The code for square case is being fixed. Due to the difference of compilers, my implementation for square case works for gcc from Clang but it goes wrong under gcc from GNU.
3. **If you have any suggestions or questions, please drop me an email.**
