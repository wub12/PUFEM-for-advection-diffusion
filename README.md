# PUFEM-for-advection-diffusion
PUFEM implementation for advection diffusion equation (beta)s

Prerequisites
=================

1. gcc (GNU or Clang)
2. Eigen
3. python

Linux & Mac OS X (Yosemite / El Capitan / Sierra)
---------------------

1. Please install gcc from GNU or Clang
2. Please specify location of Eigen in makefile
3. then type $make$
4. If you want to visualize the error or result, you can compile and run plot.py to get pictures.

Cases:

1. Two folders "CIRC_qr24_AD0" "CIRC_qr48_AD0" are independent. "CIRC_qr24_AD0" is for quadrature rule of order 24. "CIRC_qr48_AD0" is for quadrature rule of order 48. For tests in both cases, please do $make$ and compilation in corresponding folders.
