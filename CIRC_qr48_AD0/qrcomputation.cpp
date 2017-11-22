#include <iostream>
#include <fstream>
#include <iomanip>

// g++ main.cpp -I./../../eigenlib -o main

#include <Eigen/Dense>
#include "guassquad.hpp"

int main() {
	int order = 24;
	QuadRule qr;
	gaussq(order, qr);
	std::cout << "NODES:" << qr.nodes.size() << std::endl;
	for (int i = 0; i < order; ++i) {
		std::cout << std::setprecision(16) << qr.nodes(i) << ", ";
	}

	std::cout << std::endl << std::endl << "WEIGHTS:" << std::endl;
	for (int i = 0; i < order; ++i) {
		std::cout << std::setprecision(16) << qr.weights(i) << ", ";
	}
	return 0;
}