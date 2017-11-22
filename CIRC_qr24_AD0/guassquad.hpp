# ifndef GUASSQUAD_HPP
# define GUASSQUAD_HPP

#include <iostream>
#include <cassert>
#include <Eigen/Dense>


// Structure for storing Weights and Nodes for a Quad rule
struct QuadRule {
  Eigen::VectorXd nodes, weights;
};

// Compute Weights and Nodes for 1D GaussQuad
void gaussq ( const unsigned n , QuadRule& qr );



///////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////


void gaussq ( const unsigned n , QuadRule& qr ) {
  // Initialize bidiagonal matrix Jn 
  Eigen::MatrixXd M = Eigen::MatrixXd::Zero(n, n);
  for(unsigned i = 1; i < n; ++i) {
    const double b = i / std::sqrt(4. * i * i - 1.);
    M(i, i-1) = b;
    M(i-1, i) = b;
  }

  // EIGEN’s built-in helper class for eigenvalue problems
  // (use method for symmetric matrices, exploiting the structure) 
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eig(M);

  qr.nodes = eig.eigenvalues();
  qr.weights = 2 * eig.eigenvectors().topRows<1>().array().pow(2);
}


void qrprint ( const QuadRule& qr ) {
  std::cout << "------------------------------------" << std::endl
            << "Coordinates of Nodes:" << std::endl
            << qr.nodes << std::endl
            << "------------------------------------" << std::endl
            << "Information about Weights:" << std::endl
            << qr.weights << std::endl ;
}

/*

// Affine transformation for evaluating 2D numerical integral
std::pair<PMatrix, double> evaquad2d ( const std::vector<double> a,
                                       const unsigned order, const QuadRule qr) {

  assert( a.size() == 4 && a[0] < a[1] && a[2] < a[3] && 
         "Lower boundary of 2D interval must be smaller than upper boundary of interval!");

  // Affine Transformation of Coordinates in 2D
  double c1 = (a[0] + a[1])/2, s1 = (a[1] - a[0])/2, c2 = (a[2] + a[3])/2, s2 = (a[3] - a[2])/2;
  Eigen::VectorXd unitv = Eigen::VectorXd::Ones(order);
  Eigen::VectorXd coor1_t = qr.nodes * s1 + c1 * unitv;
  Eigen::VectorXd coor2_t = qr.nodes * s2 + c2 * unitv;

  double qrint = 0;
  PMatrix pm;
  pm.reserve(order);

  for(int i = 0; i < order; ++i) {
    pm.push_back( std::vector< PointC >() );
    for(int j = 0; j < order; ++j) {
    	pm.at(i).push_back(PointC(coor1_t(i), coor2_t(j)));
    }
  }

  return std::make_pair( pm, s1*s2 );
}
*/

#endif