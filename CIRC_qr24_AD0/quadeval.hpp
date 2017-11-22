# ifndef QUADEVAL_HPP
# define QUADEVAL_HPP

#include <iostream>
#include <Eigen/Dense>
#include <cmath>
#include <algorithm>
#include <cstdlib>
#include <tuple>

#include "guassquad.hpp"

/////////////////////////////////////////////////
///////////// Parameters to modify //////////////
/////////////////////////////////////////////////
// Patch number is NPatch X NPatch
const int NPatch   = 10;

// Number of time steps
const int iteration = 5000;

// Regularizing parameter
const double beta = 3600;

// Size of advection
const double VELO  = 1.0;

/////////////////////////////////////////////////
/////////////////////////////////////////////////









const int qorder = 24;           // Quagrature rule for 1D 2D integral  
const int NB = 4;
typedef Eigen::Matrix<double,qorder,qorder> MatrixNd;
typedef Eigen::Matrix<double,qorder,1> VectorNd;
typedef Eigen::Matrix<double,NB,1> VectorInt;
typedef Eigen::Matrix<double,2*qorder,1> Vector2Nd;
typedef Eigen::Matrix<double,NB,3> MatrixPhi;
typedef Eigen::Matrix<double,NB,NB> MatrixInt2D;


const MatrixNd  MatrixNdZero  = Eigen::MatrixXd::Zero(qorder,qorder);
const Vector2Nd Vector2NdZero = Eigen::MatrixXd::Zero(2*qorder,1);


std::tuple< MatrixInt2D, MatrixInt2D, VectorInt, VectorInt>
    zero_tuple = std::make_tuple(Eigen::MatrixXd::Zero(NB,NB), Eigen::MatrixXd::Zero(NB,NB),
                                 Eigen::MatrixXd::Zero(NB,1),  Eigen::MatrixXd::Zero(NB,1));


const double PI = 3.1415926535897932384;
const double err_tol = 5*std::pow(10, -12);
const double epsilon = 1.0;
const double boollap = 1.0;

const double due   = 1.0;
const double tstep = due / iteration;     // Width of timestepping


const double scale = 0.998;
const double patch_radius = 1.0 / (NPatch-1);

const double xinf = 0.0;
const double xsup = 1.0;
const double yinf = 0.0;
const double ysup = 1.0;
const double CENX = 0.5;
const double CENY = 0.5;


// This function can be only used for refinement
int pointindomain(double xcoor, double ycoor) {
    if ( std::pow(xcoor-0.5, 2) + std::pow(ycoor-0.5, 2) < 0.25 + err_tol )
        return 1;
    else
        return 0;
}

// This function can be only used for refinement
int patchindomain(double xbinf, double xbsup, double ybinf, double ybsup) {
    int num = pointindomain(xbinf, ybinf) + pointindomain(xbsup, ybinf)
            + pointindomain(xbinf, ybsup) + pointindomain(xbsup, ybsup);
    return num;
}




// QuadRule qr;
// gaussq(18, qr);

// Legendre polynomials
// l_1(x) = 1
// l_1(x) = x
// l_1(x) = (3 * x^2 - 1) / 2
// at Gaussian quadrature points


VectorNd ZeroVecNd ( (VectorNd() << Eigen::VectorXd::Zero(qorder) ).finished() );
Vector2Nd ZeroVec2Nd ( (Vector2Nd() << Eigen::VectorXd::Zero(2*qorder) ).finished() );

VectorNd OneVecNd( (VectorNd() << Eigen::VectorXd::Ones(qorder) ).finished() );
Vector2Nd OneVec2Nd( (Vector2Nd() << Eigen::VectorXd::Ones(2*qorder) ).finished() );

MatrixNd OneMatNd( (MatrixNd() << Eigen::MatrixXd::Ones(qorder,qorder) ).finished() );

Vector2Nd OneZeroVec2Nd( (Vector2Nd() << Eigen::VectorXd::Ones(qorder), Eigen::VectorXd::Zero(qorder) ).finished() );
Vector2Nd ZeroOneVec2Nd( (Vector2Nd() << Eigen::VectorXd::Zero(qorder), Eigen::VectorXd::Ones(qorder) ).finished() );



VectorNd Node1d( (VectorNd() << -0.9951872199970219, -0.9747285559713099, -0.9382745520027325, -0.886415527004401, 
                                -0.8200019859739016, -0.7401241915785535, -0.6480936519369757, -0.5454214713888392, 
                                -0.4337935076260449, -0.3150426796961629, -0.1911188674736157, -0.06405689286260542, 
                                0.06405689286260563, 0.1911188674736163, 0.3150426796961633, 0.4337935076260449, 
                                0.5454214713888396, 0.6480936519369751, 0.740124191578555, 0.8200019859739024, 
                                0.8864155270044016, 0.9382745520027327, 0.9747285559713118, 0.9951872199970235 ).finished() );
 
Vector2Nd NN1d( (Vector2Nd() << ( Node1d - OneVecNd )/2, ( Node1d + OneVecNd )/2 ).finished() );
Vector2Nd LineNode1( (Vector2Nd() << ( Node1d - OneVecNd )/2, Eigen::VectorXd::Zero(qorder) ).finished() );
Vector2Nd LineNode2( (Vector2Nd() << Eigen::VectorXd::Zero(qorder), ( Node1d + OneVecNd )/2 ).finished() );
// VectorNd LineNode1( (VectorNd() << -0.9815606342467184, -0.904117256370475, -0.7699026741943048, -0.5873179542866171, -0.3678314989981802, -0.1252334085114692, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 ).finished());
// VectorNd LineNode2( (VectorNd() << 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.1252334085114689, 0.3678314989981796, 0.5873179542866168, 0.769902674194304, 0.904117256370473, 0.9815606342467168 ).finished());
// VectorNd Node1d( (VectorNd() << qr.nodes ).finished());



// l_2(x,y) = 1
MatrixNd basis1( (MatrixNd() << OneVecNd * OneVecNd.transpose() ).finished());

// l_2(x,y) = x
MatrixNd basis2( (MatrixNd() << OneVecNd * Node1d.transpose() ).finished());

// l_2(x,y) = y
MatrixNd basis3( (MatrixNd() << basis2.transpose() ).finished());

// l_2(x,y) = x*y
MatrixNd basis4( (MatrixNd() << basis2.cwiseProduct(basis3) ).finished());





// nodes & weights of 1D and 2D Guassian quadrature rules

VectorNd weight1d ( ( VectorNd() << 0.0123412297999884, 0.02853138862893307, 0.04427743881742088, 0.05929858491543619, 
                                    0.07334648141107943, 0.08619016153195282, 0.09761865210411416, 0.1074442701159656, 
                                    0.1155056680537255, 0.1216704729278029, 0.1258374563468277, 0.1279381953467519, 
                                    0.1279381953467522, 0.1258374563468285, 0.1216704729278034, 0.1155056680537255, 
                                    0.1074442701159656, 0.09761865210411301, 0.08619016153195372, 0.07334648141108006, 
                                    0.05929858491543541, 0.04427743881741945, 0.028531388628936, 0.01234122979998833 ).finished());
// VectorNd weight1d( (VectorNd() << qr.weights ).finished());
MatrixNd weight2d( (MatrixNd() <<  weight1d * weight1d.transpose() ).finished());

Vector2Nd WW1d( (Vector2Nd() <<  weight1d/2, weight1d/2 ).finished());

MatrixNd GridXNd( (MatrixNd() << basis2 ).finished());

MatrixNd GridYNd( (MatrixNd() << basis3 ).finished());



//const Eigen::MatrixXd L1 = { 1.0 };

// Function initial condition --- u0(x,y)
double initialcondition(double x, double y);



// Function of Velocity Field --- \beta(x,t)
std::pair<double, double> velocity(double x, double y, double t);

// Function f in oringinal equation --- f(x,y)
double funf(double x, double y);

// Bump function ( Special requirement for arguements --- Check below )
// double bump(double x, double y);


// Compute 1D equidistant mesh of size n between interval [a,b]
// Eigen::VectorXd equidistantmesh(int n, double a, double b);


///////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////


// Take point coordinate (x,y) as input, return u0(x,y)
double initialcondition(double x, double y) {

    double v = std::exp(x) * ( std::pow(x-0.5,2) + std::pow(y-0.5,2) - 0.25 );
    // double temp = 1;
    return v ;
}

double theo_solution(double x, double y) {
    if(std::pow(x-0.5, 2) + std::pow(y-0.5, 2) > 0.25)
        return 0;
    else {
        double v = std::exp(due+x) * ( std::pow(x-0.5,2) + std::pow(y-0.5,2) - 0.25 );
        return v;
    }
}

double bump1d(double x) {
    if(std::fabs(x) >= 1 - err_tol ) 
        return 0.0;
    else
        return std::exp( - 1 / ( 1 - std::pow(x,2) ) + 1 );
}

// Take point coordinate (x,y) and time t as input, return \beta(x,y,t)
std::pair<double, double> velocity(double x, double y) {
    //std::pair<double, double> v( -std::sin(x*PI) * std::cos(y*PI), std::sin(y*PI) * std::cos(x*PI) );

    std::pair<double, double> v(VELO, VELO);
    return v;
}


// Take point coordinate (x,y) as input, return f(x,y)
double funf1(double x, double y) {
    double v = std::exp(x) * ( VELO * ( std::pow(x-0.5,2) + std::pow(y-0.5,2) - 0.25 + 2*x + 2*y -2 ) - 4 * x - 2 );
    return v;
}



// Bump function take point coordinate (x,y)
// original coordinates must be affinely transformed to [-1,1] X [-1,1]
std::tuple<double, double, double> bump(double x, double y) {
    //assert( std::fabs(x) <= 1 && std::fabs(y) <= 1 && 
                 //"Coordinate must belong to [0,1] X [0,1] when doing bump function evaluation!");
    if(std::fabs(x) >= 1 - err_tol || std::fabs(y) >= 1 - err_tol ) 
        return std::make_tuple(0.0, 0.0, 0.0);
    else {
        double one_xsq, one_ysq, bmpv, dxbmp, dybmp;
        one_xsq = 1.0 - std::pow(x,2);
        one_ysq = 1.0 - std::pow(y,2);
        bmpv    = std::exp( - 1/one_xsq - 1/one_ysq + 2 );
        dxbmp   = - 2.0 * x * bmpv / std::pow(one_xsq,2);
        dybmp   = - 2.0 * y * bmpv / std::pow(one_ysq,2);
        return std::make_tuple(bmpv, dxbmp, dybmp);
    }
}



#endif
