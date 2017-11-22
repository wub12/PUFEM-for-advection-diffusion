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
const int NPatch   = 40;

// Number of time steps
const int iteration = 5000;

// Regularizing parameter
const double beta = 3600;

// Size of advection
const double VELO  = 10.0;

/////////////////////////////////////////////////
/////////////////////////////////////////////////









const int qorder = 48;           // Quagrature rule for 1D 2D integral  
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
    if( std::pow(xcoor-0.5, 2) + std::pow(ycoor-0.5, 2) < 0.25 + err_tol )
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



VectorNd Node1d( (VectorNd() << -0.9987710072524263, -0.993530172266351, -0.9841245837228263, -0.9705915925462468, 
                                -0.9529877031604305, -0.9313866907065541, -0.9058791367155687, -0.8765720202742467, 
                                -0.8435882616243935, -0.8070662040294428, -0.7671590325157407, -0.7240341309238134, 
                                -0.6778723796326642, -0.6288673967765143, -0.5772247260839727, -0.5231609747222333, 
                                -0.4669029047509585, -0.4086864819907171, -0.348755886292161, -0.2873624873554558, 
                                -0.2247637903946886, -0.1612223560688914, -0.09700469920946253, -0.03238017096286944, 
                                0.0323801709628693, 0.09700469920946292, 0.1612223560688918, 0.2247637903946888, 
                                0.2873624873554562, 0.3487558862921614, 0.408686481990717, 0.4669029047509586, 
                                0.5231609747222339, 0.5772247260839738, 0.6288673967765136, 0.6778723796326634, 
                                0.7240341309238152, 0.7671590325157405, 0.8070662040294427, 0.843588261624395, 
                                0.8765720202742491, 0.9058791367155715, 0.9313866907065548, 0.952987703160432, 
                                0.9705915925462469, 0.9841245837228277, 0.9935301722663513, 0.998771007252427 ).finished() );
 
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

VectorNd weight1d ( ( VectorNd() << 0.003153346052308471, 0.007327553901276586, 0.01147723457923283, 0.01557931572294468, 
                                    0.01961616045735411, 0.02357076083932367, 0.02742650970835679, 0.03116722783279804, 
                                    0.03477722256477041, 0.0382413510658313, 0.04154508294346439, 0.04467456085669397, 
                                    0.04761665849249046, 0.05035903555385458, 0.05289018948519328, 0.05519950369998431, 
                                    0.05727729210040331, 0.05911483969839584, 0.06070443916589393, 0.06203942315989253, 
                                    0.06311419228625466, 0.0639242385846481, 0.06446616443595006, 0.06473769681268396, 
                                    0.06473769681268372, 0.06446616443594988, 0.06392423858464863, 0.06311419228625326, 
                                    0.06203942315989219, 0.06070443916589538, 0.05911483969839662, 0.05727729210040256, 
                                    0.05519950369998501, 0.05289018948519544, 0.05035903555385333, 0.04761665849248822, 
                                    0.0446745608566946, 0.04154508294346598, 0.03824135106582855, 0.03477722256477007, 
                                    0.03116722783279783, 0.02742650970836087, 0.02357076083932251, 0.01961616045735541, 
                                    0.01557931572294324, 0.01147723457923388, 0.007327553901277332, 0.003153346052307157 ).finished());
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
