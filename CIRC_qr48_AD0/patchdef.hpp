# ifndef PATCHDEF_HPP
# define PATCHDEF_HPP

// #include <Eigen/Dense>
// using namespace Eigen;
#include <iostream>
#include <vector>
#include <cassert>
#include <Eigen/Dense>
#include <Eigen/SparseCore>

#include "quadeval.hpp"

// This class is designed for 2-dimensional cases
class Patch;

// Type to store coordinate of a single point in 2D-plane
typedef std::pair<double, double> PointC;

// Type to store coordinates of a vector of points in 2D-plane
typedef std::vector< PointC > Vec;

// transform absolute coordinates into local cooordinates of i(th) patch
std::pair<MatrixNd, MatrixNd> coortrans( MatrixNd x, MatrixNd y, Patch pat);

// transform local coordinates into absolute cooordinates of i(th) patch
std::pair<Vector2Nd, Vector2Nd> coortrans_v( Vector2Nd x, Vector2Nd y, double xc, double yc, double xu,double yu);

std::pair<MatrixNd, MatrixNd> coortrans_m( MatrixNd x, MatrixNd y, double xc, double yc, double xu,double yu);

std::pair<MatrixNd, MatrixNd> pullback( MatrixNd x, MatrixNd y );


///////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////


class Patch {

public:
    // Construct Patch by passing information on the boundary and evaluation points
    // IMPORTANT: coordinates in boundary should be ordered as following:
    // The x and y coordinates of first point must be STRICTLY SMALLER than the last point
    Patch(Vec boundary, Vec evapoints, std::pair<Vector2Nd, Vector2Nd> bquadp, Vector2Nd nx, Vector2Nd ny, double diaml );
    // Function to visit boundary data
    Vec get_bou();
    // Get dimension of local approximation space
    size_t get_dim();
    // Function to visit evaluation points
    Vec get_evap();

    Vector2Nd get_xline();
    Vector2Nd get_yline();
    Vector2Nd get_xnormv();
    Vector2Nd get_ynormv();
    MatrixNd  get_xgrid();
    MatrixNd  get_ygrid();
    // Print information of this patch
    void print_info();
    // Judge whether a point A belongs to Patch or not
    bool belong(PointC A);
    // Change x_line, y_line, x_norvec, y_norvec
    void changebound_quad(std::pair<Vector2Nd, Vector2Nd> quad, Vector2Nd nx, Vector2Nd ny);


    // x,y coordinates of 2D Guassian quadrature points
    MatrixNd x_grid;
    MatrixNd y_grid;

    // x,y coordinates of MAPPED 2D Guassian quadrature points
    MatrixNd x_mapped_grid;
    MatrixNd y_mapped_grid;

    // x,y coordinates of 1D Guassian quadrature points on boundary
    Vector2Nd x_line;
    Vector2Nd y_line;
    Vector2Nd x_norvec;
    Vector2Nd y_norvec;
    double diam_line;



private:
    // Boundary of this Patch, in the form of polygon, is a Vector of Points
    Vec bou;
    // Local Freedom of this Patch
    size_t locdim;
    // Evaluation points of this patch
    Vec evap;

};



///////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////



// IMPORTANT: coordinates in boundary should be ordered as following:
// The x and y coordinates of first point must be STRICTLY SMALLER than the last point
Patch::Patch(Vec boundary, Vec evapoints, std::pair<Vector2Nd, Vector2Nd> bquadp, Vector2Nd nx, Vector2Nd ny, double diaml ) {
    assert( boundary.size() == 4 && 
         "Patches are not rectangulars!");
    assert( boundary.front().first < boundary.back().first && boundary.front().second < boundary.back().second && 
         "Coordinates of Corners are not ordered!");
    bou = boundary;
    evap = evapoints;
    locdim = evapoints.size();

    x_grid = ( boundary.back().first - boundary.front().first )/2 * GridXNd + ( boundary.back().first + boundary.front().first )/2 * OneMatNd;
    y_grid = ( boundary.back().second - boundary.front().second )/2 * GridYNd + ( boundary.back().second + boundary.front().second )/2 * OneMatNd;

    for(int i = 0; i<qorder; i++) {
        for(int j = 0; j<qorder; j++) {
            std::pair<double, double> vtp0 = velocity(x_grid(i,j), y_grid(i,j));
            std::pair<double, double> vtp1 = velocity(x_grid(i,j) - vtp0.first * tstep, y_grid(i,j) - vtp0.second * tstep);
            x_mapped_grid(i,j) = x_grid(i,j) - (vtp0.first + vtp1.first) * tstep/2;
            y_mapped_grid(i,j) = y_grid(i,j) - (vtp0.second + vtp1.second) * tstep/2;

            // x_mapped_grid(i,j) = x_grid(i,j) - vtp0.first * tstep;
            // y_mapped_grid(i,j) = y_grid(i,j) - vtp0.second * tstep;
        }
    }

    x_line = bquadp.first;
    y_line = bquadp.second;
    x_norvec = nx;
    y_norvec = ny;
    diam_line = diaml;
}


Vec     Patch::get_bou()    { return bou; }
size_t  Patch::get_dim()    { return locdim; }
Vec     Patch::get_evap()   { return evap; }
Vector2Nd Patch::get_xline(){ return x_line; }
Vector2Nd Patch::get_yline(){ return y_line; }
Vector2Nd Patch::get_xnormv(){ return x_norvec; }
Vector2Nd Patch::get_ynormv(){ return y_norvec; }
MatrixNd  Patch::get_xgrid(){ return x_grid; }
MatrixNd  Patch::get_ygrid(){ return y_grid; }

void Patch::changebound_quad(std::pair<Vector2Nd, Vector2Nd> quad, Vector2Nd nx, Vector2Nd ny) {
    x_line = quad.first;
    y_line = quad.second;
    x_norvec = nx;
    y_norvec = ny;
}


bool Patch::belong(PointC A) {
    if ( A.first >= bou.front().first + err_tol && A.first <= bou.back().first - err_tol &&
    	 A.second >= bou.front().second + err_tol && A.second <= bou.back().second - err_tol )
        return true;
    else
    	return false;
}


void Patch::print_info() {
    std::cout << "------------------------------" << std::endl;
    std::cout << "Boundary of Patch:" << std::endl;
    for (auto i = this->bou.begin(); i != this->bou.end(); ++i)
    std::cout << i->first << ' ' << i->second << std::endl;

    std::cout << "Degree of Freedom of Patch:" << std::endl << locdim << std::endl;

    std::cout << "Evaluation Points of Patch:" << std::endl;
    for (auto i = this->evap.begin(); i != this->evap.end(); ++i)
        std::cout << i->first << ' ' << i->second << std::endl;
    /*
    std::cout << std::left << std::setw(16) << "Bou. Quad" << std::left << std::setw(16) << "Nor. Vec" << std::endl;
    for (int i = 0; i < x_line.size(); ++i) {
        std::cout << std::left << "(" << std::setprecision(4) << std::setw(4) << x_line(i) << "," << std::setprecision(4) << std::setw(4) << y_line(i) << ")";
        std::cout << std::left << "(" << std::setprecision(4) << std::setw(4) << x_norvec(i) << "," << std::setprecision(4) << std::setw(4) << y_norvec(i) << ")";
        std::cout << std::endl;
    }
    */

    /*
    std::cout << "Mesh grid for numerical integral:" << std::endl;
    for (int i = 0; i < x_grid.rows(); ++i) {
        for (int j = 0; j < x_grid.rows(); ++j) {
            std::cout << "(" << x_grid(i,j) << "," << y_grid(i,j) << ") ";
        }
        std::cout << std::endl;
    }

    std::cout << "Mapped Mesh grid for numerical integral:" << std::endl;
    for (int i = 0; i < x_mapped_grid.rows(); ++i) {
        for (int j = 0; j < x_mapped_grid.rows(); ++j) {
            std::cout << "(" << x_mapped_grid(i,j) << "," << y_mapped_grid(i,j) << ") ";
        }
        std::cout << std::endl;
    }
    */
}


std::pair<MatrixNd, MatrixNd> coortrans( MatrixNd x, MatrixNd y, Patch pat) {
    std::pair<MatrixNd, MatrixNd> tmp;
    tmp.first = (x - OneMatNd * pat.get_evap().front().first) / (pat.get_evap().front().first - pat.get_bou().front().first);
    tmp.second = (y - OneMatNd * pat.get_evap().front().second) / (pat.get_evap().front().second - pat.get_bou().front().second);
    return tmp;
}

std::pair<MatrixNd, MatrixNd> coortrans_m( MatrixNd x, MatrixNd y, double xc, double yc, double xu,double yu) {
	// reference interval to real interval
	// Known: center and upper boundary point
    std::pair<MatrixNd, MatrixNd> tmp;
    tmp.first = x * (xu - xc) + OneMatNd * xc;
    tmp.second = y * (yu - yc) + OneMatNd * yc;
    return tmp;
}

std::pair<Vector2Nd, Vector2Nd> coortrans_v( Vector2Nd x, Vector2Nd y, double xc, double yc, double xu, double yu) {
	// reference interval to real interval
	// Known: center and upper boundary point
    std::pair<Vector2Nd, Vector2Nd> tmp;
    tmp.first = x * (xu - xc) + OneVec2Nd * xc;
    tmp.second = y * (yu - yc) + OneVec2Nd * yc;
    return tmp;
}

std::pair<MatrixNd, MatrixNd> pullback( MatrixNd x_grid, MatrixNd y_grid ) {
	MatrixNd x_mapped_grid;
    MatrixNd y_mapped_grid;
	for(int i = 0; i<qorder; i++) {
        for(int j = 0; j<qorder; j++) {
            std::pair<double, double> vtp0 = velocity(x_grid(i,j), y_grid(i,j));
            std::pair<double, double> vtp1 = velocity(x_grid(i,j) - vtp0.first * tstep, y_grid(i,j) - vtp0.second * tstep);
            x_mapped_grid(i,j) = x_grid(i,j) - (vtp0.first + vtp1.first) * tstep/2;
            y_mapped_grid(i,j) = y_grid(i,j) - (vtp0.second + vtp1.second) * tstep/2;

            // x_mapped_grid(i,j) = x_grid(i,j) - vtp0.first * tstep;
            // y_mapped_grid(i,j) = y_grid(i,j) - vtp0.second * tstep;
        }
    }
    return std::pair<MatrixNd, MatrixNd>( x_mapped_grid, y_mapped_grid );
}




# endif
