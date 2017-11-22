# ifndef RECTDOMAIN_HPP
# define RECTDOMAIN_HPP

#include <iostream>
#include <vector>
#include <cassert>
#include "patchdef.hpp"
#include "quadeval.hpp"


// Here DEFINE CLASS Patchset containing :
// Information of all patches
// Overlapping relationship between all patches
// Necessary functions to do integral evaluation tasks

// RESTRICTIONS:
// 1. Patches must be rectangles
// 2. Patches must have some heights and widths

class PatchSet {

public:

    // Constrct Patchset by passing a list of patches and a list of corresponding information
    PatchSet(std::vector<Patch> pt, std::vector< std::vector<int> > rela, std::vector< bool > onboundary);

    // Function to get length of PatchSet
    int getLeng();

    // Function to get Vector of Patches
    std::vector<Patch> gpt();

    // Function to get length of PatchSet
    std::vector< std::vector<int> > gre();

    // Function to get On Boundary information
    std::vector< bool > getonbound();

    // Function Phi (Particle of Unity) --- phi_i(x,y) --- (x,y) absolute coordinate
    std::tuple<double, double, double> PhiPhiD(double x, double y, int i);

    MatrixPhi PhiMat(double x, double y, int i);

    // Function to evaluate numetical integral --- integral( u * v ) over SIGMA
    void QEVAL_uv ( Eigen::SparseMatrix<double> &M, 
    				Eigen::SparseMatrix<double> &M_mapped, 
    				Eigen::VectorXd &rhs, 
    				Eigen::VectorXd &frhs);

    void evaluate( Eigen::VectorXd c, int Ntest, Eigen::MatrixXd &r );

//   // Function to evaluate numetical integral --- integral( grand(u) * grad(v) ) over SIGMA
//   double QEVAL_gugv(int i, int j);

//   // Function to evaluate numetical integral --- integral( u * grad(v) * n ) over boundary(whole DOMAIN)
//   double QEVAL_ugv(int i, int j);
    
    // Function to Print ALL Information of PatchSet
    void PS_print();

private:
    // Number of Patches
    int Leng;
    int sqrtleng;
    // Patches
    std::vector<Patch> VEC_P;
    // Relation vector: i(th) vector stores such elements that overlap with i(th) patch
    std::vector< std::vector<int> > RELA_P;
    // Check if on the boundary: i(th) element stores bool if this element is on the boundary
    std::vector< bool > OnBOUN;
    // Radius of square patch
    double radius;

    std::tuple < MatrixInt2D, MatrixInt2D, VectorInt, VectorInt > 
        numericalint2d ( double xbinf, double xbsup, double ybinf, double ybsup, int i, int j );

    std::tuple < MatrixInt2D, MatrixInt2D, VectorInt, VectorInt > 
        boundary_numericalint2d( double xbinf, double xbsup, double ybinf, double ybsup, int i, int j );

    MatrixInt2D numericalint1d ( int i, int j );

};



///////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////




PatchSet::PatchSet( std::vector<Patch> pt, 
					std::vector< std::vector<int> > rela, 
					std::vector< bool > onboundary) 
{
    assert( pt.size() == rela.size() && 
            "Sizes of Patches and Relations do NOT MATCH!");
    radius = ( pt.front().get_bou().back().first - pt.front().get_bou().front().first ) / 2;

    Leng = pt.size();
    VEC_P = pt;
    RELA_P = rela;
    OnBOUN = onboundary;

    assert( radius > 0 && "Radius is incorrectly set!");
}



int 							PatchSet::getLeng()		{ return Leng; }
std::vector<Patch> 				PatchSet::gpt() 		{ return VEC_P; }
std::vector<std::vector<int>> 	PatchSet::gre() 		{ return RELA_P; }
std::vector< bool > 			PatchSet::getonbound()	{ return OnBOUN; }



// Take point coordinate (x,y) and number of patch as input, return phi_i(x,y)
// Do affine transformation and then pass transformed coordinates to Bump Function
std::tuple<double, double, double> PatchSet::PhiPhiD(double x, double y, int i) {
    PointC cen = VEC_P[i].get_evap().front();
    PointC cor = VEC_P[i].get_bou().front();
    double Xaff = (x - cen.first) / std::fabs( cor.first - cen.first );
    double Yaff = (y - cen.second) / std::fabs( cor.second - cen.second );

    if( std::fabs(Xaff) >= 1 - err_tol || std::fabs(Yaff) >= 1 - err_tol
        || std::pow(x-0.5,2) + std::pow(y-0.5,2) >= 0.25 + err_tol ) {
        return std::make_tuple(0.0, 0.0, 0.0);
    }

    else {
        double phisum = 0.0;
        double dxphisum = 0.0;
        double dyphisum = 0.0;

        std::tuple<double, double, double> iphi = std::make_tuple(0.0, 0.0, 0.0);

        for( auto it = RELA_P[i].begin(); it != RELA_P[i].end(); ++it ) {
            int relatedp = *it;
            PointC center = VEC_P[relatedp].get_evap().front();
            PointC corner = VEC_P[relatedp].get_bou().front();
            double Xaffined = (x - center.first) / std::fabs( corner.first - center.first );
            double Yaffined = (y - center.second) / std::fabs( corner.second - center.second );
            if( std::fabs(Xaffined) <= 1 - err_tol && std::fabs(Yaffined) <= 1 - err_tol ) {
                std::tuple<double, double, double> btmp = bump(Xaffined, Yaffined);
                phisum      += std::get<0>(btmp);
                dxphisum    += std::get<1>(btmp);
                dyphisum    += std::get<2>(btmp);
                if ( relatedp == i ) {
                    iphi = btmp;
                }
            }
        }

        //assert(   phisum >= 100*err_tol && dxphisum >= 100*err_tol &&
        //        dyphisum >= 100*err_tol && "Sums of PHI function are wrong!");
        //std::cout << "i : " << i << std::endl;
        if( phisum < 0.005) {
            std::cout << "i : " << i << std::endl
                      << "x : " << x << std::endl
                      << "y : " << y << std::endl
                      << "Xaff : " << Xaff << std::endl
                      << "Yaff : " << Yaff << std::endl
                      << "phisum : " << std::get<0>(bump(Xaff, Yaff)) << std::endl;
        }
        assert(   phisum >= 0.02 && "Sums of PHI function are wrong!");

        double phi   =   std::get<0>(iphi) / phisum;
        double dxphi = ( std::get<1>(iphi) * phisum - std::get<0>(iphi) * dxphisum ) / std::pow(phisum,2);
        double dyphi = ( std::get<2>(iphi) * phisum - std::get<0>(iphi) * dyphisum ) / std::pow(phisum,2);
        return std::make_tuple(phi, dxphi, dyphi);
    }
}


MatrixPhi PatchSet::PhiMat(double x, double y, int i) {

    MatrixPhi mp = Eigen::MatrixXd::Zero(4,3);
    std::tuple<double, double, double> iphi = PhiPhiD(x, y, i);

    PointC cen = VEC_P[i].get_evap().front();
    PointC cor = VEC_P[i].get_bou().front();
    double xsca = 1.0 / std::fabs( cor.first - cen.first );
    double ysca = 1.0 / std::fabs( cor.second - cen.second );
    double Xaff = (x - cen.first)  * xsca;
    double Yaff = (y - cen.second) * ysca;

    mp(0,0) = std::get<0>(iphi);
    mp(0,1) = std::get<1>(iphi) * xsca;
    mp(0,2) = std::get<2>(iphi) * ysca;

    mp(1,0) = Xaff * mp(0,0);
    mp(1,1) = xsca * mp(0,0) + Xaff * mp(0,1);
    mp(1,2) = Xaff * mp(0,2);

    mp(2,0) = Yaff * mp(0,0);
    mp(2,1) = Yaff * mp(0,1);
    mp(2,2) = ysca * mp(0,0) + Yaff * mp(0,2);

    mp(3,0) = Xaff * mp(2,0);
    mp(3,1) = Yaff * mp(1,1);
    mp(3,2) = Xaff * mp(2,2);

    return mp;
}




std::tuple < MatrixInt2D, MatrixInt2D, VectorInt, VectorInt> PatchSet::numericalint2d ( 
    double xbinf, double xbsup, double ybinf, double ybsup, int i, int j ) 
{
    double coeff = ( xbsup - xbinf ) * ( ybsup - ybinf ) / 4.0;
    assert( coeff > err_tol && "WRONG integral domain!");

    MatrixNd ix_grid, iy_grid, ixmg, iymg;

    ix_grid = ( xbsup - xbinf )/2 * GridXNd + ( xbsup + xbinf )/2 * OneMatNd;
    iy_grid = ( ybsup - ybinf )/2 * GridYNd + ( ybsup + ybinf )/2 * OneMatNd;

    std::pair<MatrixNd, MatrixNd> imgrid = pullback(ix_grid, iy_grid);
    ixmg = imgrid.first;
    iymg = imgrid.second;

    std::vector<MatrixNd> phi_vi = {MatrixNdZero, MatrixNdZero, MatrixNdZero, MatrixNdZero}, 
                          phi_vj = {MatrixNdZero, MatrixNdZero, MatrixNdZero, MatrixNdZero}, 
                          phi_mvj= {MatrixNdZero, MatrixNdZero, MatrixNdZero, MatrixNdZero};

    std::vector<MatrixNd> phi_grad_ix = {MatrixNdZero, MatrixNdZero, MatrixNdZero, MatrixNdZero}, 
                          phi_grad_iy = {MatrixNdZero, MatrixNdZero, MatrixNdZero, MatrixNdZero},
                          phi_grad_jx = {MatrixNdZero, MatrixNdZero, MatrixNdZero, MatrixNdZero}, 
                          phi_grad_jy = {MatrixNdZero, MatrixNdZero, MatrixNdZero, MatrixNdZero};

    for(int row = 0; row < qorder; ++row) {
        for(int col = 0; col < qorder; ++col) {
            // std::cout << "X and Y: (" << VEC_P.at(i).x_grid(row,col) << "," << VEC_P.at(i).y_grid(row,col) << ")" << std::endl; 
            MatrixPhi tpPM = PhiMat(ix_grid(row,col), iy_grid(row,col), i);
            for(int k = 0; k < NB; ++k) {
                phi_vi.at(k)(row, col) = tpPM(k,0);
                phi_grad_ix.at(k)(row, col) = tpPM(k,1);
                phi_grad_iy.at(k)(row, col) = tpPM(k,2);
            }
        }
    }

    for(int row = 0; row < qorder; ++row) {
        for(int col = 0; col < qorder; ++col) {
            // std::cout << "X and Y: (" << VEC_P.at(i).x_grid(row,col) << "," << VEC_P.at(i).y_grid(row,col) << ")" << std::endl; 
            MatrixPhi tpPM  = PhiMat(ix_grid(row,col), iy_grid(row,col), j);
            MatrixPhi mtpPM = PhiMat(ixmg(row,col), iymg(row,col), j);
            for(int k = 0; k < NB; ++k) {
                phi_mvj.at(k)(row, col) = mtpPM(k,0);
                phi_vj.at(k)(row, col) = tpPM(k,0);
                phi_grad_jx.at(k)(row, col) = tpPM(k,1);
                phi_grad_jy.at(k)(row, col) = tpPM(k,2);
            }
        }
    }

    MatrixInt2D mass, mass_map;

    for(int id = 0; id < NB; ++id) {
        for(int jd = 0; jd < NB; ++jd) {
            MatrixNd tphi_vi = phi_vi.at(id);
            MatrixNd tphi_mvj = phi_mvj.at(jd);
            MatrixNd phi_ = tphi_vi.cwiseProduct(tphi_mvj) / tstep;
            phi_ = phi_.cwiseProduct(weight2d);
            double vvv = phi_.sum() * coeff;

            mass_map(id, jd) = vvv;

            MatrixNd tphi_vj = phi_vj.at(jd);
            MatrixNd tphi_grad_ix = phi_grad_ix.at(id);
            MatrixNd tphi_grad_jx = phi_grad_jx.at(jd);
            MatrixNd tphi_grad_iy = phi_grad_iy.at(id);
            MatrixNd tphi_grad_jy = phi_grad_jy.at(jd);
            phi_ = boollap * tphi_vi.cwiseProduct(tphi_vj) / tstep 
                 + epsilon * ( tphi_grad_ix.cwiseProduct(tphi_grad_jx) 
                             + tphi_grad_iy.cwiseProduct(tphi_grad_jy) );
            phi_ = phi_.cwiseProduct(weight2d);
            vvv = phi_.sum() * coeff;
            
            mass(id, jd) = vvv;
        }
    }

    VectorInt vv0, f1;

    MatrixNd v0_value, fv1;
    for(int row = 0; row < qorder; ++row) {
        for(int col = 0; col < qorder; ++col) {
            v0_value(row,col) = initialcondition( ixmg(row,col), iymg(row,col) );
            fv1(row,col) = funf1( ix_grid(row,col), iy_grid(row,col) );
        }
    }

    for(int id = 0; id < NB; ++id) {
        MatrixNd tphi_vi = phi_vi.at(id);
        MatrixNd temp_m = v0_value.cwiseProduct(tphi_vi) / tstep;
        temp_m = temp_m.cwiseProduct(weight2d);
        vv0(id)= temp_m.sum() * coeff;
        
        temp_m = fv1.cwiseProduct(tphi_vi);
        temp_m = temp_m.cwiseProduct(weight2d);
        f1(id) = temp_m.sum() * coeff;
    }

    return std::make_tuple(mass, mass_map, vv0, f1);
}


std::tuple < MatrixInt2D, MatrixInt2D, VectorInt, VectorInt> PatchSet::boundary_numericalint2d ( 
    double xbinf, double xbsup, double ybinf, double ybsup, int i, int j ) {
    if ( patchindomain( xbinf, xbsup, ybinf, ybsup ) == 4 ) 
        return numericalint2d ( xbinf, xbsup, ybinf, ybsup, i, j );
    else if( std::abs(xbsup - xbinf) <= patch_radius*patch_radius || patchindomain( xbinf, xbsup, ybinf, ybsup ) == 0 )
        return zero_tuple;
    else {
        double xbmid = ( xbinf + xbsup ) / 2;
        double ybmid = ( ybinf + ybsup ) / 2;
        std::tuple < MatrixInt2D, MatrixInt2D, VectorInt, VectorInt>
            term1 = boundary_numericalint2d( xbinf, xbmid, ybinf, ybmid, i, j );
        std::tuple < MatrixInt2D, MatrixInt2D, VectorInt, VectorInt>
            term2 = boundary_numericalint2d( xbinf, xbmid, ybmid, ybsup, i, j );
        std::tuple < MatrixInt2D, MatrixInt2D, VectorInt, VectorInt>
            term3 = boundary_numericalint2d( xbmid, xbsup, ybinf, ybmid, i, j );
        std::tuple < MatrixInt2D, MatrixInt2D, VectorInt, VectorInt>
            term4 = boundary_numericalint2d( xbmid, xbsup, ybmid, ybsup, i, j );
        return std::make_tuple( std::get<0>(term1) + std::get<0>(term2) + std::get<0>(term3) + std::get<0>(term4) ,
                                std::get<1>(term1) + std::get<1>(term2) + std::get<1>(term3) + std::get<1>(term4) ,
                                std::get<2>(term1) + std::get<2>(term2) + std::get<2>(term3) + std::get<2>(term4) ,
                                std::get<3>(term1) + std::get<3>(term2) + std::get<3>(term3) + std::get<3>(term4) );
    }
}


MatrixInt2D PatchSet::numericalint1d ( int i, int j ) 
{
    Vector2Nd x_line_tmp, y_line_tmp;
    //x_line_tmp = VEC_P.at(i).x_line;
    //y_line_tmp = VEC_P.at(i).y_line;
    x_line_tmp = VEC_P[i].get_xline();
    y_line_tmp = VEC_P[i].get_yline();
    //
    std::vector<Vector2Nd> phib_vi = {Vector2NdZero, Vector2NdZero, Vector2NdZero, Vector2NdZero}, 
                           phib_vj = {Vector2NdZero, Vector2NdZero, Vector2NdZero, Vector2NdZero};
    std::vector<Vector2Nd> phib_grad_ix = {Vector2NdZero, Vector2NdZero, Vector2NdZero, Vector2NdZero}, 
                           phib_grad_iy = {Vector2NdZero, Vector2NdZero, Vector2NdZero, Vector2NdZero}, 
                           phib_grad_jx = {Vector2NdZero, Vector2NdZero, Vector2NdZero, Vector2NdZero}, 
                           phib_grad_jy = {Vector2NdZero, Vector2NdZero, Vector2NdZero, Vector2NdZero};

    for(int row = 0; row < 2*qorder; ++row) {
        MatrixPhi tpPM = PhiMat(x_line_tmp(row), y_line_tmp(row), i);
        for(int k = 0; k < NB; ++k) {
            phib_vi.at(k)(row) = tpPM(k,0);
            phib_grad_ix.at(k)(row) = tpPM(k,1);
            phib_grad_iy.at(k)(row) = tpPM(k,2);
        }
    }
    
    for(int row = 0; row < 2*qorder; ++row) {
        MatrixPhi tpPM = PhiMat(x_line_tmp(row), y_line_tmp(row), j);
        for(int k = 0; k < NB; ++k) {
            phib_vj.at(k)(row) = tpPM(k,0);
            phib_grad_jx.at(k)(row) = tpPM(k,1);
            phib_grad_jy.at(k)(row) = tpPM(k,2);
        }
    }

    MatrixInt2D mass;

    for(int id = 0; id < NB; ++id) {
        Vector2Nd phib_grad_i = phib_grad_ix.at(id).cwiseProduct(VEC_P.at(i).get_xnormv()) 
                              + phib_grad_iy.at(id).cwiseProduct(VEC_P.at(i).get_ynormv());
        Vector2Nd tphib_vi = phib_vi.at(id);
        
        for(int jd = 0; jd < NB; ++jd) {

            Vector2Nd tphib_vj = phib_vj.at(jd);
            Vector2Nd phib_grad_j = phib_grad_jx.at(jd).cwiseProduct(VEC_P.at(i).get_xnormv()) 
                                  + phib_grad_jy.at(jd).cwiseProduct(VEC_P.at(i).get_ynormv());
            Vector2Nd phib_grad   = - epsilon * ( tphib_vi.cwiseProduct(phib_grad_j) 
                                                + tphib_vj.cwiseProduct(phib_grad_i) ) 
                                    + beta * tphib_vi.cwiseProduct(tphib_vj);
                      phib_grad   = phib_grad.cwiseProduct(WW1d);

            mass(id, jd) = phib_grad.sum() * VEC_P.at(i).diam_line;
        }
    }

    return mass;
}




void PatchSet::QEVAL_uv(Eigen::SparseMatrix<double> &M, 
						Eigen::SparseMatrix<double> &M_mapped, 
						Eigen::VectorXd &rhs, 
                        Eigen::VectorXd &frhs) {

    assert( rhs.size() == 4 * Leng && 
                 "Size of RHS vector input is not correct!");

    std::vector< Eigen::Triplet<double> > tripletList;
    std::vector< Eigen::Triplet<double> > tripletList_mapped;
    tripletList.reserve( Leng * 16 * 9 );
    tripletList_mapped.reserve ( Leng * 16 * 9 );

    for( int i = 0; i < Leng; i++ ) {
        for( auto it = RELA_P[i].begin(); it != RELA_P[i].end(); ++it ) {
            int j = *it;

            // Compute integral \int u*v on patch i WHEN patch i is on boundary
            if ( OnBOUN.at(i) && OnBOUN.at(j) )
            {
            	double xbinf, xbsup, ybinf, ybsup;
            	// xbsup = std::min(VEC_P[i].get_bou().back().first,	xsup);
            	// xbinf = std::max(VEC_P[i].get_bou().front().first,	xinf);
            	// ybsup = std::min(VEC_P[i].get_bou().back().second,	ysup);
            	// ybinf = std::max(VEC_P[i].get_bou().front().second,	yinf);

                xbsup = VEC_P[i].get_bou().back().first;
                xbinf = VEC_P[i].get_bou().front().first;
                ybsup = VEC_P[i].get_bou().back().second;
                ybinf = VEC_P[i].get_bou().front().second;

                std::tuple < MatrixInt2D, MatrixInt2D, VectorInt, VectorInt> 
                    item = boundary_numericalint2d ( xbinf, xbsup, ybinf, ybsup, i, j );

                MatrixInt2D mass1d_local = numericalint1d (i, j);

                MatrixInt2D mass_local      = std::get<0>(item) + mass1d_local;
                /*
                if(i == j) {
                    for(int inde = 0; inde<NB; inde++) {
                        std::cout << mass_local(inde,inde)-mass1d_local(inde,inde) << " , " << mass1d_local(inde,inde) << std::endl;
                    } 
                }
                */
                MatrixInt2D mass_map_local  = std::get<1>(item);
                rhs.segment(4*i, 4)         = std::get<2>(item);
                frhs.segment(4*i, 4)       = std::get<3>(item);


                for(int id = 0; id < NB; ++id) {
                    for(int jd = 0; jd < NB; ++jd) {
                        tripletList_mapped.push_back(Eigen::Triplet<double>(4*i+id, 4*j+jd, mass_map_local(id, jd)));
                        tripletList.push_back(Eigen::Triplet<double>(4*i+id, 4*j+jd, mass_local(id, jd)));
                    }
                }
            }

            // Compute integral \int u*v on patch i WHEN patch i is on boundary
            else 
            {
                std::tuple < MatrixInt2D, MatrixInt2D, VectorInt, VectorInt > 
                    item = numericalint2d ( VEC_P[i].get_bou().front().first, VEC_P[i].get_bou().back().first,
                                            VEC_P[i].get_bou().front().second, VEC_P[i].get_bou().back().second, 
                                            i, j );

                MatrixInt2D mass_local      = std::get<0>(item);
                MatrixInt2D mass_map_local  = std::get<1>(item);
                rhs.segment(4*i, 4)         = std::get<2>(item);
                frhs.segment(4*i, 4)       = std::get<3>(item);
                /*
                if(i == j) {
                    std::cout << mass_local.diagonal() << std::endl;
                }
                */
                for(int id = 0; id < NB; ++id) {
                    for(int jd = 0; jd < NB; ++jd) {
                        tripletList_mapped.push_back(Eigen::Triplet<double>(4*i+id, 4*j+jd, mass_map_local(id, jd)));
                        tripletList.push_back(Eigen::Triplet<double>(4*i+id, 4*j+jd, mass_local(id, jd)));
                    }
                }
            }
        }
    }

    M.setFromTriplets(tripletList.begin(), tripletList.end());
    M.makeCompressed();
    M_mapped.setFromTriplets(tripletList_mapped.begin(), tripletList_mapped.end());
    M_mapped.makeCompressed();
}

void PatchSet::evaluate( Eigen::VectorXd c, int Ntest, Eigen::MatrixXd &r ) {
    assert( c.size() == 4*Leng &&
            "Sizes of Coefficients and Patchset do NOT MATCH!");

    double rtmp = 1.0/(Ntest-1);
    PointC coor;
    r = Eigen::MatrixXd::Zero(Ntest, Ntest);

    for(int j = 0; j < Ntest; j++) {
        for(int k = 0; k < Ntest; k++) {
            coor.first = (double)j * rtmp;
            coor.second = (double)k * rtmp;
            for( int i = 0; i < Leng; i++ ) {
                if( VEC_P[i].belong(coor) && pointindomain(coor.first, coor.second) ) {
                    MatrixPhi mp = PhiMat(coor.first, coor.second, i);
                    r(j,k) += c.segment(4*i, 4).transpose() * mp.col(0);
                }   
            }
        }
    }
}


void PatchSet::PS_print() {
    for(int i = 0; i < Leng; ++i) {
        VEC_P[i].print_info();
        std::cout << "Relations:" << std::endl;
            for(auto it = RELA_P[i].begin(); it != RELA_P[i].end(); ++it) {
                std::cout << *it << ' ';
            }
        std::cout << std::endl
                << "On Boundary information:" << std::endl 
                << OnBOUN[i] << std::endl;
    }

    std::cout << std::endl;

    std::cout << "Number of Patches: " << Leng << std::endl << std::endl
              << "Radius of Patches: " << radius << std::endl << std::endl;

    std::cout << std::endl;
}



# endif