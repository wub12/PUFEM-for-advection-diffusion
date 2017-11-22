# ifndef GENE_CIRC_HPP
# define GENE_CIRC_HPP

#include <ctime>
#include "rectdomain.hpp"

// Here we generate a set of Patches for Rectangular Domain [0,1]X[0,1]
// Number of Patches is np*np
// Freedom of each patch is 1 && Evaluation point is center of Patch

std::pair<double, double> circ_intersect(double x1, double x2, double y1, double y2);

PatchSet circ_generator(int np) {
    assert( np >= 2 && "np must be LARGER than 1!!!");
    // If np = 2, check boundary integral in PatchSet::QEVAL_ugv in RECTDOMAIN_HPP

    std::vector<Patch> res;
    std::vector< std::vector<int> > relation;
    std::vector< bool > onboundary;
    res.reserve(np*np);
    relation.reserve(np*np);

    double rad;
    // Notice that 1.0 to make sure decimal is created
    rad = 1.0 / (np-1);
    double radwid = rad * scale;
    //std::cout << "rad: " << rad << std::endl;

    for( int i=0; i<np; ++i ) {
        for( int j=0; j<np; ++j) {
            // Generate coordinates of Points on Boundary and Evaluation Points
            int numin = patchindomain( i*rad-radwid, i*rad+radwid, j*rad-radwid, j*rad+radwid);
            if ( numin >= 1 &&  std::sqrt( std::pow(i*rad-0.5,2) + std::pow(j*rad-0.5,2) ) < 0.5+radwid*0.4 ) {
                double dist = std::sqrt( std::pow(i*rad-0.5,2) + std::pow(j*rad-0.5,2) );
                double modi = 1.0;
                /*
                if ( numin <= 3 )
                    std::cout << dist << std::endl;
                    
                if ( dist > 0.5 && numin <= 3 )
                {
                    modi = 0.5/dist;
                    // std::cout << modi << std::endl;
                }
                */
                std::pair<double, double> cen( 0.5 + modi*(i*rad-0.5), 0.5 + modi*(j*rad-0.5) );
                std::pair<double, double> p1( cen.first - radwid, cen.second - radwid ), 
                                          p2( cen.first - radwid, cen.second + radwid ), 
                                          p3( cen.first + radwid, cen.second - radwid ), 
                                          p4( cen.first + radwid, cen.second + radwid );
                Vec boutemp {p1, p2, p3, p4};
                Vec evaptemp {cen};
    
                // Set ON BOUNDARY information
                bool onbtemp;
                if( numin <= 3 ) {
                    onbtemp = true;
                    onboundary.push_back(onbtemp);

                    Vector2Nd xline, yline, xnor, ynor;
                    std::pair<double, double> arcv = circ_intersect( cen.first - radwid,  cen.first + radwid, 
                                                                     cen.second - radwid, cen.second + radwid );

                    Vector2Nd arcnode;
                    arcnode.segment(0, qorder)      = Node1d * (arcv.second-arcv.first)/2 + OneVecNd * (arcv.second+arcv.first)/2;
                    arcnode.segment(qorder, qorder) = arcnode.segment(0, qorder);

                    for(int i=0; i<2*qorder; i++) {
                        xnor(i)  = std::cos(arcnode(i));
                        xline(i) = 0.5 + 0.5 * xnor(i) ;
                        ynor(i)  = std::sin(arcnode(i));
                        yline(i) = 0.5 + 0.5 * ynor(i) ;
                    }

                    double arclength = (arcv.second - arcv.first) * 0.5;
                    Patch pattemp(boutemp, evaptemp, std::make_pair(xline,yline), xnor, ynor, arclength);
                    res.push_back(pattemp);
                }
                else {
                    onbtemp = false;
                    onboundary.push_back(onbtemp);
                    std::pair<Vector2Nd, Vector2Nd> tmp(ZeroVec2Nd, ZeroVec2Nd);
                    Patch pattemp(boutemp, evaptemp, tmp, ZeroVec2Nd, ZeroVec2Nd, 0.0);
                    res.push_back(pattemp);
                }
            }
        }
    }

    // Build Relationship Vector
    for( int pt=0; pt<res.size(); pt++ ) {
        std::vector<int> overlap;
        std::pair<double, double> aimp = res.at(pt).get_evap().front();
        for ( int k=0; k<res.size(); k++ ) {
            std::pair<double, double> intersp = res.at(k).get_evap().front();
            if( std::abs(aimp.first-intersp.first) < rad*1.4 && std::abs(aimp.second-intersp.second) < rad*1.4 )
                overlap.push_back(k);
        }
        relation.push_back(overlap);
    }
    
    std::time_t t = time(0);
    struct tm * now = localtime( & t );
    std::cout << "Components of patchset information is prepared!" << std::endl
              << "Done at " << asctime(now) << ">>>>>>>>>>>>>>>>>>>>" << std::endl;

    PatchSet infopatch ( res, relation, onboundary );

    t = time(0);
    now = localtime( & t );
    std::cout << "Patchset is built!" << std::endl
              << "Done at " << asctime(now) << ">>>>>>>>>>>>>>>>>>>>" << std::endl;
    return infopatch;
}


std::pair<double, double> circ_intersect(double x1, double x2, double y1, double y2) {
    // Inputs should be "(x1,y1) left -lower point"
    //                  "(x2,y2) right-upper point"
    assert( x1 < x2 && y1 < y2 && "Input is wrong!" );
    double PI = 3.1415926;
    double xtemp, xtemp_, ytemp, ytemp_;
    Eigen::VectorXd set(16);
    int i = 0;
    if( x1 < xsup-err_tol && x1 > xinf+err_tol ) {
        ytemp = std::sqrt( 0.25 - std::pow(x1-0.5,2) ) + 0.5;
        ytemp_ = -std::sqrt( 0.25 - std::pow(x1-0.5,2) ) + 0.5;
        if( y1 < ytemp-err_tol && y2 > ytemp+err_tol ) {
            set(i++) = x1; set(i++) = ytemp;
        }
        else if( y1 < ytemp_-err_tol && y2 > ytemp_+err_tol ) {
            set(i++) = x1; set(i++) = ytemp_;
        }
    }

    if( x2 < xsup-err_tol && x2 >= xinf+err_tol ) {
        ytemp = std::sqrt( 0.25 - std::pow(x2-0.5,2) ) + 0.5;
        ytemp_ = -std::sqrt( 0.25 - std::pow(x2-0.5,2) ) + 0.5;
        if( y1 < ytemp-err_tol && y2 > ytemp+err_tol ) {
            set(i++) = x2; set(i++) = ytemp;
        }
        else if( y1 < ytemp_-err_tol && y2 > ytemp_+err_tol ) {
            set(i++) = x2; set(i++) = ytemp_;
        }
    }

    if( y1 < ysup-err_tol && y1 > yinf+err_tol ) {
        xtemp = std::sqrt( 0.25 - std::pow(y1-0.5,2) ) + 0.5;
        xtemp_ = -std::sqrt( 0.25 - std::pow(y1-0.5,2) ) + 0.5;
        if( x1 < xtemp-err_tol && x2 > xtemp+err_tol ) {
            set(i++) = xtemp; set(i++) = y1;
        }
        else if( x1 < xtemp_-err_tol && x2 > xtemp_+err_tol ) {
            set(i++) = xtemp_; set(i++) = y1;
        }
    }

    if( y2 < ysup-err_tol && y2 > yinf+err_tol ) {
        xtemp = std::sqrt( 0.25 - std::pow(y2-0.5,2) ) + 0.5;
        xtemp_ = -std::sqrt( 0.25 - std::pow(y2-0.5,2) ) + 0.5;
        if( x1 < xtemp-err_tol && x2 > xtemp+err_tol ) {
            set(i++) = xtemp; set(i++) = y2;
        }
        else if( x1 < xtemp_-err_tol && x2 > xtemp_+err_tol ) {
            set(i++) = xtemp_; set(i++) = y2;
        }
    }

    if(i < 4){
        std::cout << x1 << std::endl << x2 <<
        std::endl << y1 << std::endl << y2 << std::endl;
    }
    assert( i>=4 && "Input is wrong!" );


    double t1 = std::atan( (set(1)-0.5) / (set(0)-0.5) );
    double t2 = std::atan( (set(3)-0.5) / (set(2)-0.5) );
    /*
    double  t1 = std::min( tt1, tt2 );
    double  t2 = std::max( tt1, tt2 );
    
    for(int k=2; k<i/2; k++) {
        double ttt1 = std::atan( (set(2*k+1)-0.5) / (set(2*k)-0.5) );
        t1 = std::min( t1, ttt1 );
        t2 = std
    }
    */

    if( set(0)-0.5<0 )
        t1 += PI;
    if( set(2)-0.5<0 )
        t2 += PI;

    if( t1-t2 > PI )
        t2 = t2+2*PI;
    if( t2-t1 > PI )
        t1 = t1+2*PI;

    return std::make_pair(std::min(t1,t2), std::max(t1,t2));
}

#endif