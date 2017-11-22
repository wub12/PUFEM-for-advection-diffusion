# ifndef GENE_RECT_HPP
# define GENE_RECT_HPP

#include <ctime>
#include "rectdomain.hpp"

// Here we generate a set of Patches for Rectangular Domain [0,1]X[0,1]
// Number of Patches is np*np
// Freedom of each patch is 1 && Evaluation point is center of Patch



PatchSet rect_generator(int np) {
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
            
            std::pair<double, double> cen( i*rad, j*rad );
            std::pair<double, double> p1( cen.first - radwid, cen.second - radwid ), 
                                      p2( cen.first - radwid, cen.second + radwid ), 
                                      p3( cen.first + radwid, cen.second - radwid ), 
                                      p4( cen.first + radwid, cen.second + radwid );

            Vec boutemp {p1, p2, p3, p4};
            Vec evaptemp {cen};


            // Set ON BOUNDARY information
            bool onbtemp;
            if( i==0 || j==0 || i==np-1 || j==np-1 ) {
                onbtemp = true;
                VectorNd xline, yline, xnor, ynor;

                if( i == 0 && j == 0 ) {
                    Vector2Nd x_orig, y_orig;
                    x_orig = -1.0 * LineNode1;
                    y_orig = LineNode2;
                    std::pair<Vector2Nd, Vector2Nd> tmp = coortrans_v( x_orig, y_orig, cen.first, cen.second, p4.first, p4.second );
                    Patch pattemp(boutemp, evaptemp, tmp, -1.0 * ZeroOneVec2Nd, -1.0 * OneZeroVec2Nd, 2 * radwid);
                    res.push_back(pattemp);
                }
                else if( i == 0 && j == np-1 ) {
                    Vector2Nd x_orig, y_orig;
                    x_orig = -1.0 * LineNode1;
                    y_orig = -1.0 * LineNode2;
                    std::pair<Vector2Nd, Vector2Nd> tmp = coortrans_v( x_orig, y_orig, cen.first, cen.second, p4.first, p4.second );
                    Patch pattemp(boutemp, evaptemp, tmp, -1.0 * ZeroOneVec2Nd, OneZeroVec2Nd, 2 * radwid);
                    res.push_back(pattemp);
                }
                else if( i == np-1 && j == 0 ) {
                    Vector2Nd x_orig, y_orig;
                    x_orig = LineNode1;
                    y_orig = LineNode2;
                    std::pair<Vector2Nd, Vector2Nd> tmp = coortrans_v( x_orig, y_orig, cen.first, cen.second, p4.first, p4.second );
                    Patch pattemp(boutemp, evaptemp, tmp, ZeroOneVec2Nd, -1.0 * OneZeroVec2Nd, 2 * radwid);
                    res.push_back(pattemp);
                }
                else if( i == np-1 && j == np-1 ) {
                    Vector2Nd x_orig, y_orig;
                    x_orig = LineNode1;
                    y_orig = -1.0 * LineNode2;
                    std::pair<Vector2Nd, Vector2Nd> tmp = coortrans_v( x_orig, y_orig, cen.first, cen.second, p4.first, p4.second );
                    Patch pattemp(boutemp, evaptemp, tmp, ZeroOneVec2Nd, OneZeroVec2Nd, 2 * radwid);
                    res.push_back(pattemp);
                }
                else if( i == 0 || i == np-1 ) {
                    std::pair<Vector2Nd, Vector2Nd> tmp = coortrans_v( ZeroVec2Nd, NN1d, cen.first, cen.second, p4.first, p4.second );
                    if( i == 0 ) {
                        Patch pattemp(boutemp, evaptemp, tmp, -1. * OneVec2Nd, ZeroVec2Nd, 2 * radwid);
                        res.push_back(pattemp);
                    }
                    else {
                        Patch pattemp(boutemp, evaptemp, tmp, OneVec2Nd, ZeroVec2Nd, 2 * radwid);
                        res.push_back(pattemp);
                    }
                }
                else if( j == 0 || j == np-1  ) {
                    std::pair<Vector2Nd, Vector2Nd> tmp = coortrans_v( NN1d, ZeroVec2Nd, cen.first, cen.second, p4.first, p4.second );
                    if( j == 0 ) {
                        Patch pattemp(boutemp, evaptemp, tmp, ZeroVec2Nd, -1. * OneVec2Nd, 2 * radwid);
                        res.push_back(pattemp);
                    }
                    else {
                        Patch pattemp(boutemp, evaptemp, tmp, ZeroVec2Nd, OneVec2Nd, 2 * radwid);
                        res.push_back(pattemp);
                    }
                }
                else
                    assert( 1 == 2 && "$local$ variable is wrong when compute Phi" );
            }
            else {
                onbtemp = false;
                std::pair<Vector2Nd, Vector2Nd> tmp(ZeroVec2Nd, ZeroVec2Nd);
                Patch pattemp(boutemp, evaptemp, tmp, ZeroVec2Nd, ZeroVec2Nd, 0.0);
                res.push_back(pattemp);
            }

            onboundary.push_back(onbtemp);



            // Build Relationship Vector
            std::vector<int> overlap;
            overlap.insert( overlap.end(), { (i-1)*np + j, i*np + j-1, i*np + j+1, (i+1)*np + j, 
                                                                             i*np + j,
                                                                             (i-1)*np + j-1, (i-1)*np + j+1, (i+1)*np + j-1, (i+1)*np + j+1 } );
            if(i==0 && j==0) {
                overlap.erase(overlap.begin()+5, overlap.begin()+8);
                overlap.erase(overlap.begin(), overlap.begin()+2);
            }
            else if(i==0 && j==np-1) {
                overlap.erase(overlap.begin()+8);
                overlap.erase(overlap.begin()+5, overlap.begin()+7);
                overlap.erase(overlap.begin()+2);
                overlap.erase(overlap.begin());
            }
            else if(i==np-1 && j==0) {
                overlap.resize(7);
                overlap.erase(overlap.begin()+5);
                overlap.erase(overlap.begin()+3);
                overlap.erase(overlap.begin()+1);
            }
            else if(i==np-1 && j==np-1) {
                overlap.resize(6);
                overlap.erase(overlap.begin()+3);
                overlap.erase(overlap.begin()+2);
            }
            else if(i==0) {
                overlap.erase(overlap.begin()+5, overlap.begin()+7);
                overlap.erase(overlap.begin());
            }
            else if(i==np-1) {
                overlap.resize(7);
                overlap.erase(overlap.begin()+3);
            }
            else if(j==0) {
                overlap.erase(overlap.begin()+7);
                overlap.erase(overlap.begin()+5);
                overlap.erase(overlap.begin()+1);
            }
            else if(j==np-1) {
                overlap.erase(overlap.begin()+8);
                overlap.erase(overlap.begin()+6);
                overlap.erase(overlap.begin()+2);
            }
            // Push the overlap relation to relation vector
            relation.push_back(overlap);
        }
    }
    
    std::time_t t = time(0);
    struct tm * now = localtime( & t );
    std::cout << "Components of patchset information is prepared!" << std::endl
              << "Done at " << asctime(now) << ">>>>>>>>>>>>>>>>>>>>" << std::endl;

    PatchSet infopatch ( res, relation, onboundary);
    t = time(0);
    now = localtime( & t );
    std::cout << "Patchset is built!" << std::endl
              << "Done at " << asctime(now) << ">>>>>>>>>>>>>>>>>>>>" << std::endl;
    return infopatch;
}

#endif