#include <iostream>
#include <fstream>
#include <iomanip>
#include <ctime>

// g++ main.cpp -I./../../eigenlib -o main
// g++ main.cpp -O3 -funroll-loops -march=native -mfpmath=sse -o main_enhance

#include <Eigen/Dense>
#include <Eigen/SparseLU>
#include <Eigen/SparseCholesky>
#include "quadeval.hpp"
#include "gene_circ.hpp"






int main()
{
    // Number of Patches is NPatch*NPatch
    // Width of Patches is 1/(NPatch-1)

    std::time_t t = time(0);
    struct tm * now = localtime( & t );
    std::cout << "Started!" << std::endl
              << "Done at " << asctime(now) << ">>>>>>>>>>>>>>>>>>>>" << std::endl;

    // Generate PATCHSET
    PatchSet re = circ_generator(NPatch);
    int np = re.getLeng();
    
    std::cout << "NUMBER of patches:" << np << std::endl;
    //re.PS_print();

    
    t = time(0);
    now = localtime( & t );
    std::cout << "PatchSet generated!" << std::endl
              << "Done at " << asctime(now) << ">>>>>>>>>>>>>>>>>>>>" << std::endl;

    // Compute matrix for the bilinear form A
    Eigen::VectorXd v = Eigen::VectorXd::Zero(np*4);
    Eigen::VectorXd w1 = Eigen::VectorXd::Zero(np*4);
    Eigen::VectorXd w2 = Eigen::VectorXd::Zero(np*4);
    Eigen::VectorXd sol = Eigen::VectorXd::Zero(np*4);

    Eigen::SparseMatrix<double> M(4*np, 4*np), M_mapped(4*np, 4*np);
    re.QEVAL_uv( M, M_mapped, v, w1);

    Eigen::SparseMatrix<double> A, B;
    A = (M + Eigen::SparseMatrix<double>(M.transpose()))/2;
    // A = M;

    t = time(0);
    now = localtime( & t );
    std::cout << "Matrix computed!" << std::endl
              << "Done at " << asctime(now) << ">>>>>>>>>>>>>>>>>>>>" << std::endl;

    //std::cout << A << std::endl;
    //std::cout << M_mapped << std::endl;

    // Eigen::SimplicialLDLT< Eigen::SparseMatrix<double> > solver;
    Eigen::SparseLU< Eigen::SparseMatrix<double> > solver;
    solver.analyzePattern(A);
    solver.factorize(A);
    t = time(0);
    now = localtime( & t );
    std::cout << "LU decomposition of matrix computed!" << std::endl
              << "Done at " << asctime(now) << ">>>>>>>>>>>>>>>>>>>>" << std::endl;

    // Solve the first step
    sol = solver.solve( v + w1 );
    //std::cout << sol << std::endl;
    
    for (int i = 1; i < iteration; ++i)
    {
        double tt = i * tstep;
        v = M_mapped * sol + w1 * exp(tt);
        sol = solver.solve(v);
        if(i % 10000 == 0) {
            t = time(0);
            now = localtime( & t );
            std::cout << i <<" st step have finished" << std::endl
                      << "Done at " << asctime(now) << ">>>>>>>>>>>>>>>>>>>>" << std::endl;
        }
    }
    

    int Ntest = 320;
    Eigen::MatrixXd sol_g(Ntest, Ntest), err_g(Ntest, Ntest);
    double l2err = 0;
    double l2norm = 0;
    re.evaluate(sol, Ntest, sol_g);
    
    double rtmp = 1.0/(Ntest-1);
    for(int i = 0; i < Ntest; ++i) {
        for(int j = 0; j < Ntest; ++j) {
            double xc = (double)i * rtmp;
            double yc = (double)j * rtmp;
            double ts = theo_solution(xc,yc);
            err_g(i,j) = sol_g(i,j) - ts;   
            //err_g(i,j) = sol_g(i,j) - initialcondition(xc,1-yc);
            if( pointindomain(xc, yc) ) {
                l2err += std::pow( err_g(i,j), 2 );
                l2norm += std::pow( ts, 2 );
            }
            //l2norm += std::pow( initialcondition(xc,1-yc), 2 );
        }
    }


    l2err = std::sqrt( l2err / (double)(Ntest*Ntest) );
    l2norm = std::sqrt( l2norm / (double)(Ntest*Ntest) );
    std::cout << "Numb. of Patch = " << NPatch << " * " << NPatch << std::endl;
    std::cout << "Orde. of Quadr = " << qorder << std::endl;
    std::cout << "Beta  coeff.   = " << beta << std::endl;
    std::cout << "Overlap coeff  = " << scale << std::endl;
    std::cout << "Abso. L2_error = " << l2err << std::endl;
    std::cout << "Abso. L2_norm  = " << l2norm << std::endl;
    std::cout << "Rela. L2_error = " << l2err / l2norm << std::endl;
    

    std::ofstream file_res, file_err;
    file_res.open("result.txt");
    file_err.open("error.txt");
    for(int i=0; i<Ntest; ++i) {
        for(int j=0; j<Ntest; ++j) {
            file_res << std::setprecision(8) << sol_g(j, Ntest-i-1);
            file_res << ",";
            file_err << std::setprecision(8) << err_g(j, Ntest-i-1);
            file_err << ",";
        }
        file_res << "\n";
        file_err << "\n";
    }

    


    file_res.close();
    file_err.close();

    
    return 0;
}