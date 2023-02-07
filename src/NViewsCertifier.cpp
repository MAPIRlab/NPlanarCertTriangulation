#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <sstream>

// include types
#include "NViewsTypes.h"

// include utils
#include "NViewsUtils.h"

// include headers
#include "NViewsCertifier.h"

// for eigendecomposition
// #include <eigen3/Eigen/Dense>
#include <Eigen/Eigenvalues> 

#include <chrono>  // timer

using namespace std::chrono;



namespace NPlanarTrian
{



double NViewCertClass::computeMult(const Eigen::VectorXd & sol, 
                                   Eigen::VectorXd & sol_mult)
{
        double f_sol = sol.dot(sol); 
       
        // this one is already good
        Eigen::VectorXd b(M_ + 1); 
        b << sol, -f_sol;              
 
        Eigen::MatrixXd A; 
        
        A.resize( (int) (M_ + 1), (int) (constr_red_1_.size() * 2)); 
    
        A.setZero(); 
        
        auto start_t_mult = high_resolution_clock::now();
        for(int i=0; i < constr_red_1_.size(); i++)
        {
                Vector2 p1 = sol.block<2,1>(constr_red_1_[i].id_1, 0); 
                Vector2 p2 = sol.block<2,1>(constr_red_1_[i].id_2, 0); 
               
               
                // Constraint #1
                A.block<2,1>(constr_red_1_[i].id_1, i * 2) = 0.5 * (constr_red_1_[i].H * p2 + constr_red_1_[i].Hp2);
                A.block<2,1>(constr_red_1_[i].id_2, i * 2) = 0.5 * (constr_red_1_[i].H.transpose() * p1 + constr_red_1_[i].Hp1);                
                
                A(sol.size(), i * 2) = constr_red_1_[i].b + 0.5 * (constr_red_1_[i].Hp2.dot(p1) + constr_red_1_[i].Hp1.dot(p2));
                
                // Constraint #2
                A.block<2,1>(constr_red_1_[i].id_1, i * 2 + 1) = 0.5 * (constr_red_2_[i].H * p2 + constr_red_2_[i].Hp2);
                A.block<2,1>(constr_red_1_[i].id_2, i * 2 + 1) = 0.5 * (constr_red_2_[i].H.transpose() * p1 + constr_red_2_[i].Hp1);                
                
                A(sol.size(), i * 2 + 1) = constr_red_2_[i].b + 0.5 * (constr_red_2_[i].Hp2.dot(p1) + constr_red_2_[i].Hp1.dot(p2));
                
        }
        
       
        auto time_mult = duration_cast<nanoseconds>(high_resolution_clock::now() - start_t_mult);
        // std::cout << "[MULT2] A:\n" << A << std::endl;
        // A is sparse
        // solve system A*l = b
        
              //   std::cout << "Size A ref: " << A.rows() << "," << A.cols() << std::endl; 
        // Eigen::JacobiSVD<Eigen::MatrixXd> svd(A);
        // std::cout << "[MULT] Singular values A:\n" << svd.singularValues() << std::endl;
        // std::cout << "[MULT] Matrix constr: " << (double) time_mult.count() << std::endl; 
        double error_sol = 1;
        auto start_t_mult2 = high_resolution_clock::now();
        if (A.rows() >= A.cols())
                error_sol = solveLinearSystemMinNorm(A.transpose() * A, A.transpose()*b, sol_mult);
        else
                error_sol = solveLinearSystemMinNorm(A, b, sol_mult);
        
        auto time_mult2 = duration_cast<nanoseconds>(high_resolution_clock::now() - start_t_mult2);
       
        
        // std::cout << "[MULT] Matrix solve: " << (double) time_mult2.count() << std::endl; 
        // std::cout << "[MULT] Matrix: " << (double) time_mult.count() << std::endl;
        // std::cout << "[MULT] Min norm: " << (double) time_mult2.count() << std::endl;
        
        
        return error_sol;
}




double NViewCertClass::computeHessian(const double cost_sol, 
                                      const Eigen::VectorXd & mult, 
                                      Eigen::MatrixXd & Hess)                     
{
        auto start_t_mult = high_resolution_clock::now();
        Hess = Eigen::MatrixXd::Identity(M_+1, M_+1); 
        Hess(M_, M_) = - cost_sol;
        int N = constr_red_1_.size();
        
        for (int i=0; i <N; i++)
        {
                int id_1 = constr_red_1_[i].id_1; 
                int id_2 = constr_red_1_[i].id_2;
                
                // Constraint #1
                Hess(M_, M_)                   -=       mult(i * 2) * constr_red_1_[i].b;   
                
                Hess.block<2,2>(id_1, id_2)    -= 0.5 * mult(i * 2) * constr_red_1_[i].H; 
                // the symmetric part
                Hess.block<2,2>(id_2, id_1)    -= 0.5 * mult(i * 2) * constr_red_1_[i].H.transpose();  
                
                Hess.block<2,1>(id_1, M_)      -= 0.5 * mult(i * 2) * constr_red_1_[i].Hp2; 
                // symmetry
                Hess.block<1,2>(M_, id_1)      -= 0.5 * mult(i * 2) * constr_red_1_[i].Hp2.transpose();  
                Hess.block<2,1>(id_2, M_)      -= 0.5 * mult(i * 2) * constr_red_1_[i].Hp1; 
                // symmetry
                Hess.block<1,2>(M_, id_2)      -= 0.5 * mult(i * 2) * constr_red_1_[i].Hp1.transpose();   
                
                
                // Constraint #2
                Hess(M_, M_)                   -=       mult(i * 2 + 1) * constr_red_2_[i].b;   
                
                Hess.block<2,2>(id_1, id_2)    -= 0.5 * mult(i * 2 + 1) * constr_red_2_[i].H; 
                // the symmetric part
                Hess.block<2,2>(id_2, id_1)    -= 0.5 * mult(i * 2 + 1) * constr_red_2_[i].H.transpose();  
                
                Hess.block<2,1>(id_1, M_)      -= 0.5 * mult(i * 2 + 1) * constr_red_2_[i].Hp2; 
                // symmetry
                Hess.block<1,2>(M_, id_1)      -= 0.5 * mult(i * 2 + 1) * constr_red_2_[i].Hp2.transpose();  
                Hess.block<2,1>(id_2, M_)      -= 0.5 * mult(i * 2 + 1) * constr_red_2_[i].Hp1; 
                // symmetry
                Hess.block<1,2>(M_, id_2)      -= 0.5 * mult(i * 2 + 1) * constr_red_2_[i].Hp1.transpose();   
        }
        
        
       
       
       
        auto time_mult = duration_cast<nanoseconds>(high_resolution_clock::now() - start_t_mult);
        // std::cout << "Hessian:\n" << Hess << std::endl;
        auto start_t_mult2 = high_resolution_clock::now();
        Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eigen_solver_M(Hess, Eigen::DecompositionOptions::EigenvaluesOnly);
                
        Eigen::VectorXd eigen_hess = eigen_solver_M.eigenvalues().real();  
        // std::cout << "Eigenvalues Hessian:\n" << eigen_hess << std::endl; 
        // std::cout << "info hessian: " << eigen_solver_M.info() << std::endl; 
        // std::cout << "Min eigenvalue Hessian:\n" << eigen_hess(0) << std::endl; 
        auto time_mult2 = duration_cast<nanoseconds>(high_resolution_clock::now() - start_t_mult2);
        
        if (debug_)
        {
            std::cout << "[HESS] Time matrix: " << (double) time_mult.count() << std::endl; 
            std::cout << "[HESS] Time EIG: " << (double) time_mult2.count() << std::endl;
        }
        
        
        return eigen_hess(0);
}

NCertRes NViewCertClass::checkOptimality(const Eigen::VectorXd & sol)
{
        /* Main function:
                1. compute multipliers
                2. compute hessian
        */
        
        // 1. Compute multipliers
        Eigen::VectorXd sol_mult; 
        auto start_t_mult = high_resolution_clock::now();
        double error_lin = computeMult(sol, sol_mult); 
        auto time_mult = duration_cast<nanoseconds>(high_resolution_clock::now() - start_t_mult);
   
      
        // 2. Form and check Hessian
        double cost_sol = sol.dot(sol); 
        
        Eigen::MatrixXd Hess;
        auto start_t_hess = high_resolution_clock::now();
        
        // reduced 
        double min_eigen = computeHessian(cost_sol, sol_mult, Hess); 
        
        auto time_hess = duration_cast<nanoseconds>(high_resolution_clock::now() - start_t_hess);
        

        if (debug_)
        {
            std::cout << "[OPT] Error linear system for multipliers: " << error_lin << std::endl; 
            std::cout << "[OPT] Minimum eigenvalue Hessian: " << min_eigen << std::endl;        
        }
        
        NCertRes res; 
        res.min_eig = min_eigen; 
        res.mult = sol_mult; 
        res.Hess = Hess; 
        res.time_mult = (double) time_mult.count(); 
        res.time_hess = (double) time_hess.count();
        return res; 

}


void NViewCertClass::printResult(const NCertRes & res)
{
        std::cout << "---------------------\n|       CERTIFIER     |\n---------------------\n";
        std::cout << "Minimum eigenvalue Hessian: " << res.min_eig << std::endl; 
        std::cout << "Time for multipliers [nanosecs]: " << res.time_mult << std::endl;
        std::cout << "Time for Hessian [nanosecs]: " << res.time_hess << std::endl;

}
}   // end of namespace NViewtrian
