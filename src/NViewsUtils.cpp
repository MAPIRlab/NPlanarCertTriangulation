#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <sstream>

// include types
#include "NViewsTypes.h"

// include header
#include "NViewsUtils.h"

#define EIGEN_USE_LAPACKE_STRICT
// for eigendecomposition
#include <eigen3/Eigen/Dense>
#include <Eigen/Eigenvalues> 

#include <chrono>  // timer



using namespace std::chrono;
using namespace std;


namespace NPlanarTrian
{


// Check constraints with reduced constraints
double checkConstraints(const Eigen::VectorXd & sol, 
                        std::vector<Constr2View> & constr_1, 
                        std::vector<Constr2View> & constr_2, 
                        Eigen::VectorXd & val_constr_1, 
                        Eigen::VectorXd & val_constr_2,
                        double & max_constr_val, 
                        double & sq_constr_val)
{
        // Evaluate the constraints at the given solution
        int M = constr_1.size(); 
        val_constr_1.resize(M); 
        val_constr_2.resize(M); 
        val_constr_1.setZero(); 
        val_constr_2.setZero();
        
        double tot_constr = 0.0; 
        max_constr_val = -10.0;
        sq_constr_val = 0.0;
        for (int i=0; i < M; i++)
        {
                // double val_i = w.dot(constr[i] * w);
                Vector2 p1 = sol.block<2,1>(constr_1[i].id_1, 0);
                Vector2 p2 = sol.block<2,1>(constr_1[i].id_2, 0);
                
                
                double val_i_1 = constr_1[i].b + p1.transpose() * constr_1[i].H * p2 + p1.dot(constr_1[i].Hp2) + p2.dot(constr_1[i].Hp1); 
                double val_i_2 = constr_2[i].b + p1.transpose() * constr_2[i].H * p2 + p1.dot(constr_2[i].Hp2) + p2.dot(constr_2[i].Hp1);                 
                
                
                val_constr_1(i) = val_i_1;
                val_constr_2(i) = val_i_2; 
                double abs_val_i_1 = (val_i_1 > 0) ? val_i_1 : -val_i_1;
                double abs_val_i_2 = (val_i_2 > 0) ? val_i_2 : -val_i_2;
                tot_constr += abs_val_i_1;
                tot_constr += abs_val_i_2; 
                
                sq_constr_val += abs_val_i_1 * abs_val_i_1;
                sq_constr_val += abs_val_i_2 * abs_val_i_2;
                
                // std::cout << "Tot constr: " << tot_constr << std::endl;
                if (abs_val_i_1 > max_constr_val)
                        max_constr_val = abs_val_i_1;
                if (abs_val_i_2 > max_constr_val)
                        max_constr_val = abs_val_i_2;
        }
        // std::cout << "Squared norm of the constraints: " << sq_constr_val << std::endl;

        return tot_constr;
}

        

// solve linear system minimum norm
double solveLinearSystemMinNorm(const Eigen::MatrixXd & A, 
                                const Eigen::VectorXd & b,
                                Eigen::VectorXd & sol)
{
     
        double tol_rank = 1e-05;  // 1e-03 works fine with simple approach
        
        Eigen::VectorXd y_sol;
        
        // std::cout << "[MIN NORM] Dimension A: " << A.rows() << " x " << A.cols() << std::endl;
        auto start_t_init = high_resolution_clock::now();
        Eigen::CompleteOrthogonalDecomposition<Eigen::MatrixXd> cod(A.rows(),
                                                A.cols());
        cod.setThreshold(tol_rank);
        cod.compute(A);
        auto time_init = duration_cast<nanoseconds>(high_resolution_clock::now() - start_t_init);
        
        // std::cout << "[MIN NORM] Rank A: " << cod.rank() << std::endl; 
        // std::cout << "Size A: " << A.rows() << "," << A.cols() << std::endl;
        // std::cout << "[LIN] Decom A: " << (double) time_init.count() << std::endl; 
        
        sol.resize(A.cols(), 1); 
        
        auto start_t_init2 = high_resolution_clock::now();
        sol = cod.solve(b);
        
        auto time_init2 = duration_cast<nanoseconds>(high_resolution_clock::now() - start_t_init2);
        // std::cout << "[LIN] Solve system: " << (double) time_init2.count() << std::endl;
        
        Eigen::VectorXd e_lin = A * sol - b;
        // std::cout << "error LS: " << e_lin.squaredNorm() << std::endl;
        double error_sol = e_lin.squaredNorm();
        
        // sol = y_sol;
        
        
        
        return error_sol;
}        


double solveLinearSystem(const Eigen::MatrixXd & A,  const Eigen::VectorXd & b, Eigen::VectorXd & sol)
{
    Eigen::VectorXd y_sol(A.cols());
     
    y_sol = A.colPivHouseholderQr().solve(b); 
     Eigen::VectorXd e_lin = A * y_sol - b;
        // std::cout << "error LS: " << e_lin.squaredNorm() << std::endl;
        double error_sol = e_lin.squaredNorm();
        sol.resize(A.cols(), 1); 
        sol = y_sol;
        
        
        
        return error_sol;

}

       

double triangulateNPoint(const std::vector<Eigen::Matrix<double, 3, 4>> & proj_s,
                         const std::vector<Eigen::Vector3d> & obs_s, 
                         Eigen::Vector3d & P_3d, 
                         Eigen::VectorXd & depths)
{

        std::vector<Matrix4> P4 = {};
        int N_cams = proj_s.size(); 
        P4.reserve(N_cams); 
        
        for (int i=0; i < N_cams; i++)
        {
                Eigen::Matrix4d Pi = Eigen::Matrix4d::Identity();
                
                Pi.block<3,4>(0, 0) = proj_s[i];
                P4.push_back(Pi);
        }
        
        return (triangulateNPoint(P4, obs_s, P_3d, depths));

}             
// triangulate point
double triangulateNPoint(const std::vector<Matrix4> & proj_s,
                         const std::vector<Eigen::Vector3d> & obs_s, 
                         Eigen::Vector3d & P_3d, 
                         Eigen::VectorXd & depths)
{

        int N_cams = proj_s.size();
        Eigen::MatrixXd P(N_cams * 3, N_cams + 4);
        P.setZero();
        for (int i=0; i < N_cams; i++)
        {
                const int ri = (i)*3;
                // int rf = (i-1)*3 + 2; 
                // std::cout << "Projection:\n" << proj_s[i] << std::endl;
                P.block<3,4>(ri, 0) = proj_s[i].block<3,4>(0,0);
                // std::cout << "obs:\n" << obs_s[i] << std::endl;
                P.block<3,1>(ri, 4+i) = obs_s[i]; 
        }
        
        // std::cout << "Matrix P:\n" << P << std::endl; 
        // std::cout << "Size P: " << P.rows() << " x " << P.cols() << std::endl; 
        Eigen::JacobiSVD<Eigen::MatrixXd> svd(P, Eigen::ComputeFullV);
        // std::cout << "Singular values P:\n" << svd.singularValues() << std::endl; 
        P_3d.setZero(); 
        // std::cout << "Matrix V:\n" << svd.matrixV() << std::endl;
        Eigen::Vector4d P_hom = svd.matrixV().topRightCorner(4, 1);
        // std::cout << "Solution X:\n" << P_hom << std::endl;
        P_3d = P_hom.topRows(3) / P_hom(3);
        // std::cout << "3d point:\n" << P_3d << std::endl;
                
        depths.resize(N_cams); 
        depths.setZero();
        depths = svd.matrixV().bottomRightCorner(N_cams, 1);
        
        // std::cout << "depths:\n" << depths << std::endl;
        int last_id = std::min(N_cams * 3, N_cams + 4) - 1; 
        // std::cout << "Last index: " << last_id << std::endl;
        
        double error_lin = svd.singularValues()(last_id); 
        // std::cout << "Error lin: " << error_lin << std::endl;

        return error_lin;
}



}   // end of namespace NViewtrian

