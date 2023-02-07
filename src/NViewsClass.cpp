#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <sstream>

// include types
#include "NViewsTypes.h"
// utils
#include "NViewsUtils.h"
// header
#include "NViewsClass.h"

#include <chrono>  // timer

using namespace std::chrono;

namespace NPlanarTrian
{
// create problem matrices
void NViewsClass::createProblemMatrices(const std::vector<PairObj> & obj,
                                        const int N_cams)
{
  
        constr_red_1_.empty(); 
        constr_red_2_.empty(); 
        
        // size of the graph
        M_ = obj.size();  // each link in the graph gives 2 constraints
        N_cams_ = N_cams;
        
        for (int i=0; i<M_; i++)
        {       
                // std::cout << "Going for index " << i << std::endl; 
                Matrix3 H = obj[i].H; 
                Vector3 p1 = obj[i].p1; 
                Vector3 p2 = obj[i].p2;      
                
                /*
                Vector3 qq = H * p1; 
                qq /= qq(2); 
                std::cout << "p2:\n" << p2 << std::endl; 
                std::cout << "H * p1:\n" << qq << std::endl; 
                */
                
                // Constraint #1
                Matrix3 T1 = Matrix3::Zero(); 
                T1 << 0, 0, 0, 0, 0, -1, 0, 1, 0;
                Matrix3 H1 = T1.transpose() * H;               
                
                Vector3 Hp2 = H1.transpose() * p2; 
                Vector3 Hp1 = H1 * p1; 
                
                double e_ep = p2.dot(H1 * p1);
                
                int id_1 = (obj[i].id1) * 2; 
                int id_2 = (obj[i].id2) * 2; 
                
                // std::cout << "Error ep 1: " << e_ep << std::endl;
                
                
               
                // create struct
                Constr2View ci; 
                ci.b = e_ep; 
                ci.Hp2 = Hp2.topRows(2); 
                ci.Hp1 = Hp1.topRows(2); 
                ci.H = H1.transpose().topLeftCorner(2,2); 
                ci.id_1 = id_1; 
                ci.id_2 = id_2;
                constr_red_1_.push_back(ci);

                // Constraint #2
                Matrix3 T2 = Matrix3::Zero(); 
                T2 << 0, 0, 1, 0, 0, 0, -1, 0, 0;
                Matrix3 H2 = T2.transpose() * H;               
                
                Hp2 = H2.transpose() * p2; 
                Hp1 = H2 * p1; 
                
                e_ep = p2.dot(H2 * p1);
                                
                // std::cout << "Error ep 2: " << e_ep << std::endl;
                
               
               
                // create struct
                Constr2View ci_2; 
                ci_2.b = e_ep; 
                ci_2.Hp2 = Hp2.topRows(2); 
                ci_2.Hp1 = Hp1.topRows(2); 
                ci_2.H = H2.transpose().topLeftCorner(2,2); 
                ci_2.id_1 = id_1; 
                ci_2.id_2 = id_2;
                constr_red_2_.push_back(ci_2);

        }
        
        constr_are_created_ = true;

}

// correction: init
double NViewsClass::initCorrection(Eigen::VectorXd & sol_init, 
                                   Eigen::MatrixXd & A, 
                                   Eigen::VectorXd & b)
{
        
        A.resize(M_ * 2, 2 * N_cams_); 
        b.resize(M_ * 2);
        
        A.setZero(); 
        b.setZero(); 
        
        auto start_t_init = high_resolution_clock::now();
        
        for (int i=0; i < M_; i++)
        {                
                
                // first constraint
                A.block<1,2>(i * 2, constr_red_1_[i].id_1) = constr_red_1_[i].Hp2.transpose(); 
                A.block<1,2>(i * 2, constr_red_1_[i].id_2) = constr_red_1_[i].Hp1.transpose();
                b(i * 2) = constr_red_1_[i].b;
                
                // first constraint
                A.block<1,2>(i * 2 + 1, constr_red_2_[i].id_1) = constr_red_2_[i].Hp2.transpose(); 
                A.block<1,2>(i * 2 + 1, constr_red_2_[i].id_2) = constr_red_2_[i].Hp1.transpose();                
                b(i * 2 + 1) = constr_red_2_[i].b;
                // std::cout << "constr_red_2_[i].b:\n" << constr_red_2_[i].b << std::endl;
                
        }
        
       
        // std::cout << "Vector b:\n" << b << std::endl; 
        auto time_init = duration_cast<nanoseconds>(high_resolution_clock::now() - start_t_init);
        
        // std::cout << "[INIT] Matrix constr: " << (double) time_init.count() << std::endl; 
        
        auto start_t_init2 = high_resolution_clock::now();
        // sol_init.resize(2 * N_cams_);
        double error_sol = 1; 
        
        if (A.cols() <= A.rows())
                error_sol = solveLinearSystemMinNorm(A.transpose()*A, - A.transpose()*b, sol_init);
                
        else
                error_sol = solveLinearSystemMinNorm(A, - b, sol_init);
       
        // double error_sol = solveLinearSystem( A, - b, sol_init);
        
        
        auto time_init2 = duration_cast<nanoseconds>(high_resolution_clock::now() - start_t_init2);
        
   
        // std::cout << "[INIT] Matrix solve: " << (double) time_init2.count() << std::endl; 
        
        // std::cout << "sol:\n" << sol_init << std::endl; 
        return error_sol;
        
}



// correction ref
// NOTE: lower value of min hessian
double NViewsClass::refineCorrection(const Eigen::MatrixXd & A0, 
                                     const Eigen::VectorXd & b0, 
                                     const Eigen::VectorXd & sol_init, 
                                     Eigen::VectorXd & sol_ref)
{

        Eigen::MatrixXd C(M_ * 2, 2 * N_cams_); 
        Eigen::VectorXd e(M_ * 2);
        C = A0; 
        e = b0; 
        
        
        auto start_t_init = high_resolution_clock::now();
   
        for (int i=0; i < M_ ; i++)
        {
               
                
                Vector2 p1 = sol_init.block<2,1>(constr_red_1_[i].id_1, 0);
                Vector2 p2 = sol_init.block<2,1>(constr_red_1_[i].id_2, 0);
                
                
                // constraint #1
                C.block<1,2>(i * 2, constr_red_1_[i].id_1) += p2.transpose() * constr_red_1_[i].H.transpose(); 
                C.block<1,2>(i * 2, constr_red_1_[i].id_2) += p1.transpose() * constr_red_1_[i].H; 
                               
                e(i * 2) -= p1.transpose() * constr_red_1_[i].H * p2;
                
                
                // constraint #2
                C.block<1,2>(i * 2 + 1, constr_red_1_[i].id_1) += p2.transpose() * constr_red_2_[i].H.transpose(); 
                C.block<1,2>(i * 2 + 1, constr_red_1_[i].id_2) += p1.transpose() * constr_red_2_[i].H; 
                               
                e(i * 2 + 1) -= p1.transpose() * constr_red_2_[i].H * p2;

        }
        
        auto time_init = duration_cast<nanoseconds>(high_resolution_clock::now() - start_t_init);
       
        double error_sol_t = 1;
        
        auto start_t_init_t = high_resolution_clock::now();
        
        if (C.cols() <= C.rows())
                error_sol_t = solveLinearSystemMinNorm(C.transpose()*C, - C.transpose()*e, sol_ref);
                
        else
                error_sol_t = solveLinearSystemMinNorm(C, - e, sol_ref);
        
        auto time_init_t = duration_cast<nanoseconds>(high_resolution_clock::now() - start_t_init_t);
           
       
        return error_sol_t;
        
}   


// correction: K-th refinement
// main function
NViewsResult NViewsClass::correctObservations(NViewsOptions & options)
{       

        if (!constr_are_created_)
        {
                std::cout << "[ERRROR] Constraints are not created!\n";
                return NViewsResult();
        }
        // else        
                
        Eigen::VectorXd sol_i, sol_init;
        Eigen::MatrixXd A0; 
        Eigen::VectorXd b0; 
        
        if (options.debug)
            std::cout << "[CORR] Estimating initial guess\n"; 
            
        auto start_t_init = high_resolution_clock::now();
        double error_init = initCorrection(sol_init, A0, b0);
        auto time_init = duration_cast<nanoseconds>(high_resolution_clock::now() - start_t_init);
        
        // std::cout << "Size A ref: " << A0.rows() << "," << A0.cols() << std::endl; 
        // std::cout << "Checking constraints\n"; 
        Eigen::VectorXd val_constr_i_1, val_constr_i_2; 
        double max_constr_val_i = 10.0, sq_constr_i = 10.0;
        
        auto start_t_ref_i2 = high_resolution_clock::now();
        
        // double tot_constr_i = checkConstraints(sol_init, constr_exp_, val_constr_i, max_constr_val_i);
        double tot_constr_i = checkConstraints(sol_init, constr_red_1_, constr_red_2_, val_constr_i_1, val_constr_i_2, max_constr_val_i, sq_constr_i);
        
        auto time_ref_i2 = duration_cast<nanoseconds>(high_resolution_clock::now() - start_t_ref_i2);
        
        int k_iter = 0; 
        
        if (options.debug)
        {
            std::cout << "[CORR] Maximum value constraint for initial guess: " << max_constr_val_i << std::endl; 
            std::cout << "[CORR] L-1 norm constraints for initial guess: " << tot_constr_i << std::endl;   
        }
            
        // std::cout << "[MAIN] Time constraints: " << (double) time_ref_i2.count() << std::endl;
        
        
        
        NViewsResult res = NViewsResult(); 
        res.tot_constr_init = tot_constr_i; 
        res.max_constr_init = max_constr_val_i; 
        res.sq_constr_init = sq_constr_i;
        if (options.record_constr)
                {
                        res.rec_constr_1.empty(); 
                        res.rec_constr_1.push_back(val_constr_i_1);
                        res.rec_constr_2.empty(); 
                        res.rec_constr_2.push_back(val_constr_i_2);
                
                }
        
        if (options.debug)
        {
            std::cout << "[CORR] Refining solution\n--- --- --- --- --- --- --- --- ---\n"; 
        }
        
        sol_i = sol_init;
        double time_ref = 0;
        double diff_sol = 10; 
        double error_lin = 10; 
        // evaluate stopping condition 
        // We sop if
        // 1. we reach the maximum number of iterations
        // 2. The diff between solutions is less than threshold
        while ( (k_iter < options.max_iters) & (diff_sol > options.max_diff_sol) )
        {
 
                // keep refining
                Eigen::VectorXd sol_next = sol_i;
              
                auto start_t_ref_i = high_resolution_clock::now();
                error_lin = refineCorrection(A0, b0, sol_i, sol_next);
                auto time_ref_i = duration_cast<nanoseconds>(high_resolution_clock::now() - start_t_ref_i);
                
                diff_sol = (sol_next - sol_i).squaredNorm(); 
                
                // std::cout << "Checking constraints\n";
                tot_constr_i = checkConstraints(sol_next, constr_red_1_, constr_red_2_, val_constr_i_1, val_constr_i_2, max_constr_val_i, sq_constr_i);
                
                if (options.debug)
                {
                    std::cout << "[CORR] Iteration #" << k_iter; 
                    std::cout << "            Max constraint: " << max_constr_val_i; 
                    std::cout << "            Total constraint L1: " << tot_constr_i;
                    std::cout << "            Total constraint L2: " << sq_constr_i;
                    std::cout << "            Norm prev solution: " << sol_i.squaredNorm(); 
                    std::cout << "            Norm new solution: " << sol_next.squaredNorm() << std::endl;
                    std::cout << "            Diff between solutions: " << diff_sol << std::endl; 
                }
                
                
                if (options.record_constr)
                {
                        res.rec_constr_1.push_back(val_constr_i_1);
                        res.rec_constr_2.push_back(val_constr_i_2);
                
                }
                time_ref += (double) time_ref_i.count();
                // update iteration
                k_iter++;        
                sol_i = sol_next;
                
        }
        
       
        if (options.debug)
        {
            std::cout << "\n--- --- --- --- --- --- --- --- ---\n";
            std::cout << "[CORR] Nr. iterations: " << k_iter; 
            std::cout <<      "       Max constraint: " << max_constr_val_i; 
            std::cout << "            Total constraint L1: " << tot_constr_i;
            std::cout << "            Total constraint L2: " << sq_constr_i;
            std::cout << "            Norm intial solution: " << sol_init.squaredNorm(); 
            std::cout <<   "          Norm final solution: " << sol_i.squaredNorm() << std::endl;
            std::cout << "            Diff between solutions: " << diff_sol << std::endl; 
        }
       
        res.n_iters = k_iter;
        res.sol_init = sol_init;
        res.sol_final = sol_i;
        res.max_constr = max_constr_val_i; 
        res.tot_constr = tot_constr_i; 
        res.sq_constr = sq_constr_i; 
        res.diff_sol = diff_sol; 
        if (options.save_val_constr)
        {
            res.all_constr_1 = val_constr_i_1;
            res.all_constr_2 = val_constr_i_2;
        }
        res.time_init = (double) time_init.count();
        res.time_ref = time_ref;
        res.error_lin = error_lin;
        
        return res;
}

// Print solution
void NViewsClass::printResult(NViewsResult & res_corr)
{

        std::cout << "---------------------\n|       SOLUTION    |\n---------------------\n"; 
        std::cout << "Number iterations: " << res_corr.n_iters << std::endl;
        // std::cout << "Initial solution:\n" << res_corr.sol_init << std::endl; 
        // std::cout << "Final solution:\n" << res_corr.sol_final << std::endl; 
        std::cout << "Maximum value of constraint: "  << res_corr.max_constr << std::endl; 
        std::cout << "Norm of the constraints [L1]: " << res_corr.tot_constr << std::endl; 
        std::cout << "Norm of the constraints L2: "   <<  res_corr.sq_constr << std::endl;
        std::cout << "Difference between the last two solutions: " << res_corr.diff_sol << std::endl;
        std::cout << "Error linear system: " << res_corr.error_lin << std::endl; 
        if (res_corr.all_constr_1.size() > 0)
        {
                std::cout << "value of all the constraints #1:\n" << res_corr.all_constr_1 << std::endl;
                std::cout << "value of all the constraints #2:\n" << res_corr.all_constr_2 << std::endl;
        }
        
        if (res_corr.rec_constr_1.size() > 0)
                {
                        std::cout << "The value of the constraints for the iterations are:\n"; 
                        for (int j = 0; j < res_corr.rec_constr_1.size(); j++)
                        {
                                std::cout << "Iteration #" << j << std::endl; 
                                std::cout << "C1: " << res_corr.rec_constr_1[j].transpose() << std::endl;
                                std::cout << "C2: " << res_corr.rec_constr_2[j].transpose() << std::endl;
                        }
                
                }
         std::cout << "Time init [nanosecs]: " << res_corr.time_init << std::endl; 
         std::cout << "Time ref [nanosecs]: " << res_corr.time_ref << std::endl;

}



















}  // end of namespace
