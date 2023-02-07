#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <sstream>

// include types
#include "NViewsTypes.h"


namespace NPlanarTrian
{

        /* Struct for the result */
        struct NCertRes
        {
                double min_eig = -10;
                Eigen::VectorXd mult; 
                Eigen::MatrixXd Hess; 
                
                double time_mult = 10; 
                double time_hess = 10; 
                
                NCertRes(){};
        };  // end of struct result
        
        
        class NViewCertClass
        {
                public:
                        EIGEN_MAKE_ALIGNED_OPERATOR_NEW 
                        
                        NViewCertClass(const int M, 
                                       const std::vector<Constr2View>& constr_1, 
                                       const std::vector<Constr2View>& constr_2,  
                                       bool debug = true):
                                        M_(2*M), constr_red_1_(constr_1), 
                                        constr_red_2_(constr_2), 
                                        debug_(debug){};
                        ~NViewCertClass(){};
                        
                        /* Compute multipliers */
                        double computeMult(const Eigen::VectorXd & sol, 
                                           Eigen::VectorXd & sol_mult);
                                                                 
                                                      
                        /* Form Hessian and compute minimum eigenvalue */
                        double computeHessian(const double cost_sol, 
                                              const Eigen::VectorXd & mult, 
                                              Eigen::MatrixXd & Hess);  
                                                               
                        /* Check optimality of solution */                        
                        NCertRes checkOptimality(const Eigen::VectorXd & sol);
                        
                        /* Print results */
                        void printResult(const NCertRes & res); 
                        
                private:
                        // std::vector<Eigen::MatrixXd> constr_;  // constraints
                        std::vector<Constr2View> constr_red_1_, constr_red_2_;  // constraints 
                        int M_;  // 2 * number cameras
                        bool debug_;  // debug flag
                
        };  // end of main class



}  // end of namespace
