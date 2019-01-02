#include "tools.h"
#include <iostream>

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
   * TODO: Calculate the RMSE here.
   */
  VectorXd rmse(4);
  rmse << 0,0,0,0;

  // check the validity of the following inputs:
  //  * the estimation vector size should not be zero
  if( (estimations.size() == 0) || (ground_truth.size() == 0) ) {
     std::cerr << "ERROR: empty estimations OR ground_truth vector! .. quitting" << std::endl;
     return rmse;
  }

  //  * the estimation vector size should equal ground truth vector size
  if( estimations.size() != ground_truth.size() ) {
     std::cerr << "ERROR: estimations vector and ground_truth vector have different sizes! .. quitting" << std::endl;
     return rmse;
  }

  // accumulate squared residuals
  VectorXd squared_res(4);
  squared_res << 0, 0, 0, 0;

  for (int i=0; i < estimations.size(); ++i) {
     VectorXd c = estimations[i] - ground_truth[i];
     VectorXd s_res = c.array() * c.array();
     squared_res = squared_res + s_res;
  }

  // calculate the mean
  VectorXd mean = squared_res / estimations.size();

  // calculate the squared root
  rmse = mean.array().sqrt();

  // return the result
  return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  /**
   * TODO:
   * Calculate a Jacobian here.
   */

   MatrixXd Hj(3,4);
   // recover state parameters
   float px = x_state(0);
   float py = x_state(1);
   float vx = x_state(2);
   float vy = x_state(3);

   // check division by zero
   if( (px == 0) && (py == 0)) {
      std::cout << "ERROR: possible division by 0! returning initial Hj" << std::endl;
      return Hj;
   }

   // compute the Jacobian matrix
   float px2py2 = px*px + py*py;

   Hj(0, 0) = px / sqrt(px2py2);
   Hj(0, 1) = py / sqrt(px2py2);
   Hj(0, 2) = 0;
   Hj(0, 3) = 0;

   Hj(1, 0) = -1 * (py / px2py2);
   Hj(1, 1) = px / px2py2;
   Hj(1, 2) = 0;
   Hj(1, 3) = 0;

   Hj(2, 0) = (py*(vx*py - vy*px)) / pow(px2py2, 3/2);
   Hj(2, 1) = (px*(vy*px - vx*py)) / pow(px2py2, 3/2);
   Hj(2, 2) = px / sqrt(px2py2);
   Hj(2, 3) = py / sqrt(px2py2);

   return Hj;
}
