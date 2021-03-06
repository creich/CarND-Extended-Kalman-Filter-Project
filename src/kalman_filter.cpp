#include "kalman_filter.h"

#include <iostream>


#include <cmath>

// https://stackoverflow.com/a/11126083
// Bring the 'difference' between two angles into [-pi; pi]
template <int K, typename T>
T normalize(T rad) {
  // Copy the sign of the value in radians to the value of pi.
  T signed_pi = std::copysign(M_PI, rad);
  // Set the value of difference to the appropriate signed value between pi and -pi.
  rad = std::fmod(rad + K * signed_pi,(2 * M_PI)) - K * signed_pi;
  return rad;
}


using Eigen::MatrixXd;
using Eigen::VectorXd;

/* 
 * Please note that the Eigen library does not initialize 
 *   VectorXd or MatrixXd objects with zeros upon creation.
 */

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_ = H_in;
  R_ = R_in;
  Q_ = Q_in;
}

void KalmanFilter::Predict() {
  /**
   * TODO: predict the state
   */
  x_ = F_ * x_;
  MatrixXd Ft = F_.transpose();
  P_ = F_ * P_ * Ft + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
  /**
   * TODO: update the state by using Kalman Filter equations
   */
  VectorXd z_pred = H_ * x_;
  VectorXd y = z - z_pred;

  _commonUpdate(y);
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  /**
   * TODO: update the state by using Extended Kalman Filter equations
   */

  //MatrixXd Hj(3,4);
  // recover state parameters
  float px = x_(0);
  float py = x_(1);
  float vx = x_(2);
  float vy = x_(3);

  // check division by zero
  if( (fabs(px) < 0.0001) && (fabs(py) < 0.0001) ) {
     std::cout << "ERROR: possible division by 0! returning..." << std::endl;
     return;
  }

  // compute the Jacobian matrix
  float rho = sqrt( px*px + py*py );
  float phi = atan2(py, px);
  float rho_dot = (px*vx + py*vy) / rho;

  // y = z - h(x')  including the normalization of phi after subtraction to be in range [-pi, pi]
  VectorXd y(3);
  y << z(0)-rho, normalize<1>(z(1) - phi), z(2) - rho_dot;

  _commonUpdate(y);
}

void KalmanFilter::_commonUpdate(const VectorXd &y) {
  MatrixXd Ht = H_.transpose();
  MatrixXd S = H_ * P_ * Ht + R_;
  MatrixXd Si = S.inverse();
  MatrixXd PHt = P_ * Ht;
  MatrixXd K = PHt * Si;

  //new estimate
  x_ = x_ + (K * y);
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H_) * P_;
}
