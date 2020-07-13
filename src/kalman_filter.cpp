#include "kalman_filter.h"
#include <iostream>

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

// Ensures given input angle is within [-pi, pi] radians
// Recursively moves input angle within range, should require <=2 calls
double normaliseAngle(double angle){
  
  if (angle > M_PI){
      angle -= 2*M_PI;
      normaliseAngle(angle);
  }
  else if (angle < -M_PI){
      angle += 2*M_PI;
      normaliseAngle(angle);
  }
  return angle;
}

void KalmanFilter::Predict() {
  /**
   * TODO: predict the state
   */
  x_ = F_ * x_;
  P_ = F_ * P_ * F_.transpose() + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
  /**
   * TODO: update the state by using Kalman Filter equations
   */
   VectorXd zpred = H_ * x_; 
   VectorXd y = z - zpred;
   MatrixXd S = H_ * P_ * H_.transpose() + R_;
   MatrixXd K = P_ * H_.transpose() * S.inverse();
   
   x_ = x_ + (K * y);
   MatrixXd I = MatrixXd::Identity(x_.size(), x_.size());
   P_ = (I - K * H_) * P_;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  /**
   * TODO: update the state by using Extended Kalman Filter equations
   */
  
  // Position variables
  double px, py, vx, vy;
  px = x_[0];
  py = x_[1];
  vx = x_[2];
  vy = x_[3];
 
  // Convert from cartesian to polar
  double rho, phi, rhodot;
  // Ensure no division by 0, set rho to a minimum positive number
  rho = std::max(0.000001, hypot(px, py));
  // atan2 returns a number within [-pi, pi]  
  phi = atan2(py, px);
  rhodot = (px * vx + py * vy) / rho;
  
  // Container for polar variables
  VectorXd hx = VectorXd(3);
  hx << rho, phi, rhodot;
  
  // Calculate difference and normalise angle to ensure it is within [-pi, pi]  
  VectorXd y = z - hx;
  y(1) = normaliseAngle(y(1));
  
  MatrixXd S = H_ * P_ * H_.transpose() + R_;
  MatrixXd K = P_ * H_.transpose() * S.inverse();
  MatrixXd I = MatrixXd::Identity(x_.size(), x_.size());
  
  // Update position
  x_ = x_ + (K * y);
  P_ = (I - K * H_) * P_;
}
