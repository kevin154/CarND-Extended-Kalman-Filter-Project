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

  // Validate inputs
  if (estimations.size() == 0) throw std::invalid_argument("Estimation vector must not be empty");
  if (estimations.size() != ground_truth.size()) throw std::invalid_argument("Estimation and Ground Truth vectors must be same size");

  // Container to store result
  VectorXd rmse(4);
  rmse << 0, 0, 0, 0;
    
  // Accumulate squared residuals
  for (unsigned int i = 0; i < estimations.size(); ++i) {
      VectorXd diff = estimations[i] - ground_truth[i];
      diff = diff.array() * diff.array();
      rmse += diff;
  }

  // Calculate the mean squared residual
  rmse << rmse / estimations.size();
  
  // Calculate the square root
  rmse << rmse.array().sqrt();

  // Return the result
  return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  /**
   * TODO:
   * Calculate a Jacobian here.
   */
  
  // Get state parameters
  double px, py, vx, vy, rtSumSq, sumSq;
  px = x_state(0);
  py = x_state(1);
  vx = x_state(2);
  vy = x_state(3);

  // Avoid division by 0 by setting to a minimum small positive amount
  rtSumSq = std::max(0.001, hypot(px, py));
  sumSq = rtSumSq * rtSumSq;
  
  // Continer for holding return values
  MatrixXd Hj(3, 4);
  
  // Compute the Jacobian matrix
  Hj << px / rtSumSq, py / rtSumSq, 0, 0,
        -1 * py / sumSq, px / sumSq, 0, 0, 
        py * (vx * py - vy * px) / pow(sumSq, 1.5), px * (vy * px - vx * py) / pow(sumSq, 1.5), px / rtSumSq, py / rtSumSq; 

  return Hj;
}
