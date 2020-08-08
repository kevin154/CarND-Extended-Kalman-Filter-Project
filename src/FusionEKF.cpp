#include "FusionEKF.h"
#include <iostream>

using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::cout;
using std::endl;
using std::vector;

/**
 * Constructor.
 */
FusionEKF::FusionEKF() {
  
  is_initialized_ = false;
  previous_timestamp_ = 0;

  // initializing matrices
  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);
  H_laser_ = MatrixXd(2, 4);
  Hj_ = MatrixXd(3, 4);
  
  // measurement covariance matrix - laser
  R_laser_ << 0.0225, 0,
              0, 0.0225;

  // measurement covariance matrix - radar
  R_radar_ << 0.09, 0, 0,
              0, 0.0009, 0,
              0, 0, 0.09;

  H_laser_ << 1, 0, 0, 0,
              0, 1, 0, 0;
  
  ekf_.F_ = MatrixXd(4, 4);
  ekf_.F_ << 1, 0, 1, 0,
       0, 1, 0, 1,
       0, 0, 1, 0,
       0, 0, 0, 1;
  
  ekf_.P_ = MatrixXd(4, 4);
  ekf_.P_ << 1, 0, 1, 0,
       0, 1, 0, 1,
       0, 0, 1000, 0,
       0, 0, 0, 1000;
  
  ekf_.Q_ = MatrixXd(4, 4);
  noise_ax_ = 9;
  noise_ay_ = 9; 
}

/**
 * Destructor.
 */
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {
  /**
   * Initialization
   */
  if (!is_initialized_) {
    
    // first measurement
    ekf_.x_ = VectorXd(4);
    ekf_.x_ << 0, 0, 0, 0;
    
    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {

      double rho = measurement_pack.raw_measurements_[0];
      double phi = measurement_pack.raw_measurements_[1]; 
       
      ekf_.x_ << rho * cos(phi), 
                 rho * sin(phi), 
                 0, 
                 0;    
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      
      ekf_.x_ << measurement_pack.raw_measurements_[0], 
                 measurement_pack.raw_measurements_[1], 
                 0, 
                 0;   
    }
    
    previous_timestamp_ = measurement_pack.timestamp_;
    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }

  /**
   * Prediction
   */

  double dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;
  previous_timestamp_ = measurement_pack.timestamp_;
  double dt2 = dt * dt;
  double dt3 = 0.5 * dt2 * dt;
  double dt4 = 0.5 * dt3 * dt;
  
  ekf_.F_(0, 2) = dt;
  ekf_.F_(1, 3) = dt;
  
  ekf_.Q_ << dt4 * noise_ax_, 0, dt3 * noise_ax_, 0,
             0, dt4 * noise_ay_, 0, dt3 * noise_ay_,
             dt3 * noise_ax_, 0, dt2 * noise_ax_, 0,
             0, dt3 * noise_ay_, 0, dt2 * noise_ay_;
  
  ekf_.Predict();

  /**
   * Update
   */
  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // Radar updates
    Hj_ = tools.CalculateJacobian(ekf_.x_);
    ekf_.H_ = Hj_;
    ekf_.R_ = R_radar_;
    ekf_.UpdateEKF(measurement_pack.raw_measurements_);
  } 
  else {
    // Laser updates
    ekf_.H_ = H_laser_;
    ekf_.R_ = R_laser_;
    ekf_.Update(measurement_pack.raw_measurements_);
  }

  // print the output
  // cout << "x_ = " << ekf_.x_ << endl;
  // cout << "P_ = " << ekf_.P_ << endl;
}
