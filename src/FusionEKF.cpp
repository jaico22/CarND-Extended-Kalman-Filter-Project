#include "FusionEKF.h"
#include <iostream>
#include "Eigen/Dense"
#include "tools.h"

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
  F_ = MatrixXd(4, 4);
  P_ = MatrixXd(4, 4);
  x_init_ = VectorXd(4);
  Q_ = MatrixXd(4,4);
  //measurement covariance matrix - laser
  R_laser_ << 0.0225, 0,
  0, 0.0225;
  //measurement covariance matrix - radar
  R_radar_ << 0.09, 0, 0,
  0, 0.0009, 0,
  0, 0, 0.09;
  cout << R_radar_(0,0) << endl;
  // Measurement Matrix  - Laser
  H_laser_ << 1, 0, 0, 0,
  0, 1, 0, 0;
  // State Transition
  ekf_.F_ = MatrixXd(4,4);
  ekf_.F_ << 1, 0, 1, 0,
  0, 1, 0, 1,
  0, 0, 1, 0,
  0, 0, 0, 1;
  // Process Noise

  ekf_.P_ = MatrixXd(4,4);
  ekf_.P_ << 1, 0, 0, 0,
  0, 1, 0, 0,
  0, 0, 1000, 0,
  0, 0, 0, 1000;


  noise_ax = 9.0;
  noise_ay = 9.0;
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
    // Set matrices

    cout << "ekf_... " << endl;
    ekf_.x_ = VectorXd(4);
    ekf_.Q_ = MatrixXd(4,4);
    cout << R_radar_(0,0) << endl;
    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      // Convert from polar to cartesian
      cout << "Radar detected" << endl;
      float rho = measurement_pack.raw_measurements_(0);
      float phi = measurement_pack.raw_measurements_(1);
      float rho_dot = measurement_pack.raw_measurements_(2);
      float sinPhi = sin(phi);
      float cosPhi = cos(phi);

      cout << "Doneee" << endl;

      ekf_.x_(0) = rho * cosPhi;
      ekf_.x_(1) = rho * sinPhi;
      ekf_.x_(2) = rho_dot * cosPhi;
      ekf_.x_(3) = rho_dot * sinPhi;
    }else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      // Already in cartesian, not a whole lot to do here...
      int size_ = x_init_.size();
      cout << size_ << endl;
      cout << "LASER" << endl;
      ekf_.x_(0) = measurement_pack.raw_measurements_(0);
      ekf_.x_(1) = measurement_pack.raw_measurements_(1);
      ekf_.x_(2) = 0.0;
      ekf_.x_(3) = 0.0;
    }
    cout << "Processed" << endl;

    // done initializing, no need to predict or update
    previous_timestamp_ = measurement_pack.timestamp_;
    is_initialized_ = true;
    cout << "Initialized" << endl;
    return;
  }

  /**
  * Prediction
  */
  // Calculate time since last measurementa and conver to secodns.
  float dt = (measurement_pack.timestamp_ - previous_timestamp_)/1000000.0;
  previous_timestamp_ = measurement_pack.timestamp_;
  // Update values in state transition matrix to reflect time difference
  int size_ = ekf_.F_.size();
  cout << size_ << endl;
  ekf_.F_(0,2) = dt;
  ekf_.F_(1,3) = dt;
  // Update noise process covariance matrix
  float dt2 = pow(dt,2);
  float dt3 = pow(dt,3);
  float dt4 = pow(dt,4);
  float dt4_4 = dt4/4.0;
  float dt3_2 = dt3/2.0;
  cout << "---- Noise AX ----" << endl;
  cout << noise_ax << endl;
  ekf_.Q_ << dt4_4*noise_ax,0,dt3_2*noise_ax,0,
  0,dt4_4*noise_ay,0,dt3_2*noise_ay,
  dt3_2*noise_ax,0,dt2*noise_ax,0,
  0,dt3_2*noise_ay,0,dt2*noise_ay;
  // Predict!!
  ekf_.Predict();

  /**
  * Update
  */

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    Tools tools;
    // Calculate Jacobian and update measurement model
    Hj_ = tools.CalculateJacobian(ekf_.x_);
    ekf_.H_ = Hj_;
    // Use Radar noise covariance
    ekf_.R_ = R_radar_;
    // Update using Extended Kalman Filter
    ekf_.UpdateEKF(measurement_pack.raw_measurements_);
  } else {
    // Use Laser measurement model and noise covariance
    ekf_.H_ = H_laser_;
    ekf_.R_ = R_laser_;
    // Update using normal kalman filter
    ekf_.Update(measurement_pack.raw_measurements_);

  }

  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}
