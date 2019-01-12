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

void KalmanFilter::Predict() {
   x_ = F_ * x_;
   MatrixXd Ft = F_.transpose();
   P_ = F_ * P_ * Ft + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
   VectorXd y = z - H_ * x_;
   MatrixXd Ht = H_.transpose();
   MatrixXd S = H_ * P_ * Ht + R_;
   MatrixXd Si = S.inverse();
   MatrixXd K = P_ * Ht * Si;
   long ISize = x_.size();
   MatrixXd I = MatrixXd::Identity(ISize,ISize);
   // New State
   x_ = x_ + (K * y);
   P_ = (I - K*H_) * P_;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  // Calculate h(x')
  VectorXd h_x(3);
  float px_2 = pow(x_(0),2);
  float py_2 = pow(x_(1),2);
  float sqrtPx2Py2 = sqrt(px_2 + py_2);
  h_x(0) = sqrtPx2Py2;
  h_x(1) = atan2(x_(1),x_(0));
  if (fabs(h_x(0)) < 0.0001){
    h_x(2) = 0.0;
  } else {
    h_x(2) = (x_(0)*x_(2) + x_(1)*x_(3))/sqrtPx2Py2;
  }
  // Complete as per usual
  VectorXd y = z - h_x;
  // Adjust y such that phi is between -pi and pi
  if(y(1) > 3.14){
    y(1) = y(1) - 2.0*3.14;
  }else if (y(1) < -3.14){
    y(1) = y(1) + 2.0*3.14;
  }
  MatrixXd Ht = H_.transpose();
  MatrixXd S = H_ * P_ * Ht + R_;
  MatrixXd Si = S.inverse();
  MatrixXd K = P_ * Ht * Si;
  long ISize = x_.size();
  MatrixXd I = MatrixXd::Identity(ISize,ISize);
  // New State
  x_ = x_ + (K * y);
  P_ = (I - K * H_) * P_;
}
