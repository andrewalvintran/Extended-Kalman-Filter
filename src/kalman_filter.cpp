#include "kalman_filter.h"
#include <cmath>

using Eigen::MatrixXd;
using Eigen::VectorXd;

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
  TODO:
    * predict the state
  */
  x_ = F_ * x_;
  MatrixXd F_transpose = F_.transpose();
  P_ = F_ * P_ * F_transpose + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
  /**
  TODO:
    * update the state by using Kalman Filter equations
  */
  VectorXd z_pred = H_ * x_;
  VectorXd y = z - z_pred;
  MatrixXd H_transpose = H_.transpose();
  MatrixXd S = H_ * P_ * H_transpose + R_;
  MatrixXd S_inv = S.inverse();
  MatrixXd P_H_transpose = P_ * H_transpose;
  MatrixXd K = P_H_transpose * S_inv;

  // new estimate
  x_ = x_ + (K * y);
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H_) * P_;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  /**
  TODO:
    * update the state by using Extended Kalman Filter equations
  */
  VectorXd polar_coords(3);

  double ro = std::sqrt(x_(0)*x_(0) + x_(1)*x_(1));
  double theta = std::atan(x_(1)/x_(0));
  double ro_dot = (x_(0)*x_(2) + x_(1)*x_(3)) / ro;

  polar_coords << ro, theta, ro_dot;

  VectorXd z_pred = polar_coords;
  VectorXd y = z - z_pred;
  MatrixXd H_transpose = H_.transpose();
  MatrixXd S = H_ * P_ * H_transpose + R_;
  MatrixXd S_inv = S.inverse();
  MatrixXd P_H_transpose = P_ * H_transpose;
  MatrixXd K = P_H_transpose * S_inv;

  // new estimate
  x_ = x_ + (K * y);
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H_) * P_;
}
