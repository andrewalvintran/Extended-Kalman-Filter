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
  x_ = F_ * x_;
  MatrixXd F_transpose = F_.transpose();
  P_ = F_ * P_ * F_transpose + Q_;
}

void KalmanFilter::UpdateXAndP(const MatrixXd &K, const MatrixXd &y) {
  x_ = x_ + (K * y);
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H_) * P_;
}

MatrixXd KalmanFilter::GenerateK() const {
  MatrixXd H_transpose = H_.transpose();
  MatrixXd S = H_ * P_ * H_transpose + R_;
  MatrixXd S_inv = S.inverse();
  MatrixXd P_H_transpose = P_ * H_transpose;
  MatrixXd K = P_H_transpose * S_inv;

  return K;
}

void KalmanFilter::Update(const VectorXd &z) {
  VectorXd z_pred = H_ * x_;
  VectorXd y = z - z_pred;
  MatrixXd K = GenerateK();

  // new estimate
  UpdateXAndP(K, y);
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  VectorXd polar_coords(3);

  double ro = std::sqrt(x_(0)*x_(0) + x_(1)*x_(1));
  double theta = std::atan(x_(1)/x_(0));
  double ro_dot = (x_(0)*x_(2) + x_(1)*x_(3)) / ro;

  polar_coords << ro, theta, ro_dot;

  VectorXd z_pred = polar_coords;
  VectorXd y = z - z_pred;

  // make sure y(1) is between -pi and pi
  const double PI = std::acos(0.0);
  while (y(1) > PI || y(1) < -PI) {
    if (y(1) > PI) {
      y(1) -= 2*PI;
    } else {
      y(1) += 2*PI;
    }
  }

  MatrixXd K = GenerateK();
  // new estimate
  UpdateXAndP(K, y);
}
