#ifndef UTILS_H_
#define UTILS_H_

#include "Eigen/Dense"
#include <cmath>
#include "measurement_package.h"
#include "kalman_filter.h"

const int NOISE_AX_ = 9;
const int NOISE_AY_ = 9;

Eigen::VectorXd convertFromPolarToCartesian(const MeasurementPackage &measurements) {
  Eigen::VectorXd convertedValues(4);
  convertedValues << 1, 1, 1, 1;

  double ro = measurements.raw_measurements_(0);
  double theta = measurements.raw_measurements_(1);
  double ro_dot = measurements.raw_measurements_(2);

  const double tan_of_theta = std::tan(theta);
  double p_x = std::sqrt(ro*ro / (1 + tan_of_theta*tan_of_theta));
  double p_y = p_x * tan_of_theta;

  convertedValues(0) = p_x;
  convertedValues(1) = p_y;

  return convertedValues;
}

void updateFAndQMatrix(KalmanFilter &ekf, const long long previous, const long long current) {
  double dt = (current - previous) / 1000000.0;

  double dt_sq = dt * dt;
  double dt_cubed = dt_sq * dt;
  double dt_fourth = dt_sq * dt_sq;

  ekf.F_(0, 2) = dt;
  ekf.F_(1, 3) = dt;

  ekf.Q_ = Eigen::MatrixXd(4, 4);
  ekf.Q_ << dt_fourth/4*NOISE_AX_, 0, dt_cubed/2 * NOISE_AX_, 0,
            0, dt_fourth/4 * NOISE_AY_, 0, dt_cubed/2 * NOISE_AY_,
            dt_cubed/2 * NOISE_AX_, 0, dt_sq * NOISE_AX_, 0,
            0, dt_cubed/2 * NOISE_AY_, 0, dt_sq * NOISE_AY_;
}


#endif /* UTILS_H_ */