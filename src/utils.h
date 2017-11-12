#ifndef UTILS_H_
#define UTILS_H_

#include "Eigen/Dense"
#include <cmath>
#include "measurement_package.h"

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


#endif /* UTILS_H_ */