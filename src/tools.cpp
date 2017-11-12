#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
	VectorXd rmse(4);
	rmse << 0, 0, 0, 0;

	if (estimations.size() != ground_truth.size() || estimations.size() == 0) {
		cout << "Invalid estimation or ground truth data" << endl;
		return rmse;
	}

	for (int i = 0; i < estimations.size(); i++) {
		VectorXd residual = estimations[i] - ground_truth[i];

		residual = residual.array() * residual.array();
		rmse += residual;
	}

	rmse = rmse / estimations.size();
	rmse = rmse.array().sqrt();
	return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  MatrixXd Hj(3, 4);
  double p_x = x_state(0);
  double p_y = x_state(1);
  double v_x = x_state(2);
  double v_y = x_state(3);

  double c_1 = p_x*p_x + p_y*p_y;
  double c_2 = sqrt(c_1);
  double c_3 = (c_1*c_2);

  if (fabs(c_1) < 0.0001) {
    std::cout << "CalculateJacobian () - Error - Division by Zero" << endl;
    return Hj;
  }

  Hj << p_x/c_2, p_y/c_2, 0, 0,
         -p_y/c_1, p_x/c_1, 0, 0,
         p_y * (v_x*p_y - v_y*p_x) / c_3, p_x * (p_x*v_y - p_y*v_x) / c_3, p_x/c_2, p_y/c_2;

  return Hj;
}
