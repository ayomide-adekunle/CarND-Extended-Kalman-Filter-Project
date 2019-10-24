#include "tools.h"
#include <iostream>

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;
using std::endl;
using std::cout;
#define EP_S 0.00001

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
   * TODO: Calculate the RMSE here.
   */
  VectorXd rmse(4);
  rmse << 0,0,0,0;
  // Check if the inputs are valid
  if(estimations.size() == 0){
    cout << "The input parameter is empty" << endl;
  }
  
  else if(estimations.size() != ground_truth.size()){
    cout << "The size of the Estimations and ground truth should be same" << endl;
  }
  
  else {
    for(unsigned int j=0; j < estimations.size(); ++j){
      VectorXd residual = estimations[j] - ground_truth[j];
      residual = residual.array()*residual.array();
      rmse += residual;
    }

    // Calculate mean
    rmse = rmse / estimations.size();
    rmse = rmse.array().sqrt();
  }
  return rmse;
  
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  /**
   * TODO:
   * Calculate a Jacobian here.
   */
  // Compute the four state parameters
  float px = x_state(0);
  float py = x_state(1);
  float vx = x_state(2);
  float vy = x_state(3);
  MatrixXd Hj(3,4);
  
  // Account for unsual cases
  if (fabs(px) < EP_S and fabs(py) < EP_S){
      px = EP_S;
      py = EP_S;
  }
  
  //Compute each Jacobian elements
  float c1 = px*px+py*py;
  float c2 = sqrt(c1);
  float c3 = (c1*c2);
  
  // Compute the Jacobian matrix
  Hj << (px/c2), (py/c2), 0, 0,
       -(py/c1), (px/c1), 0, 0,
        py*(vx*py - vy*px)/c3, px*(px*vy - py*vx)/c3, px/c2, py/c2;
  return Hj;
  
}
