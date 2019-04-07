#include "kalman_filter.h"

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
  /**
   * TODO: predict the state
   */
  x_ = F_ * x_;
  MatrixXd Ft = F_.transpose();
  P_ = F_ * P_ * Ft + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
  /**
   * TODO: update the state by using Kalman Filter equations
   */
  VectorXd y = z - H_ * x_;
  MeasUpdate(y);
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  /**
   * TODO: update the state by using Extended Kalman Filter equations
   */
  double px = x_(0);
  double py = x_(1);
  double vx = x_(2);
  double vy = x_(3);
  if(fabs(px*px+py*py) < 0.0001){
    return; 
  }
  double rho = sqrt(px*px+py*py);
  double phi = atan2(py,px);
  double rho_dot = (px*vx+py*vy)/rho;
  VectorXd hx(3);
  hx << rho, phi, rho_dot;
  VectorXd y = z - hx;
  
  // map phi error to range (-pi, pi)
  while (y(1) < -M_PI || y(1) > M_PI){
    if (y(1) < -M_PI){
      y(1) += 2 * M_PI;
    }else if(M_PI < y(1)){
      y(1) -= 2 * M_PI;
    }
  }
  
  MeasUpdate(y);
}

void KalmanFilter::MeasUpdate(const Eigen::VectorXd &y){
  MatrixXd Ht = H_.transpose();
  MatrixXd S = H_ * P_ * Ht + R_;
  MatrixXd Si = S.inverse();
  MatrixXd K =  P_ * Ht * Si;
  
  // new state
  MatrixXd I(4,4);
  I << 1,0,0,0,
       0,1,0,0,
       0,0,1,0,
       0,0,0,1;
  x_ = x_ + (K * y);
  P_ = (I - K * H_) * P_;
}
