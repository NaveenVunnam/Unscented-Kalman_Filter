#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 * This is scaffolding, do not modify
 */
UKF::UKF() {
  is_initialized_ = false;
  time_us_ = 0.0;
  
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);
  

  // initial covariance matrix
  P_ = MatrixXd(5, 5);
  P_ << 0.25, 0, 0, 0, 0,
        0, 0.25, 0, 0, 0,
        0, 0, 1, 0, 0,
        0, 0, 0, 1, 0,
        0, 0, 0, 0,1;

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 3;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = M_PI/4;
  
  //DO NOT MODIFY measurement noise values below these are provided by the sensor manufacturer.
  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;
  //DO NOT MODIFY measurement noise values above these are provided by the sensor manufacturer.
  
  /**
  TODO:

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */
  n_x_ = 5;
  n_aug_ = 7;
  lambda_ = 3 - n_x_;
  
  weights_ = VectorXd(2*n_aug_ + 1);
  double weights_0 = lambda_/(lambda_ + n_aug_);  
  double weight = 0.5/(lambda_ + n_aug_);
  weights_.fill(weight);
  weights_(0) = weights_0; 
  
  Xsig_pred_ = MatrixXd(n_x_, 2*n_aug_ + 1);
  Xsig_pred_.fill(0.0);
  R_laser_ = MatrixXd(2,2);
  R_laser_ << std_laspx_ * std_laspx_, 0,
       0, std_laspy_ * std_laspy_;
  
  R_radar_ = MatrixXd(3,3);
  R_radar_ << std_radr_ * std_radr_, 0, 0,
       0, std_radphi_ * std_radphi_, 0,
       0, 0, std_radrd_ * std_radrd_;
}
UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Make sure you switch between lidar and radar
  measurements.
  */      
  if (!is_initialized_) {    
    if ((meas_package.sensor_type_ == MeasurementPackage::RADAR) && (use_radar_ == true)) {
      cout << "Measurement came from RADAR" << endl;
      double rho = meas_package.raw_measurements_[0];
      double theta = meas_package.raw_measurements_[1];
      double rhodot = meas_package.raw_measurements_[2];
      double px = rho * cos(theta);
      double py = rho * sin(theta);
      double vx = rhodot * cos(theta);
      double vy = rhodot * sin(theta);
      double v = sqrt(vx*vx + vy*vy);
      x_ << px, py, v, 0, 0;
      
    } else if ((meas_package.sensor_type_ == MeasurementPackage::LASER) && (use_laser_ == true)) {
      cout << "Measurement came from LIDAR" << endl;
      x_ << meas_package.raw_measurements_[0], meas_package.raw_measurements_[1], 0, 0, 0;
      while (x_(0) < 0.001) x_(0) = 0.001;
      while (x_(1) < 0.001)  x_(1) = 0.001;   
     
    }
    time_us_ = meas_package.timestamp_;
    is_initialized_ = true;
    return;
  }  
  double dt_ = (meas_package.timestamp_ - time_us_)/1000000.0;
  time_us_ = meas_package.timestamp_;  
  Prediction(dt_);
  
    /*****************************************************************************
    *  Update
    ****************************************************************************/

  if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
      UpdateLidar(meas_package);
    } else if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
      UpdateRadar(meas_package);
  }
}
/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double dt) {
  /**
  TODO:

  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */
  // cout << "Prediction step started " << endl;    
  // sigma point augmentation   
  VectorXd X_aug = VectorXd(7);  //Augmented mean vector
  X_aug.head(5) = x_;
  X_aug(5) = 0;  //The mean of process noise is zero
  X_aug(6) = 0;  // ''       "                     "
  // Augmented Covariance matrix  
  MatrixXd P_aug = MatrixXd(7, 7);
  P_aug.fill(0.0);
  P_aug.topLeftCorner(n_x_, n_x_) = P_;
  P_aug(5, 5) = std_a_ * std_a_;
  P_aug(6, 6) = std_yawdd_ * std_yawdd_;
  // Square root of the Matrix  
  MatrixXd L = P_aug.llt().matrixL();
  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);
  Xsig_aug.col(0) = X_aug;
  for (int i=0; i< n_aug_; i++) {
    Xsig_aug.col(i+1) = X_aug + sqrt(lambda_ + n_aug_) * L.col(i);
    Xsig_aug.col(i+1+n_aug_) = X_aug - sqrt(lambda_ + n_aug_) * L.col(i);
  }
  // Sigma Point Predicion
  for (int i=0; i< 2*n_aug_ + 1; i++) {
    // extract and assign values
    double px = Xsig_aug(0, i);
    double py = Xsig_aug(1, i);
    double v = Xsig_aug(2, i);
    double yaw = Xsig_aug(3, i);
    double yawd = Xsig_aug(4, i);
    double nu_a = Xsig_aug(5, i);
    double nu_yawd = Xsig_aug(6, i);
    //predicted state values
    double px_p, py_p;
    //avoid division by zero
    if (fabs(yawd) > 0.001) {
        px_p = px + v/yawd * ( sin (yaw + yawd*dt) - sin(yaw));
        py_p = py + v/yawd * ( cos(yaw) - cos(yaw+yawd*dt) );
    }
    else {
        px_p = px + v*dt*cos(yaw);
        py_p = py + v*dt*sin(yaw);
    }

    double v_p = v;
    double yaw_p = yaw + yawd*dt;
    double yawd_p = yawd;
    //add noise
    
    px_p = px_p + 0.5*nu_a*dt*dt * cos(yaw);
    py_p = py_p + 0.5*nu_a*dt*dt * sin(yaw);
    v_p = v_p + nu_a*dt;
    yaw_p = yaw_p + 0.5*nu_yawd*dt*dt;
    yawd_p = yawd_p + nu_yawd*dt;
    //write predicted sigma point into right column    
    Xsig_pred_(0, i) = px_p;
    Xsig_pred_(1, i) = py_p;
    Xsig_pred_(2, i) = v_p;
    Xsig_pred_(3, i) = yaw_p;
    Xsig_pred_(4, i) = yawd_p;    
  }
  // Predict state mean  
  x_ =  Xsig_pred_*weights_ ;
    
  // Predict state covariance
  P_.fill(0.0);
  for (int i=0; i< 2*n_aug_ +1; i++) {    
    VectorXd X_diff = Xsig_pred_.col(i) - x_;    
    while (X_diff(3) > M_PI) X_diff(3) -= 2.*M_PI;
    while (X_diff(3) < -M_PI) X_diff(3) += 2.*M_PI;
    P_ = P_ + weights_(i) * X_diff * X_diff.transpose();    
  }  
}
/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */
  // Here I have used the Kalman filter equations for LASER
  int n_z = 2;
  MatrixXd H_ = MatrixXd(n_z,2*n_z+1);
  H_ << 1,0,0,0,0,
        0,1,0,0,0;  
  VectorXd y = meas_package.raw_measurements_- H_ * x_;
  MatrixXd Ht = H_.transpose();
  MatrixXd S = H_*P_*Ht + R_laser_;
  MatrixXd Si = S.inverse();
  MatrixXd K = P_ * Ht * Si;
  x_ = x_ + (K*y);
  long xsize = x_.size();
  MatrixXd I = MatrixXd::Identity(xsize,xsize);
  P_ = (I - K*H_) * P_;
  double NIS_lidar_ = y.transpose() * Si * y;
  cout << "NIS_lidar: " << NIS_lidar_ << endl; 
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */
  // Matrix for Sigma points in Measurement space  
  MatrixXd Zsig = MatrixXd(3, 2 * n_aug_ + 1);
  // mean predicted measurement  
  VectorXd Zpred = VectorXd(3);
  Zpred.fill(0.0);
  MatrixXd S = MatrixXd(3, 3);  // Measurement covariance matrix
  for (int i=0; i< 2* n_aug_ + 1; i++) {
    //Extract and assign Values
    double px = Xsig_pred_(0, i);
    double py = Xsig_pred_(1, i);
    //while (px <0.0001) px = 0.0001;
    //while (py < 0.0001) py = 0.0001;
    double v = Xsig_pred_(2, i);
    double yaw = Xsig_pred_(3, i);
    // double yawd = Xsig_pred_(4, i);
    double r = px*px + py*py;   
    double rho = sqrt(r);
    
    double phi = atan2(py, px);
    //double rhod = (px*cos(yaw)*v + py*sin(yaw)*v)/rho;
    Zsig(0, i) = rho;
    Zsig(1, i) = phi;
    Zsig(2,i) = (px*cos(yaw)*v + py*sin(yaw)*v)/ std::max(0.0001, sqrt(px*px + py*py));       
  }
  for (int i=0; i<2*n_aug_+1; i++) {
    Zpred = Zpred + weights_(i) * Zsig.col(i);
  }
  // calculate innovation covariance matrix S
  S.fill(0.0);
  for (int i=0; i<2*n_aug_+1; i++) {
    MatrixXd Z_diff = Zsig.col(i)-Zpred;
    while (Z_diff(1) > M_PI) Z_diff(1) -= 2.f*M_PI;
    while (Z_diff(1) < -M_PI) Z_diff(1) += 2.f*M_PI;
    S = S + weights_(i) * Z_diff * Z_diff.transpose();
  } 
  S = S + R_radar_;  
  // Update state
  MatrixXd T = MatrixXd(n_x_, 3);
  T.fill(0.0);
  for (int i = 0; i< 2 * n_aug_ + 1; i++) {
    VectorXd X_diff = Xsig_pred_.col(i) - x_;
    while (X_diff(3) > M_PI) X_diff(3) -= 2.f*M_PI;
    while (X_diff(3) < -M_PI) X_diff(3) += 2.f*M_PI;
    MatrixXd Z_diff = Zsig.col(i) - Zpred;
    while (Z_diff(1) > M_PI) Z_diff(1) -= 2.f*M_PI;
    while (Z_diff(1) < -M_PI) Z_diff(1) += 2.f*M_PI;
    T = T + weights_(i) * X_diff * Z_diff.transpose();
  }
  MatrixXd Si = S.inverse();  
  MatrixXd K = T * Si;
  MatrixXd Kt = K.transpose();
  // Residual 
  VectorXd Z_diff = meas_package.raw_measurements_ - Zpred;
  // Normalize the angles
  while (Z_diff(1) > M_PI) Z_diff(1) -= 2.f*M_PI;
  while (Z_diff(1) < -M_PI) Z_diff(1) += 2.f*M_PI;
  // Calculate NIS
  double NIS_radar_ = Z_diff.transpose() * Si * Z_diff;
  x_ = x_ + K * Z_diff;
  P_ = P_ - K *  S * Kt;
  //cout << "RADAR update done" << endl;
  cout << "NIS_radar: " << NIS_radar_ << endl;
}
