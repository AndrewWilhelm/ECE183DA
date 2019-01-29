// Note: this uses the MatrixMath library found here: http://playground.arduino.cc/Code/MatrixMath
// The MatrixMath code needs to be added to the Arduino library folder before running

#include <MatrixMath.h>

struct kalman_filter {
  int m = 1; // number of state variables
  int n = 1; // number of sensor observations
  int p = 1; // number of inputs into system
  
  mtx_type[m][1] x; // m x 1
  mtx_type[n][1] z; // n x 1
  mtx_type[p][1] u; // p x 1
  
  
  mtx_type[m][m] A; // m x m
  mtx_type[m][p] B; // m x p
  mtx_type[n][m] C; // n x m
  mtx_type[m][n] K; // m x n the kalman gain
  mtx_type[m][m] P; // m x m the covariance of the estimation process at step k
  mtx_type[m][m] Q; // m x m covariance of the process noise (e.g. noise from the motors)
  mtx_type[n][n] R; // n x n matrix containing covariance between each of the sensors
  
  mtx_type[m][m] I; // m x m identity matrix

};

typedef struct kalman_filter KalmanFilter;

KalmanFilter kf;

void setupFilter(){
  
  // Set up I
  for(int j = 0; j < m; j++){
    for(int k = 0; k < m; k++){
      if(j == k){
        kf.I[j][k] = 1;
      } else {
        kf.I[j][k] = 0;
      }
    }
  }
  
}

void updateFilter(){
  // Time Update
  mtx_type x_minus[kf.m][1];
  mtx_type temp1_mx1[kf.m][1];
  mtx_type temp2_mx1[kf.m][1];
  Matrix.Multiply((mtx_type*)kf.A, (mtx_type*)kf.x, kf.m, kf.m, 1, (mtx_type*)temp1_mx1); //Ax
  Matrix.Multiply((mtx_type*)kf.B, (mtx_type*)kf.u, kf.m, kf.p, 1, (mtx_type*)temp2_mx1); //Bu
  Matrix.Add((mtx_type*) temp1_mx1, (mtx_type*) temp2_mx1, kf.m, 1, (mtx_type*) x_minus); // x^- = Ax + Bu

  mtx_type A_transpose[kf.m][kf.m];
  mtx_type P_minus[kf.m][kf.m];
  mtx_type temp3_mxm[kf.m][kf.m];
  mtx_type temp4_mxm[kf.m][kf.m];
  Matrix.Transpose((mtx_type*) kf.A, kf.m, kf.m, (mtx_type*) A_transpose); //A^T
  Matrix.Multiply((mtx_type*)kf.A, (mtx_type*)kf.P, kf.m, kf.m, kf.m, (mtx_type*)temp3_mxm); //AP
  Matrix.Multiply((mtx_type*)temp3_mxm, (mtx_type*)A_transpose, kf.m, kf.m, kf.m, (mtx_type*)temp4_mxm); //APA^T
  Matrix.Add((mtx_type*) temp4_mxm, (mtx_type*) kf.Q, kf.m, kf.m, (mtx_type*) P_minus); //P^- = APA^T + Q

  // Measurement Update
  mtx_type temp5_nxm[kf.n][kf.m];
  mtx_type temp6_nxn[kf.n][kf.n];
  mtx_type temp7_nxn[kf.n][kf.n];
  mtx_type temp8_mxn[kf.m][kf.n];
  mtx_type C_transpose[kf.m][kf.n];
  Matrix.Transpose((mtx_type*) kf.C, kf.n, kf.m, (mtx_type*) C_transpose); //C^T
  Matrix.Multiply((mtx_type*)kf.C, (mtx_type*)P_minus, kf.n, kf.m, kf.m, (mtx_type*)temp5_nxm); //CP
  Matrix.Multiply((mtx_type*)temp5_nxm, (mtx_type*)C_transpose, kf.n, kf.m, kf.n, (mtx_type*)temp6_nxn); //CPC^T
  Matrix.Add((mtx_type*) temp6_nxn, (mtx_type*) kf.R, kf.n, kf.n, (mtx_type*) temp7_nxn); //CPC^T + R
  Matrix.Invert((mtx_type*)temp7_nxn, kf.n);//(CPC^T + R)^-1
  Matrix.Multiply((mtx_type*C_transpose, (mtx_type*)temp7_nxn, kf.m, kf.n, kf.n, (mtx_type*)temp8_mxn); //C^T(CPC^T + R)^-1
  Matrix.Multiply((mtx_type*)P_minus, (mtx_type*)temp8_mxn, kf.m, kf.m, kf.n, (mtx_type*)kf.K); //K = C^T(CPC^T + R)^-1

  mtx_type temp9_nx1[kf.n][1];
  mtx_type temp10_nx1[kf.n][1];
  Matrix.Multiply((mtx_type*)kf.C, (mtx_type*)kf.x, kf.n, kf.m, 1, (mtx_type*) temp9_nx1); //Cx
  Matrix.Subtract((mtx_type*) kf.z, (mtx_type*) temp9_nx1, kf.n, 1, (mtx_type*) temp10_nx1); //z-Cx
  Matrix.Multiply((mtx_type*)kf.K, (mtx_type*)temp10_nx1, kf.m, kf.n, 1, (mtx_type*)temp1_mx1); //K(z-Cx)
  Matrix.Add((mtx_type*)kf.x, (mtx_type*) temp1_mx1, kf.m, 1, (mtx_type*) temp2_mx1); //x + K(z-Cx)
  Matrix.Copy((mtx_type*)temp2_mx1, (mtx_type*) temp1_mx1, kf.m, 1, (mtx_type*) kf.x); //x = x + K(z-Cx)
  

  Matrix.Multiply((mtx_type*)kf.K, (mtx_type*)kf.C, kf.m, kf.n, kf.m, (mtx_type*)temp3_mxm); //KC
  Matrix.Subtract((mtx_type*) kf.I, (mtx_type*) temp3_mxm, kf.m, kf.m, (mtx_type*) temp4_mxm); //I-KC
  Matrix.Multiply((mtx_type*)ktemp10_nx1, (mtx_type*)P_minus, kf.m, kf.m, kf.m, (mtx_type*)kf.P); //P=(I-KC)P^-
}
