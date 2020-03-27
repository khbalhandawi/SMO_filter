//=============================================================================================
// MahonyAHRS.h
//=============================================================================================
//
// Madgwick's implementation of Mayhony's AHRS algorithm.
// See: http://www.x-io.co.uk/open-source-imu-and-ahrs-algorithms/
//
// Date			Author			Notes
// 29/09/2011	SOH Madgwick    Initial release
// 02/10/2011	SOH Madgwick	Optimised for reduced CPU load
//
//=============================================================================================
#ifndef MahonyAHRS_h
#define MahonyAHRS_h
#include <math.h>
#include <Eigen/Core>
#include <Eigen/Geometry>

typedef Eigen::Matrix3d Matrix_3;
typedef Eigen::Matrix<double, 3, 1> Vector_3;
//--------------------------------------------------------------------------------------------
// Variable declaration

class Mahony {
private:
	double twoKp;		// 2 * proportional gain (Kp)
	double twoKi;		// 2 * integral gain (Ki)
	double q0, q1, q2, q3;	// quaternion of sensor frame relative to auxiliary frame
	Matrix_3 R; // rotation matrix
	Eigen::Quaterniond q; // quaternion object
	double integralFBx, integralFBy, integralFBz;  // integral error terms scaled by Ki
	double invSampleFreq;
	double roll, pitch, yaw;
	char anglesComputed;
	double invSqrt(double x);
	void computeAngles();

//-------------------------------------------------------------------------------------------
// Function declarations

public:
	Mahony();
	void begin(double Kp, double Ki, double sampleFrequency) {
		twoKp = 2.0 * Kp;
		twoKi = 2.0 * Ki;
		invSampleFreq = 1.0 / sampleFrequency; 
	}
	void update(double gx, double gy, double gz, double ax, double ay, double az, double mx, double my, double mz);
	void updateIMU(double gx, double gy, double gz, double ax, double ay, double az);
	double getRoll() {
		if (!anglesComputed) computeAngles();
		return roll * 57.29578;
	}
	double getPitch() {
		if (!anglesComputed) computeAngles();
		return pitch * 57.29578;
	}
	double getYaw() {
		if (!anglesComputed) computeAngles();
		return yaw * 57.29578;
	}
	double getRollRadians() {
		if (!anglesComputed) computeAngles();
		return roll;
	}
	double getPitchRadians() {
		if (!anglesComputed) computeAngles();
		return pitch;
	}
	double getYawRadians() {
		if (!anglesComputed) computeAngles();
		return yaw;
	}
	double getQuaternion() {
		return 	q0, q1, q2, q3;
	}
	double getKp() {
		return twoKp;
	}
	double getKi() {
		return twoKi;
	}
	double getFreq() {
		return invSampleFreq;
	}

	Matrix_3 getRotationmatrix() {
		return 	R;
	}
};

#endif
