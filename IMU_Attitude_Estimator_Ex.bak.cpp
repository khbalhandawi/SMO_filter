#include <iostream>
#include <vector>
#include "time.h"
#include <Eigen/Core>
#include <Eigen/Geometry>
#include <Eigen/LU>

#include "Convert.h"
#include "EKF_Attitude.h"
#include "Mahony_Attitude.h"
#include "ESKF_Attitude.h"
#include "MahonyAHRS.h"

#include "matplotlibcpp.h"
#include "windows.h"

#define _USE_MATH_DEFINES
#include "math.h"
#define _CRT_SECURE_NO_WARNINGS

using namespace std;
using namespace IMU;

namespace plt = matplotlibcpp;

typedef struct {
	double Kp, Ki;
} mahony_gains;

void plot_data(std::vector<double> Index, std::vector<double> x_plot, 
	std::vector<double> y_plot, std::vector<double> z_plot, std::string plot_title, int total_units)
{
	plt::figure_size(1200, 500);
	plt::title("plot_title");
	plt::subplot(1, 3, 1);
	plt::named_plot("x", Index, x_plot, "r");

	plt::xlim(0, total_units);
	plt::title("x");
	plt::legend();

	plt::subplot(1, 3, 2);
	plt::named_plot("y", Index, y_plot, "b");

	plt::xlim(0, total_units);
	plt::title("y");
	plt::legend();

	plt::subplot(1, 3, 3);
	plt::named_plot("z", Index, z_plot, "g");

	plt::xlim(0, total_units);
	plt::title("z");
	plt::legend();
}

void convert_measurement_data(Eigen::MatrixXd measurement, vector<double> &measurement_x,
	vector<double> &measurement_y, vector<double> &measurement_z)
{
	unsigned int i = 0;

	do{
		measurement_x.push_back(measurement.row(i)[0]);
		measurement_y.push_back(measurement.row(i)[1]);
		measurement_z.push_back(measurement.row(i)[2]);

		i++;
	} while (i < measurement.rows());
}

int main(int argc, char **argv)
{
	float Kp, Ki, freq;
	int total_units;
	mahony_gains gains_estimator;

	printf("\r\nEstimator Gains (Kp, Ki) : \n");
	scanf_s("%f %f", &Kp, &Ki);
	printf("\r\nIMU frequency (f) : \n");
	scanf_s("%f", &freq);
	printf("\r\nnumber of units : \n");
	scanf_s("%i", &total_units);
	gains_estimator.Kp = Kp;
	gains_estimator.Ki = Ki;
	printf("%f %f\n", gains_estimator.Kp, gains_estimator.Ki);

	Eigen::MatrixXd data, measurements, groundtruth, acceleration, gyroscope, magnetometer;
	data = IMU::readFromfile("./datasets/NAV3_data.bin");
	if (data.isZero())
		return 0;

	const int Rows = data.rows() - 1;
	measurements = data.block(0, 0, Rows, 9);
	groundtruth = data.block(0, 9, Rows, 3) * 180 / M_PI;
	acceleration = data.block(0, 0, Rows, 3);
	gyroscope = data.block(0, 3, Rows, 3);
	magnetometer = data.block(0, 6, Rows, 3);

	IMU::EKF_Attitude EKF_AHRS(true, 1 / freq);
	//IMU::Mahony_Attitude Mahony_att(Eigen::Vector2d(Kp, Ki), 1 / freq);
	Mahony mahony_ahrs;

	Eigen::Matrix<double, 12, 1> ESKF_InitVec;
	ESKF_InitVec << 1e-5*Eigen::Vector3d::Ones(), 1e-9*Eigen::Vector3d::Ones(),
		1e-3*Eigen::Vector3d::Ones(), 1e-4*Eigen::Vector3d::Ones();
	IMU::ESKF_Attitude ESKF_AHRS(ESKF_InitVec, 1 / freq);

	unsigned int i = 0;
	Eigen::MatrixXd Euler(measurements.rows(), 3), Euler1(measurements.rows(), 3), Euler2(measurements.rows(), 3);

	std::vector<double> Index, Roll, Pitch, Yaw, Roll_gt, Pitch_gt, Yaw_gt;
	std::vector<double> Roll1, Pitch1, Yaw1, Roll_gt1, Pitch_gt1, Yaw_gt1;
	std::vector<double> Roll2, Pitch2, Yaw2, Roll_gt2, Pitch_gt2, Yaw_gt2;
	TicToc tc;
	do
	{
		Eigen::MatrixXd measure;
		Eigen::Quaterniond quaternion;

		Vector_3 Euler_single;

		quaternion = EKF_AHRS.Run(measurements.row(i).transpose());
		Euler.row(i) = Quaternion_to_Euler(quaternion).transpose();

		Vector_9 measurement_row = measurements.row(i).transpose();
		float gx = measurement_row[3];
		float gy = measurement_row[4];
		float gz = measurement_row[5];
		float ax = measurement_row[0];
		float ay = measurement_row[1];
		float az = measurement_row[2];
		mahony_ahrs.updateIMU(gx, gy, gz, ax, ay, az);
		Vector_3 angles;
		angles(0) = mahony_ahrs.getRoll();
		angles(1) = mahony_ahrs.getPitch();
		angles(2) = mahony_ahrs.getYaw();
		Euler1.row(i) = angles;

		//quaternion = Mahony_att.Run(measurements.row(i).transpose());
		// Euler1.row(i) = Quaternion_to_Euler(quaternion).transpose();

		quaternion = ESKF_AHRS.Run(measurements.row(i).transpose());
		Euler2.row(i) = Quaternion_to_Euler(quaternion).transpose();


		Index.push_back(i*1.0);
		Roll.push_back(Euler.row(i)[0]);
		Pitch.push_back(Euler.row(i)[1]);
		Yaw.push_back(Euler.row(i)[2]);

		Roll1.push_back(Euler1.row(i)[0]);
		Pitch1.push_back(Euler1.row(i)[1]);
		Yaw1.push_back(Euler1.row(i)[2]);

		Roll2.push_back(Euler2.row(i)[0]);
		Pitch2.push_back(Euler2.row(i)[1]);
		Yaw2.push_back(Euler2.row(i)[2]);

		Roll_gt.push_back(groundtruth.row(i)[0]);
		Pitch_gt.push_back(groundtruth.row(i)[1]);
		Yaw_gt.push_back(groundtruth.row(i)[2]);

		i++;
	} while (i < measurements.rows());

	cout << tc.toc() << "ms" << endl;
	writeTofile(Euler, "Euler.bin");

	// EKF estimation
	plt::figure_size(1200, 500);
	plt::title("Roll pitch yaw comparision");
	plt::subplot(1, 3, 1);
	plt::named_plot("EKF", Index, Roll, "b");
	plt::named_plot("ESKF", Index, Roll2, "g");
	plt::named_plot("Groundtruth", Index, Roll_gt, "r");

	plt::xlim(0, total_units);
	plt::title("Roll");
	plt::legend();

	plt::subplot(1, 3, 2);
	plt::named_plot("EKF", Index, Pitch, "b");
	plt::named_plot("ESKF", Index, Pitch2, "g");
	plt::named_plot("Groundtruth", Index, Pitch_gt, "r");

	plt::xlim(0, total_units);
	plt::title("Pitch");
	plt::legend();

	plt::subplot(1, 3, 3);
	plt::named_plot("EKF", Index, Yaw, "b");
	plt::named_plot("ESKF", Index, Yaw2, "g");
	plt::named_plot("Groundtruth", Index, Yaw_gt, "r");

	plt::xlim(0, total_units);
	plt::title("Yaw");
	plt::legend();

	// Mahony estimation
	plt::figure_size(1200, 500);
	plt::title("Roll pitch yaw comparision");
	plt::subplot(1, 3, 1);
	plt::named_plot("Mahony", Index, Roll1, "k");
	plt::named_plot("Groundtruth", Index, Roll_gt, "r");

	plt::xlim(0, total_units);
	plt::title("Roll");
	plt::legend();

	plt::subplot(1, 3, 2);
	plt::named_plot("Mahony", Index, Pitch1, "k");
	plt::named_plot("Groundtruth", Index, Pitch_gt, "r");

	plt::xlim(0, total_units);
	plt::title("Pitch");
	plt::legend();

	plt::subplot(1, 3, 3);
	plt::named_plot("Mahony", Index, Yaw1, "k");
	plt::named_plot("Groundtruth", Index, Yaw_gt, "r");

	plt::xlim(0, total_units);
	plt::title("Yaw");
	plt::legend();

	std::vector<double> acc_x_plot, acc_y_plot, acc_z_plot;
	convert_measurement_data(acceleration, acc_x_plot, acc_y_plot, acc_z_plot);
	std::string s = "Acc";
	plot_data(Index, acc_x_plot, acc_y_plot, acc_z_plot, s, total_units);

	std::vector<double> gyr_x_plot, gyr_y_plot, gyr_z_plot;
	convert_measurement_data(gyroscope, gyr_x_plot, gyr_y_plot, gyr_z_plot);
	s = "Gyro";
	plot_data(Index, gyr_x_plot, gyr_y_plot, gyr_z_plot, s, total_units);

	std::vector<double> mag_x_plot, mag_y_plot, mag_z_plot;
	convert_measurement_data(magnetometer, mag_x_plot, mag_y_plot, mag_z_plot);
	s = "Mag";
	plot_data(Index, mag_x_plot, mag_y_plot, mag_z_plot, s, total_units);

	plt::show();

	return 0;
}

