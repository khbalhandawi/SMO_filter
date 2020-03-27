#include <iostream>
#include <vector>
#include "time.h"
#include <Eigen/Core>
#include <Eigen/Geometry>
#include <Eigen/LU>
#include "Convert.h"
#include "Mahony_Attitude.h"
#include "MahonyAHRS.h"
#include "matplotlibcpp.h"
#include "windows.h"
#include <iomanip>

#define _USE_MATH_DEFINES
#include "math.h"
#define PI 3.14159265

using namespace std;
using namespace IMU;

namespace plt = matplotlibcpp;




Eigen::MatrixXd mahony_filter(Mahony *mahony_ahrs, Vector_3 acc_data, Vector_3 gyro_data)
{

	Vector_3 angles, Euler_out, gravity_b, acc_b, acc_i_out;
	Matrix_3 Rotation_matrix;
	double g = 9.81;
	Vector_3 gravity_i;
	gravity_i << 0.0, 0.0, g;

	// Use mahony estimated attitude
	double gx = gyro_data(0); double gy = gyro_data(1); double gz = gyro_data(2);
	double ax = acc_data(0); double ay = acc_data(1); double az = acc_data(2);

	//cout << mahony_ahrs.getKp() << '|' << mahony_ahrs.getKi() << '|'  << mahony_ahrs.getFreq() << '\n';

	mahony_ahrs->updateIMU(gx, gy, gz, ax, ay, az);
	angles(0) = mahony_ahrs->getRoll()  ;
	angles(1) = mahony_ahrs->getPitch() ;
	angles(2) = mahony_ahrs->getYaw()   ;
	Euler_out = angles;
	Rotation_matrix = mahony_ahrs->getRotationmatrix(); // Rotation matrix from body to interial R_b^I

	gravity_b = Rotation_matrix.inverse() * gravity_i; // gravity in the body frame using R_I^b
	acc_b = acc_data * g - gravity_b; // acceleraion in the body frame ( - gravity)
	acc_i_out = Rotation_matrix * acc_b; // acceleraion in the interial frame

	Eigen::MatrixXd output(2,3);
	output.row(0) = acc_i_out;
	output.row(1) = Euler_out;

	return output;
}

Eigen::MatrixXd chib_filter(Vector_3 pos_est_current, Vector_3 pos_true, Vector_3 acc_i_raw, 
	Vector_3 acc_est_prev, Vector_3 st_vel_prev, Vector_3 integ_vel, Vector_3 integ_pos, 
	double  c, double  dt, double  lambda, double  q, double p)
{

	Vector_3 acc_est_out, vel_est_out, super_twist_vel_out, pos_est_out;
	Matrix_3 Rotation_matrix;

	Vector_3 pos_err, pos_err_neg, relay, acc_est, acc_est_current,
		vel_inc, vel_est, super_twist, super_twist_vel, st_vel_current, pos_inc, pos_est;
	
	// Acceleration calculations
	
	pos_err = pos_est_current - pos_true;
	
	pos_err_neg = pos_err * -1;

	relay = pos_err_neg.cwiseSign(); 

	acc_est = acc_i_raw + (c * relay);

	acc_est_current = acc_est;

		
	// Velocity calculations
	vel_inc = (acc_est_prev + acc_est_current) / 2.0;
	integ_vel = (vel_inc * dt) + integ_vel; // integrate acceleration to get velocity
	vel_est = integ_vel;

	// Lipshitchz discontinous function
	super_twist = pos_err_neg.cwiseAbs().array().pow(q / p) * pos_err_neg.array().sign(); // CHECK SIGNUM FUNTION -1,0,1 or 0,1
	// Linear function
	//super_twist = pos_err_neg;
	super_twist_vel = vel_est + (lambda * super_twist);
	st_vel_current = super_twist_vel;

	// Position calculations
	pos_inc = (st_vel_prev + st_vel_current) / 2.0;
	integ_pos = (pos_inc * dt) + integ_pos; // integrate acceleration to get velocity
	pos_est = integ_pos;

	acc_est_out = acc_est;
	vel_est_out = vel_est;
	super_twist_vel_out = super_twist_vel;
	pos_est_out = pos_est;

	Eigen::MatrixXd output(6, 3);
	output.row(0) = acc_est_out;
	output.row(1) = vel_est_out;
	output.row(2) = super_twist_vel_out;
	output.row(3) = pos_est_out;
	output.row(4) = integ_vel;
	output.row(5) = integ_pos;

	return output;
	
}



int main(int argc, char* argv[])
{

	double c = 0., lambda = 0.;
	double Kp = 2.0, Ki = 1.0;
	double freq = 1000.0;
	double g = 9.81;
	double pos_err_sq_sum = 0.0, pos_err_RMS;
	double p = 2.0, q = 1.0;

	c = atof(argv[1]);
	lambda = atof(argv[2]);

	cout << "c: " << c << " | lambda: " << lambda << " | Kp: " << Kp << " | Ki: " << Ki << "\n";

	Vector_3 pos_current, vel_current, acc_current,
		pos_est_current, vel_prev, acc_prev, st_vel_prev, st_vel_current,
		pos_true, velocity_true, acc_true, pos_inc, vel_inc, pos_err_neg, pos_err, pos_err_true, pos_err_sq, relay, super_twist;
	double time_prev;

	Eigen::MatrixXd data, measurements, groundtruth_ds, groundtruth_pos_ds,
		groundtruth_true, groundtruth_pos_true, acceleration,
		gyroscope, magnetometer, time_raw;

	Vector_3 integ_vel, integ_pos;

	integ_vel << 0.0, 0.0, 0.0; // initialize integrator
	integ_pos << 0.0, 0.0, 0.0; // initialize integrator

	std::vector<double> Index, Roll, Pitch, Yaw, Roll_gt, Pitch_gt, Yaw_gt,
		acc_i_x, acc_i_y, acc_i_z, vel_est_x, vel_est_y, vel_est_z,
		pos_est_x, pos_est_y, pos_est_z, pos_true_x, pos_true_y, pos_true_z,
		pos_ds_x, pos_ds_y, pos_ds_z, time, time_plot;


	data = IMU::readFromfile("./datasets/NAV3_data.bin");
	if (data.isZero())
		return 0;

	const int Rows = data.rows() - 1;
	measurements = data.block(0, 0, Rows, 22);
	time_raw = data.block(0, 21, Rows, 1);
	groundtruth_ds = data.block(0, 9, Rows, 3) * 180 / M_PI;
	groundtruth_pos_ds = data.block(0, 12, Rows, 3);
	groundtruth_true = data.block(0, 15, Rows, 3) * 180 / M_PI;
	groundtruth_pos_true = data.block(0, 18, Rows, 3);
	acceleration = data.block(0, 0, Rows, 3);
	gyroscope = data.block(0, 3, Rows, 3);
	magnetometer = data.block(0, 6, Rows, 3);

	Eigen::MatrixXd Euler(measurements.rows(), 3);
	Eigen::MatrixXd acc_i(measurements.rows(), 3);
	Eigen::MatrixXd pos_est(measurements.rows(), 3);
	Eigen::MatrixXd vel_est(measurements.rows(), 3);
	Eigen::MatrixXd acc_est(measurements.rows(), 3);
	Eigen::MatrixXd super_twist_vel(measurements.rows(), 3);
	Eigen::MatrixXd mahony_out;
	Eigen::MatrixXd chib_out, chib_out_j;


	Vector_3 acc_data, gyro_data;
	Vector_3 acc_i_raw;

	Mahony mahony_ahrs;
	mahony_ahrs.begin(Kp, Ki, freq);
	//double dt = (time_raw.row(i)[0] - time_prev);
	double dt = 0.001;
	unsigned int i = 0;


	//initialize pos vector
	pos_est_current = groundtruth_pos_ds.row(i).transpose();

	//initialize acc vector
	acc_data = acceleration.row(i).transpose();//changed trans
	gyro_data = gyroscope.row(i).transpose();//changed trans

	mahony_out = mahony_filter(&mahony_ahrs, acc_data, gyro_data);
	acc_i.row(i) = mahony_out.block<1, 3>(0, 0);
	Euler.row(i) = mahony_out.block<1, 3>(1, 0);
	
	acc_est.row(i) = acc_i.row(i);
	Vector_3 acc_est_prev = acc_est.row(i).transpose();

	// initialize super_twist_vel vector
	vel_inc = acc_est.row(i).transpose();
	super_twist_vel.row(i) = (vel_inc * dt);
	st_vel_prev = super_twist_vel.row(i);// get last row of the matrix
	
	
	unsigned int n_memory = 27;
	unsigned int j = 0;

	// On board IMU loop
	TicToc tc;
	do
	{
		 acc_data = acceleration.row(i).transpose();
		 gyro_data = gyroscope.row(i).transpose();
	
		

		// initialize time vector
		if (i == 0) { time_prev = 0.0; }
		else { time_prev = time_raw.row(i - 1)[0]; }

		//  Mahony section
		mahony_out = mahony_filter(&mahony_ahrs, acc_data, gyro_data);

		acc_i.row(i) = mahony_out.block<1, 3>(0, 0);
		Euler.row(i) = mahony_out.block<1, 3>(1, 0);

		//cout << setprecision(14) << "x, y, z: " << Euler.row(i)(0) << "| " << Euler.row(i)(1) << "| " << Euler.row(i)(2) << endl;

		// CH Boiko section

		acc_i_raw = acc_i.row(i).transpose();
		pos_true = groundtruth_pos_ds.row(i).transpose();

		chib_out = chib_filter( pos_est_current,  pos_true,  acc_i_raw,
			 acc_est_prev,  st_vel_prev,  integ_vel,  integ_pos,
			  c,   dt,   lambda,   q,  p);

		acc_est.row(i) = chib_out.block<1, 3>(0, 0);
		vel_est.row(i) = chib_out.block<1, 3>(1, 0);
		super_twist_vel.row(i) = chib_out.block<1, 3>(2, 0);
		pos_est.row(i) = chib_out.block<1, 3>(3, 0);
		integ_vel = chib_out.block<1, 3>(4, 0); 
		integ_pos = chib_out.block<1, 3>(5, 0); 

		//cout << "rem:" << (i + 1) % n_memory << "\n";

		if ( (i+1) % n_memory == 0)
		{
			
			//i: 0 1 2 3 4 5 6 7
			//i + 1 : 1 2 3 4 5 6 7 8
			//n_memory : = 4;
			//i + 1 % n_memory: 0 0 0 1 0 0 0 1
			
			if ( (i - n_memory) == -1)
			{ 
				j = 0;
			}
			else
			{
				j = i - n_memory;
			}

			//cout << "j:" << j << "\n";
			//cout << "i:" << i << "\n";

			//Reset the integrators
			integ_vel = vel_est.row(j).transpose();
			integ_pos = pos_est.row(j).transpose();

			//initialize for next loop iteration
			//pos vector
			pos_est_current = pos_est.row(j).transpose();
			//pos_est will be updated at end of loop!
			//acc vector
			acc_est_prev = acc_est.row(j).transpose();
			//super_twist_vel vector
			st_vel_prev = super_twist_vel.row(j).transpose();

			while (j <= i) {
				pos_true = groundtruth_pos_ds.row(i).transpose();

				acc_i_raw = acc_i.row(j).transpose();

				chib_out_j = chib_filter(pos_est_current, pos_true, acc_i_raw,
					acc_est_prev, st_vel_prev, integ_vel, integ_pos,
					c, dt, lambda, q, p);

				acc_est.row(j) = chib_out_j.block<1, 3>(0, 0);
				vel_est.row(j) = chib_out_j.block<1, 3>(1, 0);
				super_twist_vel.row(j) = chib_out_j.block<1, 3>(2, 0);
				pos_est.row(j) = chib_out_j.block<1, 3>(3, 0);
				integ_vel = chib_out_j.block<1, 3>(4, 0);
				integ_pos = chib_out_j.block<1, 3>(5, 0);

				cout << setprecision(14) << "x, y, z: " << i << "| " << j << "| " << integ_pos(0) << "| " << integ_pos(1) << "| " << integ_pos(2) << endl;


				j++;

				// initialize for next loop iteration
				// pos vector
				pos_est_current = pos_est.row(j-1).transpose(); //get last row of the matrix
				// pos_est will be updated at end of loop!
				//acc vector
				acc_est_prev = acc_est.row(j-1).transpose();
				// super_twist_vel vector
				st_vel_prev = super_twist_vel.row(j-1).transpose();// get last row of the matrix
			}

		}
		
		Index.push_back(i * 1.0);
		Roll.push_back(Euler.row(i)[0]);
		Pitch.push_back(Euler.row(i)[1]);
		Yaw.push_back(Euler.row(i)[2]);

		Roll_gt.push_back(groundtruth_true.row(i)[0]);
		Pitch_gt.push_back(groundtruth_true.row(i)[1]);
		Yaw_gt.push_back(groundtruth_true.row(i)[2]);

		acc_i_x.push_back(acc_i.row(i)[0]);
		acc_i_y.push_back(acc_i.row(i)[1]);
		acc_i_z.push_back(acc_i.row(i)[2]);

		vel_est_x.push_back(super_twist_vel.row(i)[0]);
		vel_est_y.push_back(super_twist_vel.row(i)[1]);
		vel_est_z.push_back(super_twist_vel.row(i)[2]);

		time.push_back(time_raw.row(i)[0]);

		pos_est_x.push_back(pos_est.row(i)[0]);
		pos_est_y.push_back(pos_est.row(i)[1]);
		pos_est_z.push_back(pos_est.row(i)[2]);

		pos_true_x.push_back(groundtruth_pos_true.row(i)[0]);
		pos_true_y.push_back(groundtruth_pos_true.row(i)[1]);
		pos_true_z.push_back(groundtruth_pos_true.row(i)[2]);

		pos_ds_x.push_back(groundtruth_pos_ds.row(i)[0]);
		pos_ds_y.push_back(groundtruth_pos_ds.row(i)[1]);
		pos_ds_z.push_back(groundtruth_pos_ds.row(i)[2]);

		// Compute RMS error with respect to true position
		pos_err_true = pos_est_current - groundtruth_pos_true.row(i).transpose();

		pos_err_sq = pos_err_true.array().pow(2);
		pos_err_sq_sum += pos_err_sq.array().sum();
		
		i++;

		// initialize for next loop iteration
		// pos vector
		pos_est_current = pos_est.row(i - 1).transpose(); //get last row of the matrix
		// pos_est will be updated at end of loop!
		//acc vector
		acc_est_prev = acc_est.row(i - 1).transpose();
		// super_twist_vel vector
		st_vel_prev = super_twist_vel.row(i - 1).transpose();// get last row of the matrix
	

	} while (i < measurements.rows());

	// acceleration in the intertial frame estimation
	plt::figure_size(1200, 500);

	plt::subplot(1, 3, 1);
	plt::named_plot("truth", time, pos_true_x, "-k");
	plt::named_plot("estimated", time, pos_est_x, "--r");
	plt::named_plot("downsampled", time, pos_ds_x, "-g");
	plt::legend();
	plt::title("x");

	plt::subplot(1, 3, 2);
	plt::named_plot("truth", time, pos_true_y, "-k");
	plt::named_plot("estimated", time, pos_est_y, "--r");
	plt::named_plot("downsampled", time, pos_ds_y, "-g");
	plt::legend();
	plt::title("y");

	plt::subplot(1, 3, 3);
	plt::named_plot("truth", time, pos_true_z, "-k");
	plt::named_plot("estimated", time, pos_est_z, "--r");
	plt::named_plot("downsampled", time, pos_ds_z, "-g");
	plt::legend();
	plt::title("z");
	// save figure
	const char* filename_pos = "./position_plot.pdf";
	std::cout << "saving result to " << filename_pos << std::endl;;
	plt::save(filename_pos);

	// acceleration in the intertial frame estimation
	plt::figure_size(1200, 500);

	plt::subplot(1, 3, 1);
	plt::named_plot("estimated", time, vel_est_x, "-r");
	plt::legend();
	plt::title("x");

	plt::subplot(1, 3, 2);
	plt::named_plot("estimated", time, vel_est_y, "-r");
	plt::legend();
	plt::title("y");

	plt::subplot(1, 3, 3);
	plt::named_plot("estimated", time, vel_est_z, "-r");
	plt::legend();
	plt::title("z");
	// save figure
	const char* filename_vel = "./velocity_plot.pdf";
	std::cout << "saving result to " << filename_vel << std::endl;;
	plt::save(filename_vel);

	// Mahony estimation
	plt::figure_size(1200, 500);
	plt::title("Roll pitch yaw comparision");
	plt::subplot(1, 3, 1);
	plt::named_plot("Mahony", time, Roll, "k");
	plt::named_plot("Groundtruth", time, Roll_gt, "r");

	//plt::xlim(0, total_units);
	plt::title("Roll");
	plt::legend();

	plt::subplot(1, 3, 2);
	plt::named_plot("Mahony", time, Pitch, "k");
	plt::named_plot("Groundtruth", time, Pitch_gt, "r");

	//plt::xlim(0, total_units);
	plt::title("Pitch");
	plt::legend();

	plt::subplot(1, 3, 3);
	plt::named_plot("Mahony", time, Yaw, "k");
	plt::named_plot("Groundtruth", time, Yaw_gt, "r");

	//plt::xlim(0, total_units);
	plt::title("Yaw");
	plt::legend();
	plt::show();
	// save figure
	const char* filename_att = "./attitude_plot.pdf";
	std::cout << "saving result to " << filename_att << std::endl;;
	plt::save(filename_att);

	pos_err_RMS = sqrt(pos_err_sq_sum / (i + 1.0));

	ofstream out("POS_ERR_RMS.txt");
	out.precision(9); // number of decimal places to output
	/*cout << tc.toc() << "ms" << endl;*/
	out << fixed << pos_err_RMS << endl;
	out.close();

	cout << "RMS Error: " << pos_err_RMS << endl;
	
	return 0;
}

//// Optitrack loop
//i = 0;
//std::vector<double> vel_x_true, vel_y_true, vel_z_true,
//	acc_i_x_true, acc_i_y_true, acc_i_z_true, time_optitrack_plot;
//do
//{
//	time_optitrack_plot.push_back(time_optitrack.row(i)[0]);
//	vel_x_true.push_back(groundtruth_vel.row(i)[0]);
//	vel_y_true.push_back(groundtruth_vel.row(i)[1]);
//	vel_z_true.push_back(groundtruth_vel.row(i)[2]);

//	acc_i_x_true.push_back(groundtruth_acc.row(i)[0]);
//	acc_i_y_true.push_back(groundtruth_acc.row(i)[1]);
//	acc_i_z_true.push_back(groundtruth_acc.row(i)[2]);
//	
//	i++;
//} while (i < measurements_opti.rows());


//initialize time vector
//if (i == 0) { time_prev = 0.0f; }
//else { time_prev = time.back(); }
//
//float dt = (time_raw.row(i)[0] - time_prev);
//
//// initialize velocity vector
//if (i == 0) { pos_prev = groundtruth_pos.row(i); } // get last row of the matrix
//else { pos_prev = groundtruth_pos.row(i - 1); }
//
//pos_current = groundtruth_pos.row(i);
//velocity_true = (pos_current - pos_prev) / dt;
//vel_x_true.push_back(velocity_true[0]);
//vel_y_true.push_back(velocity_true[1]);
//vel_z_true.push_back(velocity_true[2]);
//
//// initialize acc vector
//if (i == 0) { vel_prev << vel_x_true[i], vel_y_true[i], vel_z_true[i]; } // get last row of the matrix
//else { vel_prev << vel_x_true[i - 1], vel_y_true[i - 1], vel_z_true[i - 1]; } // get last row of the matrix
//
//vel_current << vel_x_true[i], vel_y_true[i], vel_z_true[i]; // get last row of the matrix
//acc_true = (vel_current - vel_prev) / dt;
//
//acc_i_x_true.push_back(acc_true[0] / 100000.0f);
//acc_i_y_true.push_back(acc_true[1] / 100000.0f);
//acc_i_z_true.push_back(acc_true[2] / 100000.0f);