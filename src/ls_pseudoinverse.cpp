#include "ls_pseudoinverse.hpp"

namespace ls_pseudoinverse
{

	LS_Pseudo::LS_Pseudo(adap_parameters_estimator::DOFS _dof)
	{
		dof = _dof;
	}

	LS_Pseudo::~LS_Pseudo()
	{
	}

/*	void LS_Pseudo::convert_states(std::queue<base::samples::Joints> &queueOfForces, std::queue<base::samples::RigidBodyState> &queueOfRBS, base::MatrixXd &states)
	{
		if (queueOfForces.size() != queueOfRBS.size())
			std::cout << std::endl << "Vector of forces and velocities should have the same size "<< std::endl << std::endl;
		else
		{
			states.resize(4,queueOfRBS.size());

			base::samples::Joints force;
			base::samples::RigidBodyState vel;

			for (int i=0; i<queueOfRBS.size(); i++)
			{
				force	= queueOfForces.front();
				vel		= queueOfRBS.front();

				queueOfForces.pop();
				queueOfRBS.pop();

				states(0,i) = force.elements[dof].effort;
				states(3,i) = 1;
				if (dof < 3)
				{
					states(1,i) = vel.velocity[dof] * fabs(double(vel.velocity[dof]));
					states(2,i) = vel.velocity[dof];
				}
				else
				{
					states(1,i) = vel.angular_velocity[dof] * fabs(double(vel.angular_velocity[dof]));
					states(2,i) = vel.angular_velocity[dof];
				}

				queueOfForces.push(force);
				queueOfRBS.push(vel);
			}
		}
	}
	

	void LS_Pseudo::convert_acceleration(std::queue<base::samples::RigidBodyAcceleration> &queueOfRBA, base::MatrixXd &acceleration)
	{
		acceleration.resize(1,queueOfRBA.size());

		base::samples::RigidBodyAcceleration acc;

		for( int i=0; i < queueOfRBA.size(); i++)
		{
			acc = queueOfRBA.front();

			if(dof<3)
				acceleration(0,i) = acc.acceleration(dof);

			queueOfRBA.push(acc);
		}
	}
*/

	void LS_Pseudo::convert_data(std::queue<adap_samples_input::DynamicAUV> &queueOfDyn, base::MatrixXd &states, base::MatrixXd &acceleration)
	{	// [acceleration]1xn; [states]4xn
		states.resize(4,queueOfDyn.size());
		acceleration.resize(1,queueOfDyn.size());
		adap_samples_input::DynamicAUV dyn;

		// // fist acceleration is nan
		dyn = queueOfDyn.front();
		queueOfDyn.pop();
		if((isnan)((double)dyn.rba.acceleration(dof)))
		{
			std::cout << std::endl << "acce: "<< dyn.rba.acceleration(dof) << std::endl;
			std::cout << "vel: "<< dyn.rbs.velocity[dof] << std::endl;
			std::cout << "for: "<< dyn.joints.elements[dof].effort << std::endl;
		}

		for( int i=0; i < queueOfDyn.size(); i++)
			{
				dyn = queueOfDyn.front();
				queueOfDyn.pop();

				states(0,i) = dyn.joints.elements[dof].effort;
				states(3,i) = 1;

				if(dof<3)
					{

						acceleration(0,i) = dyn.rba.acceleration(dof);
						states(1,i) = dyn.rbs.velocity[dof] * fabs(double(dyn.rbs.velocity[dof]));
						states(2,i) = dyn.rbs.velocity[dof];
						if((isnan)((double)dyn.rba.acceleration(dof)))
						{
							std::cout << "f("<< i << "): " << states(0,i) << " " << states(1,i)<< " " << states(2,i) << " " << states(3,i) << std::endl;
							std::cout << "acceleration("<< i << "): " << acceleration(0,i) << std::endl;
						}

					}
				else
					{
						acceleration(0,i) = dyn.ang_rba.acceleration(dof);
						states(1,i) = dyn.rbs.angular_velocity[dof] * fabs(double(dyn.rbs.angular_velocity[dof]));
						states(2,i) = dyn.rbs.angular_velocity[dof];
					}

				queueOfDyn.push(dyn);
			}
		std::cout << "states size " << states.rows() << "x" << states.cols() << std::endl;
		std::cout << "acceleration size " << acceleration.rows() << "x" << acceleration.cols() << std::endl;
	}


	void LS_Pseudo::convert_data2 (std::queue<adap_samples_input::DynamicAUV> &queueOfDyn, base::MatrixXd &states, base::MatrixXd &forces)
	{
		states.resize(4,queueOfDyn.size());
		forces.resize(1,queueOfDyn.size());
		adap_samples_input::DynamicAUV dyn;


		for( int i=0; i < queueOfDyn.size(); i++)
		{
			dyn = queueOfDyn.front();
			queueOfDyn.pop();

			forces(0,i) = dyn.joints.elements[dof].effort;
			states(3,i) = 1;

			if(dof<3)
			{
				states(0,i) = dyn.rba.acceleration(dof);
				states(1,i) = dyn.rbs.velocity[dof] * fabs(double(dyn.rbs.velocity[dof]));
				states(2,i) = dyn.rbs.velocity[dof];
				//std::cout << "f("<< i << "): " << states(0,i) << " " << states(1,i)<< " " << states(2,i) << " " << states(3,i) << std::endl;
			}
			else
			{
				states(0,i) = dyn.ang_rba.acceleration(dof);
				states(1,i) = dyn.rbs.angular_velocity[dof] * fabs(double(dyn.rbs.angular_velocity[dof]));
				states(2,i) = dyn.rbs.angular_velocity[dof];
			}

			queueOfDyn.push(dyn);
		}
	}

	void LS_Pseudo::pseudo(base::MatrixXd &output, base::MatrixXd &states, base::MatrixXd &parameters)
	{	// acceleration = parameters.transpose()*states; [acceleration]1xn; [parameters]4x1; [states]4xn
		// parameters = (acceleration * states.pseudoinverse()).transpose()
		// states.pseudoinverse() = states.transpose() * (states*states.transpose()).inverse()

		if(output.cols() != states.cols())
			std::cout << std::endl << "Matrices of acceleration and states should have the same size "<< std::endl << std::endl;
		else
		{
			parameters.resize(4,1);

			base::MatrixXd gram = states * states.transpose();

			for(int j=0;j<output.size(); j++)
			{
				if((isnan)((double)output(0,j)))
				{
					std::cout << std::endl << "nan in pos: "<< j << std::endl;
					std::cout << "output("<< j << "): "<< output(0,j) << std::endl;
				}
			}

			base::MatrixXd para_trans = output * states.transpose() * gram.inverse();

			parameters = para_trans.transpose();
		}
	}


	void LS_Pseudo::svd(base::MatrixXd &output, base::MatrixXd &states, base::MatrixXd &parameters)
		{
			// [states]4xn; [acceleration]1xn; [parameters]4x1;
			// A*x = b;  [x]4x1; [A]nx4; [b]nx1
			// A = U S V^* ; U[nxn]; S[4x1]=diag{lamdai}/i=1...4; V[4x4]
			parameters.resize(4,1);
			base::MatrixXd A = states.transpose();
			base::MatrixXd b = output.transpose();
			base::MatrixXd gram = A.transpose()*A;

			Eigen::JacobiSVD<Eigen::MatrixXd> svd(A, Eigen::ComputeThinU | Eigen::ComputeThinV);

			parameters = svd.solve(b);
			std::cout<< std::endl << " S "<< std::endl << svd.singularValues() << std::endl;
			std::cout<< std::endl << "svd rank A " << svd.rank() << std::endl;
		//	svd.setThreshold(1e-6);
		//	std::cout<< std::endl << "rank A w/ Threshold(1e-6)" << svd.rank() << std::endl;
		//	svd.setThreshold(1e-5);
		//	std::cout<< std::endl << "rank A w/ Threshold(1e-5)" << svd.rank() << std::endl;
		//	svd.setThreshold(1e-4);
		//	std::cout<< std::endl << "rank A w/ Threshold(1e-4)" << svd.rank() << std::endl;
		//	svd.setThreshold(1e-3);
		//	std::cout<< std::endl << "rank A w/ Threshold(1e-3)" << svd.rank() << std::endl;
		//	parameters = svd.solve(b);
		//	svd.setThreshold(1e-2);
		//	std::cout<< std::endl << "rank A w/ Threshold(1e-2)" << svd.rank() << std::endl;
		//	parameters = svd.solve(b);
		//	std::cout<< std::endl << " S "<< std::endl << svd.singularValues() << std::endl;
		//	svd.setThreshold(1e-1);
		//	std::cout<< std::endl << "rank A w/ Threshold(1e-1)" << svd.rank() << std::endl;
	//		svd.setThreshold(1e0);
	//		std::cout<< std::endl << "rank A w/ Threshold(1e0)" << svd.rank() << std::endl;
			//std::cout<< std::endl << "A*A^T "<< std::endl <<  gram << std::endl;

		//	Eigen::EigenSolver<Eigen::MatrixXd> eigensolver(gram);


			//std::cout<< std::endl << "rank A " << svd.matrixU() << std::endl;
			//std::cout<< std::endl << "V " << svd.matrixV() << std::endl;

		}

	void LS_Pseudo::pcr(base::MatrixXd &output, base::MatrixXd &states, base::MatrixXd &parameters)
	{
		// [states]4xn; [acceleration]1xn; [parameters]4x1;
		// A*x = b;  [x]4x1; [A]nx4; [b]nx1
		// A = U S V^* ; U[nxn]; S[4x1]=diag{lamdai}/i=1...4; V[4x4]
		parameters.resize(4,1);
		base::MatrixXd A = states.transpose();
		base::MatrixXd b = output.transpose();
		base::MatrixXd gram = A.transpose()*A;
		int rows = A.rows();
		base::MatrixXd Ac;
		base::MatrixXd bc;
		Ac.resize(rows,4); bc.resize(rows,1);

		// centering: Cn = In - 1/n*11^T. ex: C2 = [1,0;0,1] - 1/2*[1,1;1,1]
		base::MatrixXd II;
		II.resize(rows,1);
		for(int i=0;i<rows;i++)
			II(i,0)=1;
		base::MatrixXd Cn = Eigen::MatrixXd::Identity(rows,rows) - 1/rows*II*II.transpose();

		Ac.col(0) = Cn*A.col(0);
		Ac.col(1) = Cn*A.col(1);
		Ac.col(2) = Cn*A.col(2);
		Ac.col(3) = Cn*A.col(3);
		bc		  = Cn*b;

		Eigen::JacobiSVD<Eigen::MatrixXd> svd(Ac, Eigen::ComputeThinU | Eigen::ComputeThinV);
		//parameters = svd.solve(bc);
		base::VectorXd s = svd.singularValues();
		base::MatrixXd S = s.asDiagonal();
		base::MatrixXd V = svd.matrixV();
		base::MatrixXd U = svd.matrixU();
		std::cout<< std::endl << " S "<< std::endl << S << std::endl;
		std::cout<< std::endl << " V "<< std::endl << V << std::endl;


		base::MatrixXd Wk;
		base::MatrixXd Vk;
		Wk.resize(rows,3);
		Vk.resize(4,3);

		Vk.col(0) = V.col(0);
		Vk.col(1) = V.col(1);
		Vk.col(2) = V.col(2);

		Wk = Ac*Vk;
		base::MatrixXd gramW = Wk.transpose()*Wk;
		base::MatrixXd pcr = gramW.inverse()*Wk.transpose()*bc;

		parameters = Vk*pcr;
	}

	void LS_Pseudo::convert_param(base::MatrixXd &parameters, base::MatrixXd &real_parameters)
	{
		real_parameters.resize(4,1);
		if (parameters(0,0)!=0)
		{
			std::cout<<std::endl<< "convert parameters " <<std::endl <<std::endl;
			real_parameters(0,0) = 1/parameters(0,0);
			real_parameters(1,0) = -parameters(1,0)/parameters(0,0);
			real_parameters(2,0) = -parameters(2,0)/parameters(0,0);
			real_parameters(3,0) = -parameters(3,0)/parameters(0,0);
		}
		else
			std::cout << std::endl << "Could not transform into real parameters"<< std::endl << std::endl;

	}

	void LS_Pseudo::ls_solution (std::queue<adap_samples_input::DynamicAUV> &queueOfDyn, base::MatrixXd &real_parameters, double &error)
	{
		static base::MatrixXd parameters;
		static base::MatrixXd states;
		static base::MatrixXd acceleration;
		static base::MatrixXd forces;

		parameters.resize(4,1);
		states.resize(4,queueOfDyn.size());
		acceleration.resize(1,queueOfDyn.size());
		forces.resize(1,queueOfDyn.size());


		convert_data(queueOfDyn, states, acceleration);
		pseudo(acceleration, states, parameters);

		convert_param(parameters, real_parameters);
		relativ_error(acceleration, states, parameters, error);

		std::cout<<std::endl<< "lamped parameters by pseudo_inverse: " <<std::endl << parameters <<std::endl;
		std::cout<<std::endl<< "parameters by pseudo_inverse: " <<std::endl << real_parameters <<std::endl;
		//std::cout<<std::endl<< "error pseudo: "<< error <<std::endl<< std::endl;

	//	convert_data2(queueOfDyn, states, forces);
	//	pseudo(forces, states, parameters);

	//	std::cout<<std::endl<< "parameters by pseudo_inverse2: " <<std::endl << parameters <<std::endl;


		svd(acceleration, states, parameters);
		convert_param(parameters, real_parameters);
	//	relativ_error(acceleration, states, parameters, error);

		std::cout<<std::endl<< "parameters by Singular Value Decomposition (SVD) : " <<std::endl << real_parameters <<std::endl;
		//std::cout<<std::endl<< "error svd: "<< error <<std::endl<< std::endl;

		pcr(acceleration, states, parameters);
	//	relativ_error(acceleration, states, parameters, error);
		convert_param(parameters, real_parameters);
		std::cout<<std::endl<< "parameters by Principal Component Regression (PCR): " <<std::endl << real_parameters <<std::endl;
		//std::cout<<std::endl<< "error pcr: "<< error <<std::endl<< std::endl;


		convert_data2(queueOfDyn, states, forces);
		pseudo(forces, states, parameters);

		std::cout<<std::endl<< "parameters by pseudo_inverse2: " <<std::endl << parameters <<std::endl;

		svd(forces, states, parameters);
		std::cout<<std::endl<< "parameters by SVD_2 : " <<std::endl << real_parameters <<std::endl;

		pcr(forces, states, parameters);
		std::cout<<std::endl<< "parameters by PCR_2: " <<std::endl << real_parameters <<std::endl;


	}

	void LS_Pseudo::relativ_error(base::MatrixXd &output, base::MatrixXd &states, base::MatrixXd &parameters, double &error)
	{	// [states]4xn; [acceleration]1xn; [parameters]4x1;
		// A*x = b; [x]4x1; [A]nx4; [b]nx1
		// relative_error = (A*x - b).norm() / b.norm()
		base::MatrixXd A = states.transpose();
		base::MatrixXd b = output.transpose();
		error = (A*parameters - b).norm() / b.norm();
	}



}
