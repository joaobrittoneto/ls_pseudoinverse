#ifndef _LS_PSEUDOINVERSE_HPP_
#define _LS_PSEUDOINVERSE_HPP_

#include <iostream>
#include <vector>
#include <queue>
#include <math.h>
#include "base/samples/RigidBodyState.hpp"
#include "base/samples/RigidBodyAcceleration.hpp"
#include "base/samples/Joints.hpp"
#include "base/Eigen.hpp"
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/SVD>
#include "adap_parameters_estimator/adap_dataTypes.h"
#include "adap_samples_input/samples_dataType.h"

namespace ls_pseudoinverse
{
	class LS_Pseudo
	{
		public: 
			
			LS_Pseudo(adap_parameters_estimator::DOFS _dof);
			~LS_Pseudo();
			//void convert_states (std::queue<base::samples::Joints> &queueOfForces, std::queue<base::samples::RigidBodyState> &queueOfRBS, base::MatrixXd &states);
			//void convert_acceleration (std::queue<base::samples::RigidBodyAcceleration> &queueOfRBA, base::MatrixXd &acceleration);
			void convert_data (std::queue<adap_samples_input::DynamicAUV> &queueOfDyn, base::MatrixXd &states, base::MatrixXd &acceleration);
			void convert_data2 (std::queue<adap_samples_input::DynamicAUV> &queueOfDyn, base::MatrixXd &states, base::MatrixXd &forces);
			void pseudo (base::MatrixXd &output, base::MatrixXd &states, base::MatrixXd &parameters);
			void svd (base::MatrixXd &output, base::MatrixXd &states, base::MatrixXd &parameters);
			void pcr (base::MatrixXd &output, base::MatrixXd &states, base::MatrixXd &parameters);
			void convert_param (base::MatrixXd &parameters, base::MatrixXd &real_parameters);
			void ls_solution (std::queue<adap_samples_input::DynamicAUV> &queueOfDyn, base::MatrixXd &real_parameters, double &error);
			void relativ_error(base::MatrixXd &output, base::MatrixXd &states, base::MatrixXd &parameters, double &error);

			
		private:
			adap_parameters_estimator::DOFS dof;
		
	};

} // end namespace ls_pseudoinverse

#endif // _LS_PSEUDOINVERSE_HPP_
