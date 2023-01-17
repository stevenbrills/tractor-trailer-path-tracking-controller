#include <ct/optcon/optcon.h>

/*
Class for the LQR based path tracking control and stabilisation for a tractor-trailer assembly robot
Initialize with test vehicle parameters:
1. Trcator Wheelbase
2. Trailer Wheelbase
3. Tractor Hitch Offset
4. Forward Lookahead Radius
5. Backward Lookahead Radius
6. Velocity
7. Q Matrix
8. R Matrix
*/
class TractorTrailerController{
    
    private:
        // Vehicle Parameters
        float m_tractor_wheelbase;
        float m_trailer_wheelbase;
        float m_tractor_hitch_offset;
        float m_forward_lookahead_radius;
        float m_backward_lookahead_radius;
        float m_velocity;

        // Simulator Parameters
        int m_lookahead_inflation_steps=50;

        // LQR Controller Parameters
        float m_q_matrix;
        float m_r_matrix;

        // Derived Parameters
        float m_max_alpha_e;

//////////////////////////////////////////////////////////////////////////////
/////////////////////// Utility Methods Declared Here ////////////////////////
//////////////////////////////////////////////////////////////////////////////

        /* Method that checks whether two decimal values are within epsilon of each other */
        bool checkDoubleEqual(
            double& a,
            double& b
        ){
            double epsilon = 0.001;
            if (fabs(a-b)<=epsilon){
                return true;
            }
            return false;
        }

        /* Method that computes the orientation of the tractor given the minimal state */
        double getTractorOrientation(
            const std::vector<double>& q_current
        ){
            if(q_current[3]>=0){
                return q_current[2] + (M_PI - q_current[3]);
            }
            return (q_current[2] - (M_PI - fabs(q_current[3])));    
        }

        /* Method that computes the location of the tractor's axle center given the minimal state */
        void getTractorAxleCenter(
            std::vector<double>& q
        ){
            double tractor_theta = getTractorOrientation(q);
            q[4] = q[0] + m_trailer_wheelbase*cos(q[2]) + m_tractor_hitch_offset*cos(tractor_theta);
            q[5] = q[1] + m_trailer_wheelbase*sin(q[2]) + m_tractor_hitch_offset*sin(tractor_theta);
        }

        /* Method that wraps any input angle to -PI to +PI */
        double wrapAngle(
            const double input_angle
        ){

            if((input_angle >= (-1*M_PI)) && (input_angle <= (M_PI))){
                return input_angle;
            }

            double remainder = fmod(input_angle,M_PI);

            if (remainder < 0){
                return (M_PI + remainder);
            }
            
            return ((-1*M_PI) + remainder);
        }

        /* Method that limits the commanded steering angle to the maximum allowed steering angle */
        double clipToAlphaE(
            double angle
        ){
            if(angle<(-1*m_max_alpha_e)){
                return (-1*m_max_alpha_e);
            }

            if(angle>(1*m_max_alpha_e)){
                return (1*m_max_alpha_e);
            }

            return angle;

        }


//////////////////////////////////////////////////////////////////////////////
//////////////////// Controller Methods Declared Here ////////////////////////
//////////////////////////////////////////////////////////////////////////////

        /* Method that computes the state matrix for LQR (Matrix A) */
        double getStateMatrix(
            const double& alpha_e,
            const double& beta_e
        ){
            double r1 = m_tractor_wheelbase/tan(alpha_e);
            double psi = -1*(alpha_e/fabs(alpha_e))*atan(m_tractor_hitch_offset/fabs(r1));
            double r2 = m_tractor_hitch_offset/sin(psi);
            double r3 = m_trailer_wheelbase/sin(psi-beta_e);
            double A;

            A = m_velocity*(1/m_tractor_wheelbase)*tan(alpha_e)*(((m_tractor_hitch_offset/m_trailer_wheelbase)*((-1*sin(psi)*sin(beta_e) -1*cos(psi)*cos(beta_e))/(sin(psi)))));


            return A;
        }

        /* Method that computes the input matrix for LQR (Matrix B) */
        double getInputMatrix(
            const double& alpha_e,
            const double& beta_e
        ){
            double r1 = m_tractor_wheelbase/tan(alpha_e);
            double psi = -1*(alpha_e/fabs(alpha_e))*atan(m_tractor_hitch_offset/fabs(r1));
            double r2 = m_tractor_hitch_offset/sin(psi);
            double r3 = m_trailer_wheelbase/sin(psi-beta_e);
            double B;

            double a = (1/m_tractor_wheelbase)*(1/pow(cos(alpha_e),2));
            double b = (m_tractor_hitch_offset/m_trailer_wheelbase)*sin(psi-beta_e)/sin(psi);
            double c = tan(alpha_e)/m_tractor_wheelbase;
            double d = (m_tractor_hitch_offset/m_trailer_wheelbase)*(sin(beta_e)*m_tractor_hitch_offset*m_tractor_wheelbase);
            double e = pow(sin(alpha_e),2)*pow(sin(psi),2)*(pow(r1,2) + pow(m_tractor_hitch_offset,2));

            B = -1*m_velocity*(c*(d/e)  + a*(b-1));

            return B;
        }

        /* Method that computes the rate of change of the tractor-trailer cart angle given the steering and tractor-trailer
        cart angles */
        double getBetaDot(
            double& alpha,
            double& beta
        ){
            double r1 = m_tractor_wheelbase/tan(alpha);
            double psi = -1*(alpha/fabs(alpha))*atan(m_tractor_hitch_offset/fabs(r1));
            double r2 = m_tractor_hitch_offset/sin(psi);
            double r3 = m_trailer_wheelbase/sin(beta-psi);
            double beta_dot = m_velocity*((1/r1) + ((r2)/(r1*r3)));
            double omega_2 = ((r2)/(r1*-1*r3));

            return -1*beta_dot;
        }

        /* Method that computes the gain for the equillibrium state given the equillibrium steering and tractor-trailer
        cart angles */
        double getGain(
            const double& beta_e,
            const double& alpha_e
        ){
            const size_t state_dim = 1;
            const size_t control_dim = 1;

            ct::core::StateMatrix<state_dim> A;
            ct::core::ControlMatrix<control_dim> B;
            ct::core::StateMatrix<state_dim> Q;
            ct::core::ControlMatrix<control_dim> R;

            A(0,0) = getStateMatrix(alpha_e, beta_e);
            B(0,0) = getInputMatrix(alpha_e, beta_e);
            Q(0,0) = m_q_matrix;
            R(0,0) = m_r_matrix;

            ct::optcon::LQR<state_dim, control_dim> lqrSolver;
            ct::core::FeedbackMatrix<state_dim, control_dim> K;

            lqrSolver.compute(Q, R, A, B, K);

            return K(0,0);
        }

//////////////////////////////////////////////////////////////////////////////
///////////////////// Simulator Methods Defined Here ////////////////////////
//////////////////////////////////////////////////////////////////////////////

        /* Method that returns the tractor-trailer cart angle at equillibruim for the given equillibrium steering angle alpha_e */
        double getBetaEGivenAlpha(
            const double& alpha_e
        ){
            double r1 = m_tractor_wheelbase/tan(alpha_e);
            double r2 = sqrt(pow(r1,2) + pow(m_tractor_hitch_offset,2));
            double psi = atan(m_tractor_hitch_offset/fabs(r1));
            double theta_1 = asin(fabs(r1)/r2);
            double theta_2 = acos(m_trailer_wheelbase/r2);
            double beta_e = (alpha_e/fabs(alpha_e))*(theta_1 + theta_2);
            return beta_e;

        }

        /* Method that takes in the current state, the piecewise linear path reference, the segment start and end index of the
        particular segment on the picewise linear path, the direction flag and a return parameter that stores the computed
        intersection point if any */
        bool findIntersectionPoint(
            const std::vector<double>& q_current, 
            const std::vector<std::vector<double>>& piecewise_linear, 
            const int& id_start, 
            const int& id_goal,
            const bool& is_forward,
            std::vector<double>& intersection_point
        ){

            // if the motion is reverse, use the trailer centerd lookahead circle
            if (!is_forward){

                // Set up co-efficients for quadratic equation
                double a = pow((piecewise_linear[id_goal][0] - piecewise_linear[id_start][0]),2) + 
                pow((piecewise_linear[id_goal][1] - piecewise_linear[id_start][1]),2);

                double b = 2*(
                    ((piecewise_linear[id_start][0] - q_current[0])*(piecewise_linear[id_goal][0] - piecewise_linear[id_start][0])) + 
                    ((piecewise_linear[id_start][1] - q_current[1])*(piecewise_linear[id_goal][1] - piecewise_linear[id_start][1]))
                );

                double c = pow((piecewise_linear[id_start][0] - q_current[0]),2) + pow((piecewise_linear[id_start][1] - q_current[1]),2) - pow(m_backward_lookahead_radius,2);

                double det = pow(b,2) - (4*a*c);


                if(det>=0){

                    double t1 = ((-1*b) + sqrt(det))/(2*a);
                    double t2 = ((-1*b) - sqrt(det))/(2*a);

                    if (checkDoubleEqual(t1, t2)){


                        intersection_point[0] = (piecewise_linear[id_goal][0]*t1) + (1-t1)*piecewise_linear[id_start][0];
                        intersection_point[1] = (piecewise_linear[id_goal][1]*t1) + (1-t1)*piecewise_linear[id_start][1];
                        if((t1>=0) && (t1 <=1)){
                            return true; 
                        } //intersection detected
                        else{
                            return false;
                        }

                    }
                    else{

                        // t1 will always be greater than t2 and since the second point of the piecewise linear
                        // line will be the direction of travel, using t1 should give the point to drive to

                        double x_intersection1 = (piecewise_linear[id_goal][0]*t1) + (1-t1)*piecewise_linear[id_start][0];
                        double y_intersection1 = (piecewise_linear[id_goal][1]*t1) + (1-t1)*piecewise_linear[id_start][1];

                        double x_intersection2 = (piecewise_linear[id_goal][0]*t2) + (1-t2)*piecewise_linear[id_start][0];
                        double y_intersection2 = (piecewise_linear[id_goal][1]*t2) + (1-t2)*piecewise_linear[id_start][1];

                        intersection_point[0] = (piecewise_linear[id_goal][0]*t1) + (1-t1)*piecewise_linear[id_start][0];
                        intersection_point[1] = (piecewise_linear[id_goal][1]*t1) + (1-t1)*piecewise_linear[id_start][1];

                        if((t1>=0) && (t1 <=1)){
                            return true; 
                        } //intersection detected
                        else{
                            return false;
                        }

                    }

                }
                return false;
            }

            // If the motion is forward, use tractor lookahead circle
            else{

                // Set up co-efficients for quadratic equation
                double a = pow((piecewise_linear[id_goal][0] - piecewise_linear[id_start][0]),2) + 
                pow((piecewise_linear[id_goal][1] - piecewise_linear[id_start][1]),2);

                double b = 2*(
                    ((piecewise_linear[id_start][0] - q_current[4])*(piecewise_linear[id_goal][0] - piecewise_linear[id_start][0])) + 
                    ((piecewise_linear[id_start][1] - q_current[5])*(piecewise_linear[id_goal][1] - piecewise_linear[id_start][1]))
                );

                double c = pow((piecewise_linear[id_start][0] - q_current[4]),2) + pow((piecewise_linear[id_start][1] - q_current[5]),2) - pow(m_forward_lookahead_radius,2);

                double det = pow(b,2) - 4*a*c;

                if(det>=0){

                    double t1 = ((-1*b) + sqrt(det))/(2*a);
                    double t2 = ((-1*b) - sqrt(det))/(2*a);

                    if (checkDoubleEqual(t1, t2)){

                        if((t1>=0) && (t1 <=1)){
                        intersection_point[0] = (piecewise_linear[id_goal][0]*t1) + (1-t1)*piecewise_linear[id_start][0];
                        intersection_point[1] = (piecewise_linear[id_goal][1]*t1) + (1-t1)*piecewise_linear[id_start][1];
                        return true; // Intersection detected
                        }
                        else{
                            return false;
                        }



                    }
                    else{

                        // t1 will always be greater than t2 and since the second point of the piecewise linear
                        // line will be the direction of travel, using t1 should give the point to drive to

                        double x_intersection1 = (piecewise_linear[id_goal][0]*t1) + (1-t1)*piecewise_linear[id_start][0];
                        double y_intersection1 = (piecewise_linear[id_goal][1]*t1) + (1-t1)*piecewise_linear[id_start][1];

                        double x_intersection2 = (piecewise_linear[id_goal][0]*t2) + (1-t2)*piecewise_linear[id_start][0];
                        double y_intersection2 = (piecewise_linear[id_goal][1]*t2) + (1-t2)*piecewise_linear[id_start][1];
                        if((t1>=0) && (t1 <=1)){
                            intersection_point[0] = (piecewise_linear[id_goal][0]*t1) + (1-t1)*piecewise_linear[id_start][0];
                            intersection_point[1] = (piecewise_linear[id_goal][1]*t1) + (1-t1)*piecewise_linear[id_start][1];
                            return true;
                        }
                        else{
                            return false;
                        }
                    }
                }

                return false;

            }
            
            return false;
        }

        /* Method that takes in the current state, intersection point and the direction flag and returns
        the desired tractor-trailer cart angle */
        double getBetaDesired(
            const std::vector<double>& q_current, 
            const std::vector<double>& intersection_point,
            const bool& is_forward 
        ){

            // Transform the intersection point into the local frame of the vehicle
            double tx = -1*((q_current[0]*cos(q_current[2])) + (q_current[1]*sin(q_current[2])));
            double ty = -1*((-1*q_current[0]*sin(q_current[2])) + (q_current[1]*cos(q_current[2])));

            double projected_distance_on_axle_normal = (intersection_point[0]*cos(q_current[2])) + (intersection_point[1]*sin(q_current[2])) + tx;

            double projected_distance_on_axle = (-1*intersection_point[0]*sin(q_current[2])) + (intersection_point[1]*cos(q_current[2])) + ty;

            double theta_e = atan2(projected_distance_on_axle, projected_distance_on_axle_normal);

            double turning_circle_radius = fabs(m_backward_lookahead_radius/(2*sin(theta_e)));

            double beta_1 = atan(turning_circle_radius/m_trailer_wheelbase);

            double beta_2 = acos(m_tractor_hitch_offset/sqrt(pow(m_trailer_wheelbase,2)+pow(turning_circle_radius,2)));

            double beta_d = (theta_e/fabs(theta_e))*(beta_1+beta_2);

            return beta_d;
        }

        /* Method that returns the steering angle when in forward motion */
        double getAlphaForForwardMotion(
            const std::vector<double>& intersection_point,
            const std::vector<double>& q_current
        ){

            // Frame to be transformed into is the tractor axle frame
            // Use the axle center x,y and theta of tractor
            double theta_tractor = getTractorOrientation(q_current);

            // Transform the intersection point into the local frame of the vehicle
            double tx = -1*((q_current[4]*cos(theta_tractor)) + (q_current[5]*sin(theta_tractor)));
            double ty = -1*((-1*q_current[4]*sin(theta_tractor)) + (q_current[5]*cos(theta_tractor)));

            // Projected distances 
            double projected_distance_on_axle_normal = (intersection_point[0]*cos(theta_tractor)) + (intersection_point[1]*sin(theta_tractor)) + tx;
            double projected_distance_on_axle = (-1*intersection_point[0]*sin(theta_tractor)) + (intersection_point[1]*cos(theta_tractor)) + ty;

            double r1 = (pow(m_forward_lookahead_radius,2))/(2*fabs(projected_distance_on_axle));

            double alpha=0;

            if(projected_distance_on_axle==0){
                alpha = 0;
            }
            else{
                alpha = atan(m_tractor_wheelbase/r1)*(projected_distance_on_axle/fabs(projected_distance_on_axle));
            }

            alpha = wrapAngle(alpha);

            return alpha;
        }

        /* Method that returns the equillibrium value of the steering angle given the equillibrium tractor-trailer cart angle*/
        double getAlphaE(
            double& beta_e
        ){

            double r3 = (m_tractor_hitch_offset/cos(fabs(beta_e) - M_PI_2))*(((m_trailer_wheelbase*cos(M_PI - fabs(beta_e)))/m_tractor_hitch_offset) + 1);

            double r2 = sqrt(pow(m_trailer_wheelbase,2) + pow(r3,2));

            double r1 = sqrt(pow(r2,2) - pow(m_tractor_hitch_offset,2));

            double alpha_e = atan(m_tractor_wheelbase/r1)*(beta_e/fabs(beta_e));

            return alpha_e;
        }

        /* Method that returns the rate of change of state by taking in the current state and the steering angle */
        std::vector<double> qDot(
            const std::vector<double>& q_current,
            const double& alpha
        ){

            std::vector<double> q_dot(4,0); // Variable to store the rate of change of state

            double va=0;
            double r1=0;
            double psi=0;
            double r2=0;
            double r3=0;
            double vb=0;
            double omega_2=0;
            double trailer_theta_dot = 0;
            double hitch_turn_radius_zero_alpha = 0;
            double phi = 0;
            double beta_dot = 0;
            double omega_1 = 0;
            
            r1 = m_tractor_wheelbase/tan(alpha);
            psi = -1*(alpha/fabs(alpha))*atan(m_tractor_hitch_offset/fabs(r1));
            va = (m_velocity/fabs(m_velocity))*fabs(((m_velocity*m_tractor_hitch_offset*tan(alpha))/(m_tractor_wheelbase*sin(psi))));
            vb = va*fabs(cos(psi-q_current[3]));
            r2 = m_tractor_hitch_offset/sin(psi);
            r3 = m_trailer_wheelbase/sin(psi-q_current[3]);
            omega_1 = m_velocity*(1/r1);
            omega_2 = m_velocity*((r2)/(r1*r3));
            trailer_theta_dot = omega_2;
            beta_dot = m_velocity*(((r2)/(r1*r3)) - (1/r1));

            // Special case when alpha equals zero
            if(alpha==0){
                va = m_velocity;
                psi = 0;
                phi = M_PI_2 - (fabs(q_current[3]) - M_PI_2);
                hitch_turn_radius_zero_alpha = m_trailer_wheelbase/sin(phi);
                trailer_theta_dot = va/hitch_turn_radius_zero_alpha;
                beta_dot = trailer_theta_dot;
                vb = va*fabs(cos(psi-q_current[3]));

                if(sin(phi)==0){
                    beta_dot = 0;
                    trailer_theta_dot = 0;
                }
            }

            q_dot[0] = vb*cos(q_current[2]);
            q_dot[1] = vb*sin(q_current[2]);
            q_dot[2] = trailer_theta_dot;
            q_dot[3] = beta_dot;

            return q_dot;
        }

        /* Method that integrates the state by taking in the steering angle (alpha), the current state and the timestep size */
        std::vector<double> rk4Integrator(
        const double& alpha,
        const std::vector<double> q_current,
        const double& timestep
        ){

            std::vector<double> k1(4,0);
            std::vector<double> k2(4,0);
            std::vector<double> k3(4,0);
            std::vector<double> k4(4,0);


            k1 = qDot(q_current, alpha);
            std::vector<double> q_delta1(4,0);
            q_delta1[0] = q_current[0] + k1[0]*(timestep/2);
            q_delta1[1] = q_current[1] + k1[1]*(timestep/2);
            q_delta1[2] = q_current[2] + k1[2]*(timestep/2);
            q_delta1[3] = q_current[3] + k1[3]*(timestep/2);

            k2 = qDot(q_delta1, alpha);
            std::vector<double> q_delta2(4,0);
            q_delta2[0] = q_current[0] + k2[0]*(timestep/2);
            q_delta2[1] = q_current[1] + k2[1]*(timestep/2);
            q_delta2[2] = q_current[2] + k2[2]*(timestep/2);
            q_delta2[3] = q_current[3] + k2[3]*(timestep/2);

            k3 = qDot(q_delta2, alpha);
            std::vector<double> q_delta3(4,0);
            q_delta3[0] = q_current[0] + k3[0]*(timestep);
            q_delta3[1] = q_current[1] + k3[1]*(timestep);
            q_delta3[2] = q_current[2] + k3[2]*(timestep);
            q_delta3[3] = q_current[3] + k3[3]*(timestep);


            k4 = qDot(q_delta3, alpha);

            std::vector<double> q_next(6,0);

            q_next[0] = q_current[0] + ((k1[0]/6)+(k2[0]/3)+(k3[0]/3)+(k4[0]/6))*timestep;
            q_next[1] = q_current[1] + ((k1[1]/6)+(k2[1]/3)+(k3[1]/3)+(k4[1]/6))*timestep;
            q_next[2] = q_current[2] + ((k1[2]/6)+(k2[2]/3)+(k3[2]/3)+(k4[2]/6))*timestep;
            q_next[3] = q_current[3] + ((k1[3]/6)+(k2[3]/3)+(k3[3]/3)+(k4[3]/6))*timestep;

            q_next[2] = wrapAngle(q_next[2]);
            q_next[3] = wrapAngle(q_next[3]);

            // For every calculated future state, also populate the axle center of the tractor
            getTractorAxleCenter(q_next);

            return q_next;

        }

        /* Method that takes in the current state, the direction of travel, a return parameter that stores the intersection point
        and the entire piecewise linear path. It iterates through all the segments of the path and returns the correct intersection
        point based on the direction of travel.*/
        bool getIntersectionPointAlongPiecewiseLinear(
            const std::vector<double>& q_current, 
            const std::vector<std::vector<double>>& piecewise_linear, 
            const bool& is_forward,
            std::vector<double>& intersection_point
        ){

            bool intersection_detected_in_prev_segment = false;
            intersection_point[0] = q_current[0];
            intersection_point[1] = q_current[1];

            std::vector<double> previous_intersection_point(2,0);


            for (int i=0; i<(piecewise_linear.size()-1); i++){

                if(intersection_detected_in_prev_segment && findIntersectionPoint(q_current, piecewise_linear, i, i+1, is_forward, intersection_point)){
                    return true;
                };

                if(intersection_detected_in_prev_segment && (!findIntersectionPoint(q_current, piecewise_linear, i, i+1, is_forward, intersection_point))){
                    intersection_point = previous_intersection_point;
                    return true;
                };

                if(findIntersectionPoint(q_current, piecewise_linear, i, i+1, is_forward, intersection_point)){
                    intersection_detected_in_prev_segment = true;
                    previous_intersection_point = intersection_point;
                };

                if(findIntersectionPoint(q_current, piecewise_linear, i, i+1, is_forward, intersection_point) && (i==(piecewise_linear.size()-2))){
                    return true;
                };
            }

            return false;
        }

        /* Method that runs the integration over the path for one segment of the entire path. This method is repatedly called 
        by the public "forwardSimulator" method until the entire piecewise length path has been simulated over or the intersection
        point goes outside the lookahead circle. */
        std::vector<std::vector<double>> segmentSimulator(
            std::vector<double> q_init,
            std::vector<std::vector<double>> segment,
            bool is_forward
        ){

            // Set m_velocity sign based on direction
            if(is_forward){
                m_velocity = fabs(m_velocity);
            }
            else{
                m_velocity = -1*fabs(m_velocity);
            }

            // Initialize variables for simulation
            double beta_desired;        // Desired tractor-trailer cart angle
            double beta_e;              // Equillibrium tractor-trailer cart angle
            double alpha_e;             // Equillibrium steering angle
            double alpha;               // Steering angle
            std::vector<double> q_next; //Next state


            // Create the data structure to store the trajectory
            std::vector<std::vector<double>> trajectory;

            // Simulation parameters
            float beta_prop_gain = 0;   // gain on tractor trailer cart angle error
            double timestep = 0.001;    // integration time step

            // Check if input path is valid
            if ((segment.size()<2) || (segment.size()>2)){
                throw std::runtime_error("Size of segment in segmentwise simulator is less than or greater than 2 which is invalid. Size should exactly equal two! Segment must be defined by a pair of points.");
            }

            // Initialize the intersection point to be at the center of the rear axle
            std::vector<double> intersection_point{q_init[0], q_init[1]};

            // Set the current state to the initial state of the trailer
            std::vector<double> q_current=q_init;

            // Run the simulation till the intersection point of the look-ahead circle and piecewise linear path
            // reaches the last point along the piecewise linear path
            bool found_intersection_flag=true;  // flag for improper termination of simulation
            int while_loop_counter = 0;         // simulation step counter

            while(
                !(checkDoubleEqual(intersection_point[0], segment[segment.size()-1][0]) &&
            checkDoubleEqual(intersection_point[1], segment[segment.size()-1][1]))
            ){

                // Find the intersection points of the lookahead circle and piecewise linear path
                found_intersection_flag = getIntersectionPointAlongPiecewiseLinear(q_current, segment, is_forward, intersection_point);

                // If the new segment starts and the intersection cannot be found (precision issues for double and a result of 
                // termination critera of last segment), inflate the lookaheads by 1%
                if(while_loop_counter<=m_lookahead_inflation_steps && (!found_intersection_flag)){
                    // std::cout << "Inflation required!" << std::endl;
                    if(is_forward){
                        m_forward_lookahead_radius = m_forward_lookahead_radius*1.01;
                        found_intersection_flag = getIntersectionPointAlongPiecewiseLinear(q_current, segment, is_forward, intersection_point);
                        m_forward_lookahead_radius = m_forward_lookahead_radius/1.01;
                    }
                    else{
                        // For reverse motion, inflate by 5%
                        m_backward_lookahead_radius = m_backward_lookahead_radius*1.05;
                        // std::cout << "Inflated backward lookahead radius: " << m_backward_lookahead_radius << std::endl;
                        // std::cout << "Current trailer X: " << q_current[0] << "Current Trailer Y: " << q_current[1] << std::endl;
                        found_intersection_flag = getIntersectionPointAlongPiecewiseLinear(q_current, segment, is_forward, intersection_point);
                        // std::cout << "Dist to segment start: " << sqrt(pow(segment[0][0] - q_current[0],2)+pow(segment[0][1] - q_current[1],2)) << std::endl;
                        m_backward_lookahead_radius = m_backward_lookahead_radius/1.05;
                    }
                }

                if(!found_intersection_flag){
                    break;
                }

                // Here, use the intersection points to determine steering angle control inputs
                if(!is_forward){

                    // Compute the desired tractor-trailer cart angle using intersection
                    beta_desired = getBetaDesired(q_current, intersection_point, is_forward);

                    // Compute the tractor-trailer cart angle about which to linearize by multiplying gain on deviation
                    beta_e = beta_desired + beta_prop_gain*(beta_desired - q_current[3]);

                    // Get the value of equillibrium steering angle using the desired steering angle
                    alpha_e = getAlphaE(beta_e);

                    // Use alpha_e to compute steering angle input by multiplying by stabilising gain from LQR
                    alpha = clipToAlphaE(wrapAngle(alpha_e - getGain(beta_e, alpha_e)*(q_current[3] - beta_e)));

                    // Integrate motion through 
                    q_next = rk4Integrator(alpha, q_current, timestep);

                }
                else{

                    // Directly compute steering angle using pure-pursuit algorithm
                    alpha = clipToAlphaE(getAlphaForForwardMotion(intersection_point, q_current));
                    alpha = 1*alpha; // Multiply by gain if aggressive tracking is desired

                    // Catch invalid steering angle
                    if (std::isnan(alpha)){
                        throw std::runtime_error("Computed steering angle (alpha) was Nan");
                    }

                    // Integrate motion through 
                    q_next = rk4Integrator(alpha, q_current, timestep);

                }

                q_current = q_next;
                q_next.push_back(alpha);        // Add steering angle for debugging purposes
                trajectory.push_back(q_next);   // Append the next state into the trajectory
                while_loop_counter++;
            }

            std::cout << "Simulation complete!" << std::endl;
            if(!found_intersection_flag){
                std::cout << "Simulation cut off since no intersection point found" << std::endl;
            }

            return trajectory;
        }


    public:

        /* Getter method for forward lookahead radius parameter */
        double getForwardLookaheadRadius(){
            return this->m_forward_lookahead_radius;
        }

        /* Getter method for backward lookahead radius parameter */
        double getBackwardLookaheadRadius(){
            return this->m_backward_lookahead_radius;
        }

        /* Getter method for velocity parameter */
        double getVelocity(){
            return this->m_velocity;
        }

        /* Getter method for trailer wheelbase parameter */
        double getTrailerWheelbase(){
            return this->m_trailer_wheelbase;
        }

        /* Getter method for tractor wheelbase parameter */
        double getTractorWheelbase(){
            return this->m_tractor_wheelbase;
        }

        /* Getter method for tractor hitch offset parameter */
        double getTractorHitchOffset(){
            return this->m_tractor_hitch_offset;
        }

        /* Class Constructor: Takes in the tractor wheelbase, trailer wheelbase, tractor hitch offset, forward lookahead radius,
        backward lookahead radius, velocity and Q and R gains */
        TractorTrailerController(
        float tractor_wheelbase,
        float trailer_wheelbase,
        float tractor_hitch_offset,
        float forward_lookahead_radius, 
        float backward_lookahead_radius,
        float velocity,
        float Q,
        float R)
        {
            // Initialize simulator parameters
            this->m_tractor_wheelbase = tractor_wheelbase;
            this->m_trailer_wheelbase = trailer_wheelbase;
            this->m_tractor_hitch_offset = tractor_hitch_offset;
            this->m_forward_lookahead_radius = forward_lookahead_radius;
            this->m_backward_lookahead_radius = backward_lookahead_radius;
            this->m_velocity = velocity;

            // Initialize controller parameters
            this->m_q_matrix = Q;
            this->m_r_matrix = R;

            // Compute derived parameters
            this->m_max_alpha_e = atan(m_tractor_wheelbase/sqrt(pow(m_trailer_wheelbase,2) - pow(m_tractor_hitch_offset,2)));;
        }

        /* Method that iterates through each segment of the piecewise linear path and calls the segment simulator method on 
        each individual segment and returns the resulting trajectory from the simulation. */
        std::vector<std::vector<double>> forwardSimulator(
            std::vector<double> q_init,
            std::vector<std::vector<double>> piecewise_linear
        ){

            // Create the vector of vectors which stores the trajectory
            std::vector<std::vector<double>> final_trajectory;
            std::vector<std::vector<double>> temp_trajectory;

            // Add the initial configuration of the tractor-trailer into the trajectory
            final_trajectory.push_back(q_init);

            // Check if input piecewise linear path is valid
            if (piecewise_linear.size()<2){
                throw std::runtime_error("Piecewise linear input path is of insufficient size of elements!");
            }

            // Set the current state to the initial state of the trailer
            std::vector<double> q_current=q_init;

            // Run the simulation till the intersection point of the look-ahead circle and piecewise linear path
            // reaches the last point along the piecewise linear path
            bool found_intersection_flag=true;

            // Current sim segment is different from the piecewise linear since it does not store the direction boolean flag
            // and is only one segment of the piecewise linear
            std::vector<std::vector<double>> current_sim_segment(2,std::vector<double>(2,0));

            for (int i=0; i<piecewise_linear.size()-1; i++){

                std::cout << "Simulating segment: " << i << std::endl;

                // Select the correct segment from the piecewise linear path for simulation
                current_sim_segment[0][0] = piecewise_linear[i][0];
                current_sim_segment[0][1] = piecewise_linear[i][1];
                current_sim_segment[1][0] = piecewise_linear[i+1][0];
                current_sim_segment[1][1] = piecewise_linear[i+1][1];

                if(i==0){
                    final_trajectory = segmentSimulator(q_current, current_sim_segment, piecewise_linear[i+1][2]);
                }
                else{
                    temp_trajectory = segmentSimulator(q_current, current_sim_segment, piecewise_linear[i+1][2]);
                    final_trajectory.insert(final_trajectory.end(), temp_trajectory.begin(), temp_trajectory.end());
                }

                q_current = final_trajectory[final_trajectory.size()-1];

            }

            return final_trajectory;

        }

};