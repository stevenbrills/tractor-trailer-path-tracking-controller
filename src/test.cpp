#include <vector>
#include <stdexcept>
#include <cmath>
#include <iostream>
#include <ct/optcon/optcon.h>
#include "TractorTrailerController.hpp"


int main(){

    // Instantiate tractor-trailer simmulator object
    // Parameters in brackets are: tractor wheelbase, trailer wheelbase, tractor rear axle hitch offset, forward lookahead radius
    // backward lookahead radius, velocity, Q and R controller gains
    TractorTrailerController simulator(1.5, 5, 0.5, 4, 5, 0.3, 1, 1);

    // Create a starting state for the vehicle
    std::vector<double> q_init(6,0);
    q_init[0] = 8;
    q_init[1] = 22.5;
    q_init[2] = M_PI;
    q_init[3] = M_PI;

    // Create a piecewise linear path for tracking starting at tractor axle center
    std::vector<std::vector<double>> path(6, std::vector<double>(3,0));

    path[0][0] = 8;     // x-coordinate of point
    path[0][1] = 22.5;  // y-coordinate of point
    path[0][2] = 0;     // Flag for forward or backward motion (0=reverse, 1=forward)

    path[1][0] = 39;
    path[1][1] = 19;
    path[1][2] = 0;

    path[2][0] = 40;
    path[2][1] = 2;
    path[2][2] = 0;

    path[3][0] = 40.5;
    path[3][1] = 23;
    path[3][2] = 1;

    path[4][0] = 50;
    path[4][1] = 26;
    path[4][2] = 1;

    path[5][0] = 70;
    path[5][1] = 22.5;
    path[5][2] = 1;

    // Simulate the trajectory by calling the forwardSimulator method
    std::vector<std::vector<double>> traj = simulator.forwardSimulator(q_init, path);

    // Save the computed trajectory into a text file
	std::ofstream trajectory_file;
	trajectory_file.open("../output/control_test1.txt", std::ios::trunc);
	if (!trajectory_file.is_open()) {
		throw std::runtime_error("Cannot open file");
	}

    // File Description
	trajectory_file << "Trajectory File" << std::endl;

    // Tractor and Trailer Parameters
	trajectory_file << "Forward lookahead radius: " << simulator.getForwardLookaheadRadius() << 
    "Backward lookahead radius: " << simulator.getBackwardLookaheadRadius() << "Tractor Wheelbase: " 
    << simulator.getTractorWheelbase() << "Trailer Wheelbase: " << simulator.getTrailerWheelbase() << 
    "Tractor hitch offset: " << simulator.getTractorHitchOffset() << "Velocity" << simulator.getVelocity() 
    << std::endl;

    trajectory_file << std::endl;

    // Save the piecewise linear path
    trajectory_file << "Piecewise Linear Path" << std::endl;
    for (auto& segment : path){
        trajectory_file << segment[0] << " " << segment[1] << std::endl;
    }

    // Save the trajectory
    trajectory_file << "Trajectory" << std::endl;
    for (auto& state : traj){
        trajectory_file << state[0] << " " << state[1] << " " << state[2] << " " <<
        (M_PI - state[3]) << " " << state[4] << " " << state[5] << " " << state[6] << std::endl;
    }

    return 0;

}