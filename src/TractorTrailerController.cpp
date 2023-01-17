#include <vector>
#include <stdexcept>
#include <cmath>
#include <iostream>
#include <ct/optcon/optcon.h>
// #include "matplotlibcpp.h"
// #include "planner.hpp"
#include "TractorTrailerController.hpp"


int main(){


    TractorTrailerController simulator(1, 3, 0.3, 2, 5, 0.3, 1, 1);

    // Create a reference piecewise linear path for tracking
    std::vector<std::vector<double>> path(6,std::vector<double>(2,0));
    path[0][0] = 8;
    path[0][1] = 22.5;

    path[1][0] = 42;
    path[1][1] = 20;

    path[2][0] = 39;
    path[2][1] = 5;

    path[3][0] = 37;
    path[3][1] = 18;

    path[4][0] = 50;
    path[4][1] = 23;

    path[5][0] = 70;
    path[5][1] = 22.5;

    std::vector<double> q_init(6,0);
    q_init[0] = 8;
    q_init[1] = 22.5;
    q_init[2] = M_PI;
    q_init[3] = M_PI;

    std::vector<double> q_goal(6,0);
    q_goal[0] = 70;
    q_goal[1] = 22.5;
    q_goal[2] = M_PI;
    q_goal[3] = M_PI;

    auto trajectory = simulator.forwardSimulator(q_init, path);
    std::cout << "Completed Simulation over path!" << std::endl;

    return 0;

}

//