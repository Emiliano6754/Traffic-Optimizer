#ifndef PERIODIC_STATEH
#define PERIODIC_STATEH

#include<Eigen/Dense>

struct driver_params {
    const unsigned int &n_drivers;
    Eigen::VectorXd alpha, beta, tau;

    // Creates driver_params with n_driv drivers, all parameters filled with zeros
    driver_params(const unsigned int &n_driv);

    // Creates driver_params with fixed input parameters
    driver_params(const unsigned int &n_driv, const Eigen::VectorXd &agilities, const Eigen::VectorXd &anticipations, const Eigen::VectorXd &safety_times);
};

struct periodic_state {
    const double length, time_delta;
    Eigen::VectorXd positions, velocities, accelerations;
    
    // Creates periodic_state with all positions, velocities and accelerations initialized to zero. Should probably be changed to initialized between some random values. Or create another for that
    periodic_state(const unsigned int &n_drivers, const double &length, const double &time_delta);

    // Creates periodic_state with specified initial positions, velocities and accelerations
    periodic_state(const unsigned int &n_drivers, const double &length, const double &time_delta, const Eigen::VectorXd &pos, const Eigen::VectorXd &vel, const Eigen::VectorXd &acc);
};

// Evolves periodic state for a time step given by state.time_delta. Calculates the mod of the positions on each step, considering the periodic state 
void step_state_mod(periodic_state &state, const driver_params &driver_params, const double &length);

// Evolves periodic state for a time step given by state.time_delta. Positions are free to increase higher than state.length, to avoid modulus operations. Can't be used for extremely long time simulations
void step_state(periodic_state &state, const driver_params &driver_params);

// Evolves state for T steps and returns state variables for all the steps in the evolution as matrices, where fixed rows are fixed steps. Unless pos_modulus is True, positions are free to increase higher than state.length, to avoid modulus operations. For large T simulations, pos_modulus should be True
void evolve_state(const unsigned int &T, periodic_state &state, const driver_params &driver_params, Eigen::MatrixXd &state_positions, Eigen::MatrixXd &state_velocities, Eigen::MatrixXd &state_accelerations, const bool &pos_modulus);

#endif