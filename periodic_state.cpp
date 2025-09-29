#include "periodic_state.h"
#include<iostream>
#include<cmath>

using Eigen::VectorXd, Eigen::seq, Eigen::last, Eigen::fix;

// Checks if vec has size elements. If it doesn't, cut the elements after size or fill with zeros the missing values. Name is used to warn user of this
void check_correct_size(VectorXd &vec, const unsigned int &size, const std::string &name) {
    if (vec.size() != size) {
        std::cout << name << " vector is not the right size, cutting the last elements or filling them with zeros";
        int in_size = vec.size();
        // Should check if resize actually resizes in place
        vec.resize(size);
        if (in_size < size) {
            for (int j = in_size; j < size; j++) {
                vec[j] = 0;
            }
        }
    }
}


// Creates driver_params with n_driv drivers, all parameters filled with zeros
driver_params::driver_params(const unsigned int &n_driv) : n_drivers(n_driv) {
    alpha = VectorXd::Zero(n_driv), beta = VectorXd::Zero(n_driv), tau = VectorXd::Zero(n_driv);
}

// Creates driver_params with fixed input parameters
driver_params::driver_params(const unsigned int &n_driv, const VectorXd &agilities, const VectorXd &anticipations, const VectorXd &safety_times): n_drivers(n_driv), alpha(agilities), beta(anticipations), tau(safety_times) {
    check_correct_size(alpha, n_drivers, "alpha");
    check_correct_size(beta, n_drivers, "beta");
    check_correct_size(tau, n_drivers, "tau");
}
  
// Creates periodic_state with all positions, velocities and accelerations initialized to zero. Should probably be changed to initialized between some random values. Or create another for that
periodic_state::periodic_state(const unsigned int &n_drivers, const double &length, const double &time_delta) : length(length), time_delta(time_delta) {
    positions = VectorXd::Zero(n_drivers), velocities = VectorXd::Zero(n_drivers), accelerations = VectorXd::Zero(n_drivers);
}

// Creates periodic_state with specified initial positions, velocities and accelerations
periodic_state::periodic_state(const unsigned int &n_drivers, const double &length, const double &time_delta, const VectorXd &pos, const VectorXd &vel, const VectorXd &acc) : length(length), time_delta(time_delta), positions(pos), velocities(vel), accelerations(acc) {
    check_correct_size(positions, n_drivers, "positions");
    check_correct_size(velocities, n_drivers, "velocities");
    check_correct_size(accelerations, n_drivers, "accelerations");
}

// Evolves periodic state for a time step given by state.time_delta. Calculates the mod of the positions on each step, considering the periodic state 
void step_state_mod(periodic_state &state, const driver_params &driver_params, const double &length) {
    const VectorXd position_deltas = ( VectorXd(driver_params.n_drivers) << state.positions(seq(fix<1>, last)), state.positions(0) + length ).finished() - state.positions;
    const VectorXd lead_velocities = ( VectorXd(driver_params.n_drivers) << state.velocities(seq(fix<1>, last)), state.velocities(0) ).finished();
    state.accelerations = - (driver_params.alpha + driver_params.beta).cwiseProduct(state.velocities) + position_deltas.cwiseProduct(driver_params.alpha).cwiseQuotient(driver_params.tau) + driver_params.beta.cwiseProduct(lead_velocities);
    state.velocities += state.time_delta * state.accelerations;
    state.positions += state.time_delta * state.velocities;
    state.positions = state.positions.unaryExpr([length](const double &x) { return  x - length * std::floor(x / length); });
}

// Evolves periodic state for a time step given by state.time_delta. Positions are free to increase higher than state.length, to avoid modulus operations. Can't be used for extremely long time simulations
void step_state(periodic_state &state, const driver_params &driver_params) {
    const VectorXd position_deltas = ( VectorXd(driver_params.n_drivers) << state.positions(seq(fix<1>, last)), state.positions(0) + state.length ).finished() - state.positions;
    const VectorXd lead_velocities = ( VectorXd(driver_params.n_drivers) << state.velocities(seq(fix<1>, last)), state.velocities(0) ).finished();
    state.accelerations = - (driver_params.alpha + driver_params.beta).cwiseProduct(state.velocities) + position_deltas.cwiseProduct(driver_params.alpha).cwiseQuotient(driver_params.tau) + driver_params.beta.cwiseProduct(lead_velocities);
    state.velocities += state.time_delta * state.accelerations;
    state.positions += state.time_delta * state.velocities;
}

// Evolves state for T steps and returns state variables for all the steps in the evolution as matrices, where fixed rows are fixed steps. Unless pos_modulus is True, positions are free to increase higher than state.length, to avoid modulus operations. For large T simulations, pos_modulus should be True
void evolve_state(const unsigned int &T, periodic_state &state, const driver_params &driver_params, Eigen::MatrixXd &state_positions, Eigen::MatrixXd &state_velocities, Eigen::MatrixXd &state_accelerations, const bool &pos_modulus) {
    state_positions.resize(T, driver_params.n_drivers);
    state_velocities.resize(T, driver_params.n_drivers);
    state_accelerations.resize(T, driver_params.n_drivers);
    state_positions.row(0) = state.positions;
    state_velocities.row(0) = state.velocities;
    state_accelerations.row(0) = state.accelerations;
    if (pos_modulus) {
        for (unsigned int t = 1; t < T; t++) {
            step_state_mod(state, driver_params, state.length);
            state_positions.row(t) = state.positions;
            state_velocities.row(t) = state.velocities;
            state_accelerations.row(t) = state.accelerations;
        }
    } else {
        for (unsigned int t = 1; t < T; t++) {
            step_state(state, driver_params);
            state_positions.row(t) = state.positions;
            state_velocities.row(t) = state.velocities;
            state_accelerations.row(t) = state.accelerations;
        }
    }
}