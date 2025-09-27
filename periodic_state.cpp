#include "periodic_state.h"
#include<iostream>

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

struct driver_params {
    const unsigned int &n_drivers;
    VectorXd alpha, beta, tau;

    driver_params(const unsigned int &n_driv) : n_drivers(n_driv) {
        alpha = VectorXd::Zero(n_driv), beta = VectorXd::Zero(n_driv), tau = VectorXd::Zero(n_driv);
    }

    driver_params(const unsigned int &n_driv, const VectorXd &agilities, const VectorXd &anticipations, const double &safety_times): n_drivers(n_driv), alpha(agilities), beta(anticipations), tau(safety_times) {
        check_correct_size(alpha, n_drivers, "alpha");
        check_correct_size(beta, n_drivers, "beta");
        check_correct_size(tau, n_drivers, "tau");
    }
};

struct periodic_state {
    VectorXd positions, velocities, accelerations;
    
    periodic_state(const unsigned int &n_drivers) {
        positions = VectorXd::Zero(n_drivers), velocities = VectorXd::Zero(n_drivers), accelerations = VectorXd::Zero(n_drivers);
    }

    periodic_state(const unsigned int &n_drivers, const VectorXd &pos, const VectorXd &vel, const VectorXd &acc) : positions(pos), velocities(vel), accelerations(acc) {
        check_correct_size(positions, n_drivers, "positions");
        check_correct_size(velocities, n_drivers, "velocities");
        check_correct_size(accelerations, n_drivers, "accelerations");
    }
};

periodic_state evolve_state(periodic_state &state, const driver_params &driver_params) {
    state.accelerations(seq(0, Eigen::last-Eigen::fix<1>)) = - (driver_params.alpha + driver_params.beta).cwiseProduct(state.velocities) - driver_params.alpha.cwiseProduct(state.positions).cwiseQuotient(driver_params.tau);
}