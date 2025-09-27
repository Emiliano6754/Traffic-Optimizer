#include<Eigen/Dense>

struct driver_params {
    const unsigned int &n_drivers;
    Eigen::VectorXd alpha, beta, tau;

    // Creates driver_params with n_driv drivers, all parameters filled with zeros 
    driver_params(const unsigned int &n_driv);

    // Creates driver_params with fixed input parameters
    driver_params(const unsigned int &n_driv, const Eigen::VectorXd &agilities, const Eigen::VectorXd &anticipations, const double &safety_times);
};

struct periodic_state {
    Eigen::VectorXd positions, velocities, accelerations;

    periodic_state(const unsigned int &n_drivers, const Eigen::VectorXd &pos, const Eigen::VectorXd &vel, const Eigen::VectorXd &acc);
};
