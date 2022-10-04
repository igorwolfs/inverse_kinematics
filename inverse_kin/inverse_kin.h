#ifndef INVERSE_KIN_H
#define INVERSE_KIN_H

#include "../matrix/matrix.h"
#include <cmath>
#include <iostream>
#include <tuple>


#define DOF 4
#define MAX_STEP 1000
#define MAX_ANGLE_STEP 2*PI
#define CLAMP 70
#define EPS 0.0000001 // Value that checks when something is close to 0
#define ACCURACY 0.1 // Value that determines the accuracy of the inverse kinematics for X, Y and Z coordinates
#define MAXITER 1000 // maximum amount of iterations for the inverse kinematics algorithm

#define SGN(val) ((0 < val) - (val < 0))


const long double PI = 3.14159265358979323846264338327950288419716939937510582097494459230781640628620899862803482534211706798214808;
const long double PI_2 = PI/2;

// Define enum
enum move {t_x, t_y, t_z, r_x, r_y, r_z};
enum method {damped_ls_GS, damped_ls, t_jac, inv_jac, p_inv_jac};

using std::tuple;
using std::get;

class inverse_kin {
    private:
    Matrix dh;
    Matrix Ttool;
    Matrix Tfk;
    Matrix Tdes;
    Matrix Trajectory;
    Matrix dq;
    Matrix q;
    Matrix qt;

    bool save;

    public:
    method inverse_method;

    // constructor & destructor
    inverse_kin(Matrix &, Matrix & , bool save = true);
    ~inverse_kin(){};

    // Jacobian, forward kinematics, DH

    Matrix jacobian();
    Matrix fkine(Matrix &, tuple<int,int>);
    Matrix den_hart(int, double);
    void update_q();
    void reset_q();


    // Inverse kinematics methods
    void p_inv_jac(Matrix &, Matrix &);
    void inv_jac(Matrix &, Matrix &);
    void t_jac(Matrix &, Matrix &);
    void damped_ls(Matrix &, Matrix &);
    void damped_ls_GS(Matrix &, Matrix &);


    // Main update function
    void update();

    // Update dr and dq
    void standard_update();
    Matrix clamp_T(Matrix &);

    // helper functions
    Matrix return_column(Matrix &);
    void set_Ttool(Matrix & tool){Ttool = tool;};
    void save_path();

    // Completed methods
    void online_invert(double);
    void follow_trajectory(Matrix&);

    // performance check of the inverse kinematics method
    void performance_check();
    Matrix get_qt();
};

 // function to create matrix that moves/rotates the end effector in a certain direction.
Matrix move_matrix(double, enum move);

#endif