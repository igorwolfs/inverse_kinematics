#ifndef MAIN_H
#define MAIN_H


#include "inverse_kin/inverse_kin.h"
#include "matrix/matrix.h"

int main()
{
    // Note: T_Z, R_Z, T_X, R_X
    Matrix robot_dh = Matrix(DOF, 4, 0);
    robot_dh(0, 0) = 67;
    robot_dh(0, 3) = PI_2;
    robot_dh(1, 1) = PI_2;
    robot_dh(1, 2) = 95;
    robot_dh(2, 1) = -PI_2;
    robot_dh(2, 2) = 108;
    robot_dh(3, 1) = PI_2;
    robot_dh(3, 3) = PI_2;
    //robot_dh(4,0) = 58;
    //robot_dh(5,0) = 58;

    Matrix Ttool = Matrix(4);
    Matrix Rz = move_matrix(PI, r_z);
    Matrix Ry = move_matrix(-PI_2, r_y);

    Ttool *= Rz;
    Ttool *= Ry;

    // Create inverse kinematics object
    inverse_kin IK = inverse_kin(robot_dh, Ttool);

    Matrix Traj = Matrix(4, 0, 0);
    Matrix traj_temp;
    Matrix FORWARD = IK.fkine(Ttool, tuple<int, int>{0, DOF});

    for (int i = 0; i < 101; i++)
    {
        traj_temp = FORWARD;
        if (i <= 100)
        {
            traj_temp(2, 3) = traj_temp(2, 3) + (10 * cos(2 * PI * i / 100) - 10);
            traj_temp(1, 3) = traj_temp(1, 3) + (10 * cos(2 * PI * i / 100) - 10);
            traj_temp(0, 3) = traj_temp(0, 3) + (abs(i - 50) - 50);
        }
        Traj = Traj.append(traj_temp, 1);
    }

    // apply inverse kinematics to the trajectory
    IK.follow_trajectory(Traj);
    IK.performance_check();

 
}

#endif