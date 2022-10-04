# inverse_kinematics

Basic project solving inverse kinematics problems using 5 numerical methods
- p-inverted jacobian method
- inverted jacobian method
- transposed jacobian method
- least squares method
- damped least squares method

## How to use
1. [main.cpp] Create a matrix for the robotic arm containing the corresponding denavit hartenberg parameters.
2. [main.cpp] Create a tool matrix Ttool using rotation matrices.
3. [main.cpp] 
	- A: Create a trajectory for the end manipulator and pass it to inverse_kin::follow_trajectory
	- B: Adapt the manipulator online using the inverse_kin::online_invert-function.
