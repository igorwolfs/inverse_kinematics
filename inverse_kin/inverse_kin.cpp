#include "inverse_kin.h"


inverse_kin::inverse_kin(Matrix & robot, Matrix & T_tool, bool save_)//= true)
{
    this->dh = robot;
    Ttool = T_tool;
    dq = Matrix(DOF,1,0);
    q = Matrix(DOF,1,0);
    qt = Matrix(DOF,0,0);
    Trajectory = Matrix(4,0,0); //fkine(Ttool, tuple<int, int>{0,DOF});
    save = save_;
    inverse_method = method::damped_ls_GS;
};

// NOT CONVERGING
void inverse_kin::p_inv_jac(Matrix & jac, Matrix & dr)
{
    Matrix jac_pinv = jac.pinv_4();
    dq = jac_pinv * dr;
}
void inverse_kin::inv_jac(Matrix & jac, Matrix & dr)
{
    Matrix jac_inv = jac.inv_4();
    dq = jac_inv * dr;
}

void inverse_kin::t_jac(Matrix & jac, Matrix & dr)
{
    // Calculate factor alfa
    Matrix dr_T = dr.transpose();
    Matrix jac_T = jac.transpose();
    Matrix numerator = (dr_T * jac) * jac_T;
    numerator = (numerator * dr);
    Matrix denominator = (jac * jac_T) * dr;
    Matrix den_T = denominator.transpose();
    denominator = den_T * denominator;
    double long alfa = numerator(0,0) / denominator(0,0);

    // Transpose Jacobian multiplication
    dq = (jac_T * dr);
    dq = dq * alfa;
}

void inverse_kin::damped_ls(Matrix & jac, Matrix & dr)
{
    long double damping_factor = 0.0001;
    Matrix jac_T = jac.transpose();
    Matrix eye = Matrix(5);
    Matrix right_side = eye * (damping_factor*damping_factor);
    right_side = (jac * jac_T) + right_side;
    // BIG NOTE: if you're going to != 4 DOF's you need a matrix inversion function here that can invert not only 4 DOF's
    Matrix right_side_inv = right_side.inv_4(); 
    dq = (jac_T * right_side_inv) * dr;
}

// UNSTABLE?
void inverse_kin::damped_ls_GS(Matrix & jac, Matrix & dr)
{
    Matrix eye = Matrix(DOF);
    double damping_factor = 0.0001;
    Matrix jac_T = jac.transpose();
    Matrix right_side = eye * (damping_factor*damping_factor);
    Matrix A_eq = ((jac_T * jac) + right_side);
    Matrix b_eq = jac_T * dr;
    Matrix dq_temp = Matrix(DOF,1,0);

    for (int i=0; i < DOF; i++)
    {
        long double Ax_sum = 0;
        for (int j=0; j < DOF; j++)
        {
            Ax_sum += A_eq(i,j) * dq(j,0);
        }

        dq_temp(i,0) = (b_eq(i,0) - Ax_sum)/A_eq(i,i);
    }
    dq = dq_temp;
}

void inverse_kin::update_q()
{

    for (int i=0; i<DOF; i++)
    {
        q(i,0) += dq(i,0); // added a minus here
    }
}

void inverse_kin::reset_q()
{
    int cols_qt = qt.getCols();
    q = qt.slice(tuple<int, int>{0,DOF}, tuple<int, int>{cols_qt-1, cols_qt});
}

void inverse_kin::save_path()
{
    if (save)
    {
        qt = qt.append(q, 1);
        Trajectory = Trajectory.append(Tdes, 1);
        double diff = 0;

        Tfk = fkine(Ttool, tuple<int, int>{0,DOF});
        for (int i=0; i<3; i++)
        {
            diff += (Tfk(i,3)-Tdes(i,3))*(Tfk(i,3)-Tdes(i,3));
        }

        std::cout << "T_fk_des, number: "  << qt.getCols() << std::endl;
        std::cout << Tdes << std::endl;
        std::cout << "T_fk_real, number: " << qt.getCols() << std::endl;
        std::cout << Tfk << std::endl;

        diff = sqrt(diff); 
        std::cout << diff << std::endl;
    }

}


Matrix inverse_kin::jacobian()
{
    Matrix Sz = Matrix(4,4,0);
    Matrix eye = Matrix(4);
    Sz(0,1) = -1;
    Sz(1,0) = 1;

    // Inverse Matrix
    Matrix T_workspace = fkine(Ttool, tuple<int, int>{0,DOF});
    Matrix T_inv = T_workspace.inv_4();

    // Column calculation of Jacobian
    Matrix column1 = (T_inv * Sz) * T_workspace;
    column1 = return_column(column1);
    Matrix J = column1; // initialize Jacobian
    if (DOF == 1)
    {
        return J;
    }

    Matrix fkine_2_1 = den_hart(0, this->q(0,0));
    Matrix fkine_2_2 = fkine(Ttool, tuple<int, int>{1,DOF});
    Matrix column2 = (T_inv * fkine_2_1);
    column2 *= Sz;
    column2 *= fkine_2_2;
    column2 = return_column(column2);

    J = J.append(column2, 1); // update Jacobian
    if (DOF == 2)
    {
        return J;
    }

    Matrix fkine_3_1 = fkine(eye, tuple<int, int>{0,1});
    Matrix fkine_3_2 = fkine(Ttool, tuple<int, int>{2, DOF});
    Matrix column3 = (T_inv * fkine_3_1);
    column3 = column3 * Sz;
    column3 = column3 * fkine_3_2;
    column3 = return_column(column3);

    J = J.append(column3, 1);
    if (DOF == 3)
    {
        return J;
    }

    Matrix fkine_4_1 = fkine(eye, tuple<int, int>{0,2});
    Matrix fkine_4_2 = fkine(Ttool, tuple<int, int>{3, DOF});
    Matrix column4 = (T_inv * fkine_4_1);
    column4 *= Sz;
    column4 *= fkine_4_2;
    column4 = return_column(column4);
 
    J = J.append(column4, 1);
    if (DOF == 4)
    {
        return J;
    }

    Matrix fkine_5_1 = fkine(eye, tuple<int, int>{0,3});
    Matrix fkine_5_2 = fkine(Ttool, tuple<int, int>{4, DOF});
    Matrix column5 = (T_inv * fkine_5_1);
    column5 *= Sz;
    column5 *= fkine_5_2;
    column5 = return_column(column5);

    J = J.append(column5, 1);
    if (DOF == 5)
    {
        return J;
    }

    Matrix fkine_6_1 = fkine(eye, tuple<int, int>{0,4});
    Matrix fkine_6_2 = fkine(Ttool, tuple<int, int>{5, DOF});
    Matrix column6 = (T_inv * fkine_6_1);
    column6 *= Sz;
    column6 *= fkine_6_2;
    column6 = return_column(column6);

    J = J.append(column6, 1);
    return J;
}

// Performs forward kinematics
Matrix inverse_kin::fkine(Matrix& Ttool, tuple<int,int> dh_slice = tuple<int,int>{0,DOF})
{
    Matrix T = Matrix(4);
    for (int i = get<0>(dh_slice); i < get<1>(dh_slice); i++)
    {
        Matrix DH_temp = den_hart(i, q(i,0));
        T = T * DH_temp;
    }
    T = T * Ttool;

    return T;
}


// shows the end effector position after the translation/rotation emposed.
/*
Sequence:
    > T_z
    > R_Z
    > T_x
    > R_x
*/
Matrix inverse_kin::den_hart(int link, double q_i)
{
    double dh0 = dh(link, 0);
    double dh1 = dh(link, 1);
    double dh2 = dh(link, 2);
    double dh3 = dh(link, 3);
    Matrix m_t_z = move_matrix(dh0, t_z);
    Matrix m_r_z = move_matrix((dh1 + q_i), r_z);
    Matrix m_t_x = move_matrix(dh2, t_x);
    Matrix m_r_x = move_matrix(dh3, r_x);
    Matrix T = m_t_z * m_r_z;
    T *= m_t_x;
    T *= m_r_x;
    return T;
}

// Arbitrary rotation about axes are allowed depending on the DOF's.
Matrix inverse_kin::return_column(Matrix & A)
{
    Matrix dr = Matrix(DOF,1,0);

    for (int i = 0; i < DOF; i++)
    {
        // Translational components
        if (i < 3)
        {
            dr(i, 0) = A(i, 3);
        }
        // Rotational components
        if (i == 3)
        {
            dr(i,0) = 0.5*(A(2,1)-A(1,2));
        }
        if (i == 4)
        {
            dr(i,0) = 0.5*(A(0,2)-A(2,0));
        }
        if (i == 5)
        {
            dr(i,0) = 0.5*(A(1,0)-A(0,1));
        }
    }
    return dr;
}


Matrix move_matrix(double l, enum move dir)
{
    Matrix T = Matrix(4);
    switch(dir)
    {
        case (t_x):
            T(0,3) = (abs(l)<MAX_STEP) ? l : SGN(l)*MAX_STEP;
            break;
        case(t_y):
            T(1,3) = (abs(l)<MAX_STEP) ? l : SGN(l)*MAX_STEP;
            break;

        case(t_z):
            T(2,3) = (abs(l)<MAX_STEP) ? l : SGN(l)*MAX_STEP;
            break;
        case(r_x):
            l = (abs(l)<MAX_ANGLE_STEP) ? l : SGN(l)*MAX_ANGLE_STEP;
            T(1,1) = cos(l);
            T(1,2) = -sin(l);
            T(2,1) = sin(l);
            T(2,2) = cos(l);
            break;
        case(r_y):
            l = (abs(l)<MAX_ANGLE_STEP) ? l : SGN(l)*MAX_ANGLE_STEP;
            T(0,0) = cos(l);
            T(0,2) = sin(l);
            T(2,0) = -sin(l);
            T(2,2) = cos(l);
            break;

        case(r_z):
            l = (abs(l)<MAX_ANGLE_STEP) ? l : SGN(l)*MAX_ANGLE_STEP;
            T(0,0) = cos(l);
            T(0,1) = -sin(l);
            T(1,0) = sin(l);
            T(1,1) = cos(l);
            break;

    }
    return T;
}


void inverse_kin::online_invert(double step)
{
    std::cout << "ONLINE TRAJECTORY" << std::endl;

    char val = 'x';
    Matrix move_matrix_;

    Tdes = fkine(Ttool, tuple<int,int>{0,DOF});
    Tfk = Tdes;
    while (val != 'e')
    {
        std::cin >> val;
        switch(val)
        {
            case('2'):
            {
                std::cout << "- X" << std::endl;                
                Tdes(0,3) -= step;
                //move_matrix_ = move_matrix(-step, t_x);
                //Tdes =  Tfk * move_matrix_;
                break;
            }
            case('8'):
            {
                Tdes(0,3) += step;
                std::cout << "+ X" << std::endl;
                //move_matrix_ = move_matrix(step, t_x);
                //Tdes =  Tfk * move_matrix_;
                break;
            }
            case('4'):
            {
                Tdes(1,3) += step;
                std::cout << "+ Y" << std::endl;
                //move_matrix_ = move_matrix(step, t_y);
                //Tdes =  Tfk * move_matrix_;
                break;
            }
            case('6'):
            {
                Tdes(1,3) += -step;
                std::cout << "- Y" << std::endl;
                //move_matrix_ = move_matrix(-step, t_y);
                //Tdes =  Tfk * move_matrix_;

                break;
            }
            case('5'):
            {
                Tdes(2,3) += step;
                std::cout << "+ Z" << std::endl;
                //move_matrix_ = move_matrix(step, t_z);
                //Tdes =  Tfk * move_matrix_;
                break;
            }
            case('0'):
            {         
                Tdes(2,3) -= step;   
                std::cout << "- Z" << std::endl;
                //move_matrix_ = move_matrix(-step, t_z);
                //Tdes =  Tfk * move_matrix_;
                break;
            }
        }

        

        update();
    }
}

void inverse_kin::follow_trajectory(Matrix& Traj_new)
{
    Matrix Tdes_seq;
    Matrix Tdes_end;
    int q_cols = Traj_new.getCols()/4;
    std::cout << "amount of steps in trajectory: " << q_cols << std::endl;
    
    Tfk = fkine(Ttool, tuple<int, int>{0,DOF});
    for (int i=0; i < (q_cols); i++)
    {
        std::cout << "i: " << i << std::endl;

        Tdes_end = Traj_new.slice(tuple<int, int>{0,4}, tuple<int, int>{i*4,(i+1)*4});

        // Insert your path-creating function from Tcurrent -> Tdes_end
        Tdes_seq = clamp_T(Tdes_end);
        
        // Calculate inverse kinematics for all elements of the path function
        int q_seq = Tdes_seq.getCols()/4;
        for (int path_i = 0; path_i < q_seq; path_i ++)
        {
            std::cout << "path_i: " << path_i << std::endl;
            Tdes = Tdes_seq.slice(tuple<int, int>{0,4}, tuple<int, int>{path_i*4,(path_i+1)*4});
            update();
        }
    }
}


void inverse_kin::update()
{
    Matrix A;
    Matrix T_inv;
    Matrix difference;
    Matrix eye = Matrix(4);
    Matrix J;
    double diff;
    double diff_old;

    int iter = 0;
    int attempt = 0;
    while (true)
    {
        if (iter > MAXITER)
        {
            std::cout << "\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n" << std::endl;
            attempt += 1;
            std::cout << "difference: " << diff << std::endl;
            std::cout << "Method: " << inverse_method << std::endl;

            int next_enum = (int(inverse_method) + 1)%5;
            if (attempt > 4)
            { 
                // NOTE: YOU CAN DO OTHER THINGS HERE AS WELL!
                save_path();
                throw(std::range_error("Maximum Amount of iterations reached, try again with another method."));     
            }
            inverse_method = static_cast<method>(next_enum);

            std::cout << "Switched to method: " << inverse_method << std::endl;
            reset_q();
            iter = 0;
            std::cout << "\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n" << std::endl;
        }

        Tfk = fkine(Ttool, tuple<int,int>{0,DOF});

        diff = 0;
        // Calculate difference between the desired and forward kinematics
        for (int i=0; i<3; i++)
        {
            diff += (Tfk(i,3)-Tdes(i,3))*(Tfk(i,3)-Tdes(i,3));
        }
        diff = sqrt(diff);

        if (diff < ACCURACY)
        {
            std::cout << "difference: " << diff << std::endl;
            save_path();
            std::cout << "Amount of iterations: " << iter << std::endl;
            break;
        }

        if (abs(diff_old-diff) < EPS)
        {
            std::cout << "\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n" << std::endl;
            std::cout << "Slow convergence: " << diff_old << " -> " << diff << std::endl;
            attempt += 1;
            std::cout << "difference: " << diff << std::endl;
            std::cout << "Method: " << inverse_method << std::endl;

            int next_enum = (int(inverse_method) + 1)%5;
            if (attempt > 4)
            { 
                // NOTE: YOU CAN DO OTHER THINGS HERE AS WELL!
                save_path();
                throw(std::range_error("Convergence too slow."));     
            }
            inverse_method = static_cast<method>(next_enum);

            std::cout << "Switched to method: " << inverse_method << std::endl;
            std::cout << "\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n" << std::endl;

            reset_q();
            iter = 0;
        }

        diff_old = diff;
        iter += 1;

        standard_update();
    }
}

void inverse_kin::standard_update()
{
    Matrix eye = Matrix(4);
    Matrix T_inv = Tfk.inv_4();
    Matrix A = (T_inv * Tdes) - eye;
    Matrix dr = return_column(A);
    Matrix J = jacobian();

    switch(inverse_method)
    {
        case method::damped_ls_GS:
        {
            damped_ls_GS(J, dr);
            break;
        }
        case method::damped_ls:
        {
            damped_ls(J, dr);
            break;
        }
        case method::t_jac:
        {
            t_jac(J, dr); 
            break;
        }
        case method::inv_jac:
        {
            inv_jac(J, dr);
            break;  
        }
        case method::p_inv_jac:
        {
            p_inv_jac(J, dr);
            break; 
        }
    }
    update_q(); // adds dq to q
}


Matrix inverse_kin::clamp_T(Matrix& Tdes_end)
{

    Matrix diff =  Tdes_end - Tfk;
    Matrix diff_trans = diff.slice(tuple<int,int>{0,3},tuple<int, int>{3,4});
    double diff_norm = diff_trans.norm();

    if (diff_norm < EPS)
    {
        return Tdes_end;
    }

    double counter = diff_norm;
    double diff_norm_temp;
    double ratio = 0;


    Matrix Tdes_seq = Matrix(4,0,0);
    Matrix Tdes_step;


    while (true)
    {
        diff_norm_temp = (counter > CLAMP) ? CLAMP : counter;
        std::cout << "diff_norm: " << diff_norm_temp << std::endl;
        if (counter < EPS)
        {
            break;
        }

        ratio += (diff_norm_temp/diff_norm);
        Tdes_step = diff * ratio + Tfk;
        Tdes_seq = Tdes_seq.append(Tdes_step, 1);

        // Update counter
        counter -= diff_norm_temp;
        std::cout << "counter: " << counter << std::endl;
    }
    return Tdes_seq;
}

// Check inverse kinematics
void inverse_kin::performance_check()
{
    if (!save)
    {
        throw std::invalid_argument("Performance check impossible without saved trajectory!");
    }
    int q_cols = qt.getCols();
    int Trajectory_cols = Trajectory.getCols();
    if (q_cols == 0)
    {
        throw std::invalid_argument("Inverse kinematics not performed yet!");
    }

    std::cout << "Returning information on " << q_cols << " steps" << std::endl;

    Matrix T_fk_real;
    Matrix T_fk_des;
    Matrix residual = Matrix(1, q_cols, 0);
    double total_residual = 0;

    for (int i = 0; i < (q_cols-2); i++)
    {
        q = qt.slice(tuple<int, int>{0,DOF}, tuple<int, int>{i,i+1});
        T_fk_des = Trajectory.slice(tuple<int, int>{0,4}, tuple<int, int>{4*i,4*(i+1)});
        T_fk_real = fkine(Ttool, tuple<int, int>{0,DOF});
        
        std::cout << "Tfk_des: " << i << T_fk_des << std::endl;
        std::cout << "Tfk_real: " << i << T_fk_real << std::endl;

        double diff = 0;
        
        for (int i=0; i<3; i++)
        {
            diff += (T_fk_des(i,3)-T_fk_real(i,3))*(T_fk_des(i,3)-T_fk_real(i,3));
        }
        diff = sqrt(diff);

        residual(0,i) = diff;
        total_residual += diff;
    }
    std::cout << residual << std::endl;
    std::cout << "average deviation: " << (total_residual/q_cols) << std::endl;
}

Matrix inverse_kin::get_qt()
{
    return this->qt;
}
