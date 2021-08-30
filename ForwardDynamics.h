#pragma once

#include <fstream>
#include <iostream>
#include <stdio.h>
#include <commonmath.h>
#include <convexoptimization.h>
#include <numericalpde.h>
#include <timemanagement.h>
#include <chrono>
#include <vector>


namespace ContinuumRobotLibrary {

//Independent Parameters
const double E = 207e9;
const double rad = 0.0009;
const double total_mass = 0.040;
const Vector3d g = -9.81 * Vector3d::UnitZ();
const double L = 0.40;
const double base_to_motor = 0.0518; //m
const double dt = 5e-2;
const double alpha = 0;
const double T = 25;
const int N = 200;
const int num_tendons = 2;
typedef Array<double, num_tendons, 1> ArrayNd;
const ArrayNd compliance = ArrayNd::Constant(1e-4);
const DiagonalMatrix<double, 3> Bse(0, 0, 0);
const DiagonalMatrix<double, 3> Bbt(5e-4, 5e-4, 5e-4);
const DiagonalMatrix<double, 3> C(1e-4, 1e-4, 1e-4);
const double tendon_offset = 0.019;
const Vector3d r[num_tendons] = { tendon_offset * Vector3d::UnitX(), tendon_offset * Vector3d::UnitY() };

const double zA = 0;          //   z(m)  Ramp Motor Input:
const double zB = -0.01619;   // zA ^ ____           ______
const double t1 = 3.652;      //    |     \         /
const double t2 = 3.967;      // zB |      \_______/
const double t3 = 17.63;      //    -------------------------> t(s)
const double t4 = 17.94;      //         t1 t2    t3 t4

double zX(double t) {
    if (t > t1 && t <= t2)       return zA + (zB - zA) * (t - t1) / (t2 - t1); //Ramp lower
    else if (t > t2 && t <= t3)  return zB;                          //Low state
    else if (t > t3 && t <= t4)  return zB + (zA - zB) * (t - t3) / (t4 - t3); //Ramp higher
    else                    return zA;                          //High state
}

std::vector<double> q_val = { 0, 0 };
double z(double t) {
    return q_val[0];
}

double z1(double t) {
    return q_val[0];
}

double z2(double t) {
    return q_val[1];
}

//Dependent parameter calculations
const double c0 = TimeManagerBdfAlpha::getC0(dt, alpha);
const double G = E / (2 * 1.3);
const double A = pi * pow(rad, 2);
const double I = pi * pow(rad, 4) / 4;
const double rho = total_mass / (L * A);
const DiagonalMatrix<double, 3> Kse(G* A, G* A, E* A);
const DiagonalMatrix<double, 3> Kbt(E* I, E* I, G * 2 * I);
const DiagonalMatrix<double, 3> J(I, I, 2 * I);
const Matrix3d Kse_dense = Kse.toDenseMatrix();
const Matrix3d Kbt_dense = Kbt.toDenseMatrix();
const Matrix3d Kse_c0Bse = Kse.toDenseMatrix() + c0 * Bse.toDenseMatrix();
const Matrix3d Kbt_c0Bbt = Kbt.toDenseMatrix() + c0 * Bbt.toDenseMatrix();
const Vector3d rhoAg = rho * A * g;
const double rhoA = rho * A;
const DiagonalMatrix<double, 3> rhoJ = rho * J;
static ArrayNd tau;
static double t = 0;

//ODE describing tendon robot statics
void cosseratTendonRobotOde(VectorXd& y_s_out, Map<VectorXd> z_out, Map<VectorXd> y) {
    //Unpack state vector
    Matrix3d R = Map<Matrix3d>(&y[3]);
    Vector3d v = Map<Vector3d>(&y[12]);
    Vector3d u = Map<Vector3d>(&y[15]);

    Vector3d a = Vector3d::Zero();
    Vector3d b = Vector3d::Zero();
    Matrix3d A_plus_Kse = Kse_dense;
    Matrix3d G = Matrix3d::Zero();
    Matrix3d H_plus_Kbt = Kbt_dense;

    Map<ArrayNd> pb_s_norm(&y_s_out[24], num_tendons);

    for (int i = 0; i < num_tendons; i++) {
        Vector3d pb_si = u.cross(r[i]) + v;
        pb_s_norm(i) = pb_si.norm();
        Matrix3d A_i = -hat_squared(pb_si) * (tau(i) / pow(pb_s_norm(i), 3));
        Matrix3d G_i = -hat_postmultiply(A_i, r[i]);
        Vector3d a_i = A_i * (u.cross(pb_si));

        a += a_i;
        b += r[i].cross(a_i);
        A_plus_Kse += A_i;
        G += G_i;
        H_plus_Kbt += hat_premultiply(r[i], G_i);
    }

    Matrix6d K;
    K << A_plus_Kse, G, G.transpose(), H_plus_Kbt;

    Vector3d nb = Kse * (v - Vector3d::UnitZ());
    Vector3d mb = Kbt * u;

    Vector6d rhs;
    rhs << -u.cross(nb) - transposeMultiply(R, rhoAg) - a,
        -u.cross(mb) - v.cross(nb) - b;

    //Pack state vector derivative
    Map<Vector3d> p_s(&y_s_out[0]);
    Map<Matrix3d> R_s(&y_s_out[3]);
    Map<Vector6d> vs_and_us(&y_s_out[12]);
    Map<Vector3d> q_s = Map<Vector3d>(&y_s_out[18]);
    Map<Vector3d> w_s = Map<Vector3d>(&y_s_out[21]);

    //ODEs
    p_s = R * v;
    R_s = hat_postmultiply(R, u);
    vs_and_us = K.selfadjointView<Eigen::Upper>().llt().solve(rhs);
    q_s = Vector3d::Zero();
    w_s = Vector3d::Zero();

    //Output argument for variables with time derivatives
    z_out << v, u, Vector6d::Zero(), vs_and_us;
}

//PDE semi-discretized in time describing tendon backbone dynamics
void tendonBackbonePDE(VectorXd& y_s_out, Map<VectorXd> z_out, Map<VectorXd> y, Map<VectorXd> z_h) {
    //Unpack state vector
    Matrix3d R = Map<Matrix3d>(&y[3]);
    Vector3d v = Map<Vector3d>(&y[12]);
    Vector3d u = Map<Vector3d>(&y[15]);
    Vector3d q = Map<Vector3d>(&y[18]);
    Vector3d w = Map<Vector3d>(&y[21]);

    Vector3d a = Vector3d::Zero();
    Vector3d b = Vector3d::Zero();
    Matrix3d A_plus_Kse_c0Bse = Kse_c0Bse;
    Matrix3d G = Matrix3d::Zero();
    Matrix3d H_plus_Kbt_c0Bbt = Kbt_c0Bbt;

    Map<ArrayNd> pb_s_norm(&y_s_out[24], num_tendons);

    for (int i = 0; i < num_tendons; i++) {
        Vector3d pb_si = u.cross(r[i]) + v;
        pb_s_norm(i) = pb_si.norm();
        Matrix3d A_i = -hat_squared(pb_si) * (tau(i) / pow(pb_s_norm(i), 3));
        Matrix3d G_i = -hat_postmultiply(A_i, r[i]);
        Vector3d a_i = A_i * (u.cross(pb_si));

        a += a_i;
        b += r[i].cross(a_i);
        A_plus_Kse_c0Bse += A_i;
        G += G_i;
        H_plus_Kbt_c0Bbt += hat_premultiply(r[i], G_i);
    }

    Matrix6d K;
    K << A_plus_Kse_c0Bse, G, G.transpose(), H_plus_Kbt_c0Bbt;

    Map<Vector3d> v_h(&z_h[0]);
    Map<Vector3d> u_h(&z_h[3]);
    Map<Vector3d> q_h(&z_h[6]);
    Map<Vector3d> w_h(&z_h[9]);
    Map<Vector3d> v_sh(&z_h[12]);
    Map<Vector3d> u_sh(&z_h[15]);

    Vector3d v_t = c0 * v + v_h;
    Vector3d u_t = c0 * u + u_h;
    Vector3d q_t = c0 * q + q_h;
    Vector3d w_t = c0 * w + w_h;

    Vector3d nb = Kse * (v - Vector3d::UnitZ()) + Bse * v_t;
    Vector3d mb = Kbt * u + Bbt * u_t;

    Vector6d rhs;
    rhs << -a + rhoA * (w.cross(q) + q_t) + C * q.cwiseProduct(q.cwiseAbs()) - transposeMultiply(R, rhoAg) - u.cross(nb) - Bse * v_sh,
        -b + w.cross(rhoJ * w) + rhoJ * w_t - v.cross(nb) - u.cross(mb) - Bbt * u_sh;

    //Pack state vector derivative
    Map<Vector3d> p_s = Map<Vector3d>(&y_s_out[0]);
    Map<Matrix3d> R_s = Map<Matrix3d>(&y_s_out[3]);
    Map<Vector6d> vs_and_us(&y_s_out[12]);
    Map<Vector3d> q_s = Map<Vector3d>(&y_s_out[18]);
    Map<Vector3d> w_s = Map<Vector3d>(&y_s_out[21]);

    //ODEs
    p_s = R * v;
    R_s = hat_postmultiply(R, u);
    vs_and_us = K.selfadjointView<Eigen::Upper>().llt().solve(rhs);
    q_s = v_t - u.cross(q) + w.cross(v);
    w_s = u_t - u.cross(w);

    //Output argument for variables with time derivatives
    z_out << v, u, q, w, vs_and_us;
}

static MatrixXd Y(24 + num_tendons, N), Z(18, N), Z_h(18, N);
template<bool is_dynamic>
VectorXd obj(VectorXd guess) {
    Vector3d v0 = Kse.inverse() * guess.head(3) + Vector3d::UnitZ(); //not exactly guessing n0 due to viscoelastic constitutive equation
    Vector3d u0 = guess.segment<3>(3);
    tau = guess.tail<num_tendons>().cwiseMax(0);
    ArrayNd slack = -(guess.tail<num_tendons>().cwiseMin(0));

    VectorXd y0(24 + num_tendons);
    y0 << 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, v0, u0, Vector6d::Zero(), ArrayNd::Constant(base_to_motor);
    Y.col(0) = y0;

    //Numerically integrate the Cosserat rod equations
    if (is_dynamic) TimeBdfAlpha_SpaceEuler<tendonBackbonePDE, 24 + num_tendons, 18, N>(Y, Z, y0, L, Z_h);
    else euler<cosseratTendonRobotOde, 24 + num_tendons, 18, N>(Y, Z, y0, L);

    //Find the internal forces in the backbone prior to the final plate
    Vector3d vL = Y.block<3, 1>(12, N - 1);
    Vector3d uL = Y.block<3, 1>(15, N - 1);
    Vector3d vL_t = is_dynamic ? (c0 * vL + Z_h.block<3, 1>(0, N - 1)).eval() : Vector3d::Zero();
    Vector3d uL_t = is_dynamic ? (c0 * uL + Z_h.block<3, 1>(3, N - 1)).eval() : Vector3d::Zero();

    Vector3d nb = Kse * (vL - Vector3d::UnitZ()) + Bse * vL_t;
    Vector3d mb = Kbt * uL + Bbt * uL_t;

    //Find the equilibrium error at the tip, considering tendon forces
    Vector3d force_error = -nb;
    Vector3d moment_error = -mb;
    for (int i = 0; i < num_tendons; i++) {
        Vector3d pb_si = uL.cross(r[i]) + vL;
        Vector3d Fb_i = -tau(i) * pb_si.normalized();
        force_error += Fb_i;
        moment_error += r[i].cross(Fb_i);
    }

    //Find the length violation error
    ArrayNd integrated_lengths = Y.block<num_tendons, 1>(24, N - 1);
    //ArrayNd l_star = ArrayNd::Constant(L + base_to_motor + z(t));
    //(ArrayNd() << L + base_to_motor + z1(t), L + base_to_motor + z2(t));
    ArrayNd l_star;
    l_star << L + base_to_motor + z1(t), L + base_to_motor + z2(t);
    ArrayNd stretch = l_star * compliance * tau;
    ArrayNd length_error = integrated_lengths + slack - (l_star + stretch);

    VectorXd distal_error(6 + num_tendons);
    distal_error << force_error, moment_error, length_error;

    return distal_error;
}

class ForwardDynamics {
public:
    ForwardDynamics() :
        iter(0),
        guess(VectorXd::Zero(6 + num_tendons)),
        time_scheme(Z, Z_h, dt, alpha),
        M(static_cast<int>(T / dt)),
        tendon(3 * M, N + 1),
        tip(3, M),
        centerline(3 * M, N),
        disks(3 * M, 4 * 7)
    {
        guess = solveLevenbergMarquardt<obj<false> >(guess, 1e-12, 1500, 1e-2, 0.5, 1e-7, 1e-9);
        TimeManagerBdfAlpha time_scheme(Z, Z_h, dt, alpha);
    }

    void Next(std::vector<double> q) {
        q_val[0] = q[0];
        q_val[1] = q[1];
        guess = solveLevenbergMarquardt<obj<true> >(guess, 1e-4, 1500, 1e-2, 0.5, 1e-7, 1e-9);

        iter += 1;
        t += dt;
        time_scheme.advanceTime();

        //Store results
        tip.col(iter) = Y.block<3, 1>(0, N - 1);
        std::cout << q_val[1] << " > " << tip.col(iter).row(1) << std::endl;
        centerline.block<3, N>(3 * iter, 0) = Y.block<3, N>(0, 0);
        tendon.block<3, 1>(3 * iter, 0) = r[0] + z(t) * Vector3d::UnitZ();
        for (int j = 0; j < N; j++) {
            Matrix3d R = Map<Matrix3d>(&Y(3, j));
            tendon.block<3, 1>(3 * iter, j + 1) = centerline.block<3, 1>(3 * iter, j) + R * r[0];
        }
        for (int j = 1; j <= 7; j++) {
            int k = (N - 1) * j / 7;
            disks.block<3, 1>(3 * iter, 4 * (j - 1) + 3) = centerline.block<3, 1>(3 * iter, k);
            disks.block<3, 3>(3 * iter, 4 * (j - 1)) = Map<Matrix3d>(&Y(3, k));
        }
    }

    void Save() {
        //Save results for Blender visualization
        std::fstream file("tendon.dat", std::fstream::out);
        file << tendon;
        file.close();

        //Save results for Matlab visualization
        file = std::fstream("tip.dat", std::fstream::out);
        file << tip;
        file.close();

        file = std::fstream("disks.dat", std::fstream::out);
        file << disks;
        file.close();

        file = std::fstream("centerline.dat", std::fstream::out);
        file << centerline;
        file.close();
    }

private:
    
    //VectorXd minimization(VectorXd guess, bool is_dynamic);
    //void cosseratTendonRobotOde(VectorXd& y_s_out, Map<VectorXd> z_out, Map<VectorXd> y);
    //template<bool is_dynamic> VectorXd func(VectorXd guess);
    /*
    const double c0;
    const double G;
    const double A;
    const double I;
    const double rho;
    const DiagonalMatrix<double, 3> Kse;
    const DiagonalMatrix<double, 3> Kbt;
    const DiagonalMatrix<double, 3> J;
    const Matrix3d Kse_dense;
    const Matrix3d Kbt_dense;
    const Matrix3d Kse_c0Bse;
    const Matrix3d Kbt_c0Bbt;
    const Vector3d rhoAg;
    const double rhoA;
    const DiagonalMatrix<double, 3> rhoJ;
    static ArrayNd tau;
    double t;
    */

    int iter;
    VectorXd guess;
    TimeManagerBdfAlpha time_scheme;
    int M;
    MatrixXd tendon;
    Matrix3Xd tip;
    MatrixXd centerline;
    MatrixXd disks;

    //MatrixXd Y;
    //MatrixXd Z;
    //MatrixXd Z_h;

};
/*
VectorXd guess = VectorXd::Zero(6 + num_tendons);

int inc = 0;
int M = static_cast<int>(T / dt);
MatrixXd tendon(3 * M, N + 1);
Matrix3Xd tip(3, M);
MatrixXd centerline(3 * M, N);
MatrixXd disks(3 * M, 4 * 7);

void Init() {
    guess = solveLevenbergMarquardt<obj<false> >(guess, 1e-12, 1500, 1e-2, 0.5, 1e-7, 1e-9);


    //Solve dynamic
    TimeManagerBdfAlpha time_scheme(Z, Z_h, dt, alpha);

    for (int i = 0; i < M; i++) {
        guess = solveLevenbergMarquardt<obj<true> >(guess, 1e-12, 1500, 1e-2, 0.5, 1e-7, 1e-9);
        //std::cout << t << std::endl;

        t += dt;
        time_scheme.advanceTime();

        //Store results
        tip.col(i) = Y.block<3, 1>(0, N - 1);
        std::cout << tip.col(i).row(1) << std::endl;
        centerline.block<3, N>(3 * i, 0) = Y.block<3, N>(0, 0);
        tendon.block<3, 1>(3 * i, 0) = r[0] + z(t) * Vector3d::UnitZ();
        for (int j = 0; j < N; j++) {
            Matrix3d R = Map<Matrix3d>(&Y(3, j));
            tendon.block<3, 1>(3 * i, j + 1) = centerline.block<3, 1>(3 * i, j) + R * r[0];
        }
        for (int j = 1; j <= 7; j++) {
            int k = (N - 1) * j / 7;
            disks.block<3, 1>(3 * i, 4 * (j - 1) + 3) = centerline.block<3, 1>(3 * i, k);
            disks.block<3, 3>(3 * i, 4 * (j - 1)) = Map<Matrix3d>(&Y(3, k));
        }
    }
}


void Next(std::vector<double> q) {
    guess = solveLevenbergMarquardt<obj<true> >(guess, 1e-12, 1500, 1e-2, 0.5, 1e-7, 1e-9);
    //std::cout << t << std::endl;

    t += dt;
    time_scheme.advanceTime();

    //Store results
    tip.col(inc) = Y.block<3, 1>(0, N - 1);
    std::cout << tip.col(inc).row(1) << std::endl;
    centerline.block<3, N>(3 * inc, 0) = Y.block<3, N>(0, 0);
    tendon.block<3, 1>(3 * inc, 0) = r[0] + z(t) * Vector3d::UnitZ();
    for (int j = 0; j < N; j++) {
        Matrix3d R = Map<Matrix3d>(&Y(3, j));
        tendon.block<3, 1>(3 * inc, j + 1) = centerline.block<3, 1>(3 * inc, j) + R * r[0];
    }
    for (int j = 1; j <= 7; j++) {
        int k = (N - 1) * j / 7;
        disks.block<3, 1>(3 * inc, 4 * (j - 1) + 3) = centerline.block<3, 1>(3 * inc, k);
        disks.block<3, 3>(3 * inc, 4 * (j - 1)) = Map<Matrix3d>(&Y(3, k));
    }

    inc += 1;
}

void Save() {
    //Save results for Blender visualization
    std::fstream file("tendon.dat", std::fstream::out);
    file << tendon;
    file.close();

    //Save results for Matlab visualization
    file = std::fstream("tip.dat", std::fstream::out);
    file << tip;
    file.close();

    file = std::fstream("disks.dat", std::fstream::out);
    file << disks;
    file.close();

    file = std::fstream("centerline.dat", std::fstream::out);
    file << centerline;
    file.close();
}

int Run() {
    //Solve static

    //Solve dynamic
    TimeManagerBdfAlpha time_scheme(Z, Z_h, dt, alpha);

    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

    for (int i = 0; i < M; i++) {
        guess = solveLevenbergMarquardt<obj<true> >(guess, 1e-12, 1500, 1e-2, 0.5, 1e-7, 1e-9);
        //std::cout << t << std::endl;

        t += dt;
        time_scheme.advanceTime();

        //Store results
        tip.col(i) = Y.block<3, 1>(0, N - 1);
        std::cout << tip.col(i).row(1) << std::endl;
        centerline.block<3, N>(3 * i, 0) = Y.block<3, N>(0, 0);
        tendon.block<3, 1>(3 * i, 0) = r[0] + z(t) * Vector3d::UnitZ();
        for (int j = 0; j < N; j++) {
            Matrix3d R = Map<Matrix3d>(&Y(3, j));
            tendon.block<3, 1>(3 * i, j + 1) = centerline.block<3, 1>(3 * i, j) + R * r[0];
        }
        for (int j = 1; j <= 7; j++) {
            int k = (N - 1) * j / 7;
            disks.block<3, 1>(3 * i, 4 * (j - 1) + 3) = centerline.block<3, 1>(3 * i, k);
            disks.block<3, 3>(3 * i, 4 * (j - 1)) = Map<Matrix3d>(&Y(3, k));
        }
    }

    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();

    std::cout << "Speed = " << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() / M << "[us]" << std::endl;




    //Show the end-effector trajectory
#ifdef QT_CORE_LIB
    plot(VectorXd::LinSpaced(M, 0, t), tip.row(0), "Tip Trajectory", "t (s)", "x (m)");
#endif

    return 0;
}
*/

}