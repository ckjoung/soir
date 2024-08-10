/**@page Overview
 * This document provides documentation for the simple program convexity checker, and accompanies the paper
 *
 *
 *
 * The purpose of the program is to validate that the Gauss-Kronecker curvature of bounded component of the Levi-Civita
 * regularized Hamiltonian \f$H\f$ is positive for the values mentioned in the paper. Put in more direct language, this means
 * the following. Write \f$\Sigma:=H^{-1}(0)\f$ for the level set.
 * The program then tries to verify whether $Hess\, H|_{T \Sigma}$ is a positive definite 3 by 3-matrix.
 *
 * We use the following conventions for this Hamiltonian:
 * \f[
 * H=\frac{1}{2}\Vert w \Vert^2 + c \Vert z \Vert^2 -\frac{1-\mu}{2}+2\Vert z \Vert^2 (z_1 w_2-z_2w_1)
 * -\mu(z_1w_2+z_2w_1) -\frac{\Vert z\Vert^2}{ |2z^2-1|}.
 * \f]
 *
 * Here \f$\mu\f$ is the mass ratio, and
 * \f$c\f$ is the Jacobi energy. To clarify our conventions: the rotating Kepler problem has \f$\mu=0\f$, and
 * critical energy of the rotating Kepler problem is given by \f$c=\frac{3}{2}\f$.
 *
 *
 * We have computed the derivatives (gradient and Hessian) with Maple (CAS) and verified by hand.
 * Automatic differentiation would be cleaner since this is already present in CAPD,
 * but this manual approach is more direct and more transparent for checking.
 *
 * This code is NOT efficient (no attempt was made to reduce the dependency effect, for example
 * eg. abs in front squaring would help)
 *
 * In the paper we order the coordinates as (z_1,z_2,w_1,w_2).
 *
 * In the code we order the coordinates as (w_1,w_2,z_1,z_2).
 *
 */


#include <iostream>
#include <fstream>
#include <vector>

#include "capd/capdlib.h"
using namespace capd;
using namespace std;

/** Basic template version of the squaring function.
 *
 * @tparam T
 * @param a
 * @return
 */
template <typename T>
inline T sq(T a) {
    return a * a;
}


/** Basic operator overloading for printing vectors.
 *
 * @tparam T
 * @param out
 * @param v
 * @return
 */
template <typename T>
std::ostream &operator<<(std::ostream &out, std::vector<T> &v) {
    for( int i = 0; i < v.size(); ++i ) {
        out << v[i] << " ";
    }//
    return out;
}

/**
 * We use global variables since this is a short script only and this should improve
 * readability.
 */

Interval mu = Interval(0.0, 0.01);
Interval c = Interval(2.1, 2.1+1e-6);

/**
 * This function computes bounds on the values of the Levi-Civita Hamiltonian on the interval vector (box)
 * wz=(w_1,w_2;z_1,z_2).
 * @param wz
 * @return
 */
Interval LC_Hamiltonian(IVector &wz) {
    Interval w1 = wz[0];
    Interval w2 = wz[1];
    Interval z1 = wz[2];
    Interval z2 = wz[3];

    Interval t1 = w1*w1;
    Interval t3 = w2*w2;
    Interval t5 = z1*z1;
    Interval t6 = z2*z2;
    Interval t7 = t5+t6;
    Interval t10 = w1*z2;
    Interval t11 = z1*w2;
    Interval t20 = sq(2.0*t5-2.0*t6-1.0);
    Interval t23 = 16.0*t5*t6+t20;
    Interval t24 = sqrt(t23);
    Interval t27 = t1/2.0+t3/2.0+c*t7-1.0/2.0+mu/2.0+2.0*t7*(-t10+t11)-mu*(t10+t11)-mu*t7/t24;
    return t27;
}


/**
 *  This function computes the gradient of the Levi-Civita Hamiltonian
 *  The inner part was generated with Maple (CAS) and then checked and adapted.
 *
 * @param wz
 * @param result
 */
void LC_grad(IVector &wz, IVector &result) {
    Interval w1 = wz[0];
    Interval w2 = wz[1];
    Interval z1 = wz[2];
    Interval z2 = wz[3];

    result[0] = w1-(2.0*z1*z1+2.0*z2*z2)*z2-mu*z2;
    result[1] = w2+(2.0*z1*z1+2.0*z2*z2)*z1-mu*z1;

    Interval t8 = z1*z1;
    Interval t9 = z2*z2;
    Interval t10 = t8+t9;
    Interval t16 = 2.0*t8-2.0*t9-1.0;
    Interval t17 = t16*t16;
    Interval t20 = 16.0*t8*t9+t17;
    Interval t21 = sqrt(t20);
    Interval t36 = 2.0*c*z1+4.0*z1*(-w1*z2+z1*w2)+2.0*t10*w2-mu*w2-2.0*mu*z1/t21+mu*t10/t21/t20*(8.0*t16*z1+32.0*z1*t9)/2.0;
    result[2] = t36;

    t8 = z1*z1;
    t9 = z2*z2;
    t10 = t8+t9;
    t16 = 2.0*t8-2.0*t9-1.0;
    t17 = t16*t16;
    t20 = 16.0*t8*t9+t17;
    t21 = sqrt(t20);
    t36 = 2.0*c*z2+4.0*z2*(-w1*z2+z1*w2)-2.0*t10*w1-mu*w1-2.0*mu*z2/t21+mu*t10/t21/t20*(-8.0*t16*z2+32.0*t8*z2)/2.0;
    result[3] = t36;
}

/**
 *  This function computes the Hessian (4 by 4) of the Levi-Civita Hamiltonian
 *  The inner part was generated with Maple (CAS) and then checked and adapted.
 * @param wz
 * @return
 */
IMatrix LC_Hessian(IVector &wz) {
    Interval w1 = wz[0];
    Interval w2 = wz[1];
    Interval z1 = wz[2];
    Interval z2 = wz[3];

    IMatrix Hess(4,4);
    Hess[0][0] = 1.0;
    Hess[0][1] = 0.0;
    Hess[0][2] = -4.0*z1*z2;
    Hess[0][3] = -2.0*z1*z1-6.0*z2*z2-mu;
    Hess[1][1] = 1.0;
    Hess[1][2] = 6.0*z1*z1+2.0*z2*z2-mu;
    Hess[1][3] = 4.0*z1*z2;

    Interval t6 = z1*z1;
    Interval t8 = z2*z2;
    Interval t10 = 2.0*t6-2.0*t8-1.0;
    Interval t11 = t10*t10;
    Interval t14 = 16.0*t6*t8+t11;
    Interval t15 = sqrt(t14);
    Interval t21 = 1/t15/t14;
    Interval t26 = 8.0*t10*z1+32.0*z1*t8;
    Interval t31 = mu*(t6+t8);
    Interval t32 = t14*t14;
    Interval t35 = t26*t26;
    Interval t45 = 2.0*c-4.0*w1*z2+12.0*z1*w2-2.0*mu/t15+2.0*mu*z1*t21*t26-3.0/4.0*t31/t15/t32*t35+t31*t21*(48.0*t6+16.0*t8-8.0)/2.0;
    Hess[2][2] = t45;
    t6 = z1*z1;
    t8 = z2*z2;
    t10 = 2.0*t6-2.0*t8-1.0;
    t11 = t10*t10;
    t14 = 16.0*t6*t8+t11;
    t15 = sqrt(t14);
    Interval t17 = 1/t15/t14;
    Interval t22 = -8.0*t10*z2+32.0*t6*z2;
    Interval t30 = 8.0*t10*z1+32.0*z1*t8;
    Interval t34 = mu*(t6+t8);
    t35 = t14*t14;
    Interval t46 = -4.0*z1*w1+4.0*z2*w2+mu*z1*t17*t22+mu*z2*t17*t30-3.0/4.0*t34/t15/t35*t30*t22+16.0*t34*t17*z1*z2;
    Hess[2][3] = t46;
    t6 = z1*z1;
    t8 = z2*z2;
    t10 = 2.0*t6-2.0*t8-1.0;
    t11 = t10*t10;
    t14 = 16.0*t6*t8+t11;
    t15 = sqrt(t14);
    t21 = 1/t15/t14;
    t26 = -8.0*t10*z2+32.0*t6*z2;
    t31 = mu*(t6+t8);
    t32 = t14*t14;
    t35 = t26*t26;
    t45 = 2.0*c-12.0*w1*z2+4.0*z1*w2-2.0*mu/t15+2.0*mu*z2*t21*t26-3.0/4.0*t31/t15/t32*t35+t31*t21*(16.0*t6+48.0*t8+8.0)/2.0;
    Hess[3][3] = t45;
    for( int i=0; i < 4; i++ ) {
        for( int j = 0; j<i; j++)
            Hess[i][j] = Hess[j][i];
    }// To symmetrize; we had only inserted the upper triangular part.
    return Hess;
}

/**
 * Some helper functions for inner product and matrix vector multiplication.
 * @param v1
 * @param v2
 * @return
 */
Interval dot(IVector &v1, IVector &v2) {
    Interval result = 0;
    for( int i = 0; i < 4; i++ )
        result += v1[i] * v2[i];
    return result;
}

void Multiply(IMatrix &M, IVector &v, IVector &result) {
    for( int i = 0; i < 4; i++ ) {
        result[i] = 0;
        for( int j = 0; j < 4; j++ )
            result[i] += M[i][j] * v[j];
    }//
}//Matrix-Vector multiplication


/**
 * This function computes the tangential Hessian and the Gauss-Kronecker curve at the point
 * wz=(w_1,w_2,z_1,z_2)
 * @param wz
 * @return
 */
Interval TangentialHessian(IVector &wz) {
    IMatrix Hess = LC_Hessian(wz);
    IVector DH(4);
    LC_grad(wz, DH);
    IMatrix I(4,4), J(4,4), K(4,4); //We obtain the standard quaternions
    I[0][2] = 1;
    I[1][3] = 1;
    I[2][0] = -1;
    I[3][1] = -1;

    J[0][1] = 1;
    J[1][0] = -1;
    J[2][3] = -1;
    J[3][2] = 1;
    K = I * J;

    IVector IDH(4), JDH(4), KDH(4);
    Multiply(I, DH, IDH);
    Multiply(J, DH, JDH);
    Multiply(K, DH, KDH);


    IMatrix THess(3,3);
    IVector tmp(4);
    Multiply(Hess, IDH, tmp);

    THess[0][0] = dot(IDH, tmp);
    THess[1][0] = dot(JDH, tmp);
    THess[2][0] = dot(KDH, tmp);
    Multiply(Hess, JDH, tmp);
    THess[1][1] = dot(JDH, tmp);
    THess[2][1] = dot(KDH, tmp);
    Multiply(Hess, KDH, tmp);
    THess[2][2] = dot(KDH, tmp);

    for( int i=0; i < 3; i++ ) {
        for (int j = 0; j < i; j++)
            THess[j][i] = THess[i][j];
    } //symmetrize
    Interval det = THess[0][0]*(THess[1][1]*THess[2][2]-THess[1][2]*THess[2][1]);
    det -= THess[1][0]*(THess[0][1]*THess[2][2]-THess[0][2]*THess[2][1]);
    det += THess[2][0]*(THess[0][1]*THess[1][2]-THess[0][2]*THess[1][1]);
    return det;
}

/**
 *  This function tests whether the Gauss-Kronecker curvature is positive on a sets parametrized by
 *  the index range.
 * @param argc
 * @param argv
 * @return
 */
int main(int argc, char** argv) {
    if( argc != 7 )
        return 1;
    mu = Interval(std::stod(argv[1]), std::stod(argv[2]));
    c = Interval(std::stod(argv[3]), std::stod(argv[4]));
    double delta = std::stod(argv[5]); //0.005
    double wDelta = std::stod(argv[6]); //0.005;

    IVector xI(4);

    std::cout << "mu: " << mu << std::endl;
    std::cout << "c: " << c << std::endl;
    std::cout << "deltas: " << delta << ", " << wDelta << std::endl;


    int success = 0;
    int failures = 0;
    for( int l=-int(2.0 / wDelta); l< int(2.0 / wDelta) + 1; l++) {
        xI[0] = Interval(l*wDelta, (l+1) * wDelta);
        std::cout << "l=" << l << std::endl;
        for (int k = -int(2.0 / wDelta); k < int(2.0 / wDelta) + 1; k++) {
            xI[1] = Interval(k* wDelta, (k+1) * wDelta);
            for (int j = -int(0.6/delta); j < int(0.6/delta)+1; j++) {
                xI[3] = Interval(j*delta, (j+1) * delta);
                for (int i = 0; i < int(0.6/delta)+1; i++) {
                    xI[2] = Interval(i*delta, (i+1)*delta);
                    // The following is the criterion from Lemma 4.1.
                    if(  5 * 2*(sq(xI[2]) + sq(xI[3]) )> (3 - 2 * mu) )
                        break;
                    Interval LC = LC_Hamiltonian(xI);
                    if (!(LC < 0) && !(LC>0) ) {
                        Interval Hess = TangentialHessian(xI);
                        if ( !(Hess>0) ) {
                            failures++;
                            std::cout << i << ", " << j << ": " << LC << std::endl;
                            std::cout << xI << std::endl;
                            std::cout << 2 * ( sq(xI[2]) - sq(xI[3]) ) - mu << ", "
                                      << 4* xI[2] * xI[3] << std::endl;
                            std::cout << "Hess" << Hess << std::endl;
                        }//if
                        else {
                            success++;
                        }
                    }//Box covers part of hypersurface
                }//z1
            }//z2
        }//w2
    }//w1
    std::cout << "Successes: " << success << std::endl;
    std::cout << "Failed validations: " << failures << std::endl;
    return 0;
}
