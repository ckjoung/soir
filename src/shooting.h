#ifndef SHOOTING_H
#define SHOOTING_H
#include <iostream>
#include "capd/capdlib.h"
#include "rtbp.h"
#include "utils.h"
using namespace capd;
using namespace capd::poincare;
using namespace std;

class Shooting 
{
public:
    interval mu;
    interval c;
    int dim; // dimension of phase space
    int axis; // axis to perform shooting on (0: z1, 1: z2)
    int sign; // sign(+1 or -1) of w for fixed point locus parametrization
    IMap *vf;
    IOdeSolver *solver;
    ITimeMap *tm;
    int sectionCoord;    
    ICoordinateSection *section;
    CrossingDirection crossing;
    IPoincareMap *pm;

    Shooting(interval mu, interval c, int axis, int sign, bool includeAction=false) 
    {
        this->mu = mu;
        this->c = c;
        initSolver(includeAction); 
        setFixedPointLocus(axis, sign);
    }

    void initSolver(bool includeAction=false) 
    {
        vf = new IMap(LCRTBPStr(includeAction));
        vf->setParameter("mu", mu);
        vf->setParameter("c", c);
        dim = vf->dimension();
        solver = new IOdeSolver(*vf, 20);
        tm = new ITimeMap(*solver);
    }

    void setFixedPointLocus(int axis, int sign)
    {
        this->axis = axis;
        this->sign = sign;
    }

    void setPoincareSection(int sectionCoord, CrossingDirection crossing) 
    {
        this->sectionCoord = sectionCoord;
        this->crossing = crossing;
        section = new ICoordinateSection(dim, sectionCoord);
        pm = new IPoincareMap(*solver, *section, crossing);
    }

    IVector IC0Trace(IVector &p, int k=1) 
    {
        C0HOTripletonSet s(p);
        IVector q = (*pm)(s, k);
        return q;
    }

    void IC1Trace(IVector &p, IVector &q, IMatrix &Dphi, IMatrix &DP, int k=1) 
    {
        C1HORect2Set s(p);
        Dphi(dim, dim);
        q = (*pm)(s, Dphi, k);
        DP = (*pm).computeDP(q,Dphi);
    }

    /**
     * The functions w_quadratic_coeffs, discriminant and fixedPointLocus are helper functions
     * for parametrizing the fixed point locus. 
     * They correspond to 
     * Formulas for $w^\pm_2$ and $w^\pm_1$ in Section (2.4.2) on
     * Parametrizing the fixed point locus.
     *  
     */
    pair<interval, interval> w_quadratic_coeffs(interval z)
    {
        interval B, C; 
        if (axis == 0) // z = z1
        {
            B = 2.0 * z * z * z - mu * z;
            C = - 0.5 + c * z * z + mu / 2.0 - mu * z * z / (1 - 2.0 * z * z);
        }
        else // axis = 1 and z = z2
        {
            B = - 2.0 * z * z * z - mu * z;
            C = - 0.5 + c * z * z + mu / 2.0 - mu * z * z / (1 + 2.0 * z * z);
        }
        return make_pair(B, C);
    }

    interval discriminant(interval z) 
    {
        auto coeffs = w_quadratic_coeffs(z);
        interval B = coeffs.first; 
        interval C = coeffs.second; 
        return B * B - 2 * C; 
    }

    /**
     * Lift given z value to corresponding point in the fixed point locus for given energy.
     * 
     * @param z interval for z coordinate
    */
    IVector fixedPointLocus(interval z, int sign) 
    {
        auto coeffs = w_quadratic_coeffs(z);
        interval B = coeffs.first; 
        interval C = coeffs.second; 
        interval D = B * B - 2 * C; // discriminant

        if (D.leftBound() < 0.0)
            D = interval(0.0, std::max(0.0, D.rightBound())); 
            // intersect D with [0, inf) to get real roots only

        interval w = (- B + sign * sqrt(D)); // solve for w

        IVector p(dim); // initialized to be vector of zeros
        p[axis] = z;
        p[w_coordinate(axis)] = w;
        return p;
    }

    /**
     * Compute derivative of fixed point locus parametrization.
     * 
     * @param z interval for z coordinate
     * @return DF
    */
    IMatrix fixedPointLocusDerivative(interval z, int sign) 
    {
        auto coeffs = w_quadratic_coeffs(z);
        interval B = coeffs.first; 
        interval C = coeffs.second; 
        interval dB, dC; // coefficients and their derivatives

        if (axis == 0) // z = z1
        {
            dB = 6.0 * z * z - mu;
            dC = 2.0 * c * z - 2.0 * mu * z / sqr(1 - 2.0 * z * z);
        }
        else // axis = 1 and z = z2
        {
            dB = - 6.0 * z * z - mu;
            dC = 2.0 * c * z - 2.0 * mu * z / sqr(1 + 2.0 * z * z);
        }

        interval D = B * B - 2 * C; // discriminant
        if (D.leftBound() <= 0.0)
            throw runtime_error("Runtime error: non-positive discriminant encountered while computing derivative of fixed point locus parametrization.");

        interval dw = - dB + sign * (B * dB - dC) / sqrt(D); // solve for w

        IMatrix dF(dim, 1); // initialized to be matrix of zeros
        dF[axis][0] = interval(1.0);
        dF[w_coordinate(axis)][0] = dw;
        return dF;
    }

    /**
     * Simple search for locating the boundary; the large stepsize 1e-4 works for the 
     * mass and energy parameters considered in the paper.
     */
    double searchHillBoundary() 
    {
        // line search
        double b = 0.0;
        interval D;
        double step = 1e-4;
        do 
        {
            b += step;
            D = discriminant(b);
        } while (!(D < 0));
        
        // bisection to refine boundary
        double a = b - step;
        int numIter = 10;
        for (int i = 0; i < numIter; i++) 
        {
            double mid = (a + b) / 2;
            D = discriminant(mid);
            if (D < 0)
                b = mid;
            else
                a = mid;
        }
        return b;
    }

    /**
     * Performs non-rigorous search for symmetric orbit using Newton iterations with given initial z value.
     * 
     * @param searchInterval The interval in which the search is performed
    */
    double searchSymmetricOrbit(double z)
    {
        double delta = 1e-6;
        int numIter = 10;
        int sectionCoord = axis; // perpendicular axis
        int slopeCoord = slope_coordinate(sectionCoord); // coordinate of either w1 or w2
        setPoincareSection(sectionCoord, Both); // shoot to perpendicular axis

        // Newton iterations
        for (int i = 0; i < numIter; i++) 
        {
            IVector p0 = fixedPointLocus(z, sign);
            IVector p1 = fixedPointLocus(z + delta, sign);

            IVector q0 = IC0Trace(p0);
            IVector q1 = IC0Trace(p1);
            
            double y0 = q0[slopeCoord].mid().leftBound();
            double y1 = q1[slopeCoord].mid().leftBound();
            double derivative = (y1 - y0) / delta + 1e-16;  // linearized flow not needed in this non-rigorous search
            
            z -= y0 / derivative;   
        }
        return z;
    }

    /**
     * Performs non-rigorous search for symmetric orbit on the specified fixed point locus using a line search 
     * to find a suitable initial point for Newton iterations. 
    */
    double searchSymmetricOrbit()
    {
        // line search on entire fixed point locus
        double z = 0.0, b = searchHillBoundary();
        double step = b / 100;
        double prevSlope, currSlope = 0.0;
        int sectionCoord = axis; // perpendicular axis
        setPoincareSection(sectionCoord, Both); // shoot to perpendicular axis

        while (z + step < b)
        {
            z += step;
            prevSlope = currSlope;

            IVector p = fixedPointLocus(z, sign);
            try
            {
                IVector q = IC0Trace(p);
                currSlope = q[3 - sectionCoord].mid().leftBound();
            }
            catch(const std::exception& e)
            {
                continue;
            }
            
            if (!(prevSlope * currSlope >= 0))
            {
                try
                {
                    return searchSymmetricOrbit(z - step / 2);
                }
                catch(std::exception& e)
                {
                    continue;
                }
            }
        } 
        throw runtime_error("Runtime error: symmetric orbit search failed.");
    }

    /**
     * This implements the condition for checking retrograde or prograde orbits following 
     * Proposition 2.1 and the corresponding formula (2.1).
     */

    interval angleDerivative(IVector p)
    {
        interval z_norm_sq = p[0] * p[0] + p[1] * p[1];
        interval q1_plus_mu = 2 * (p[0] * p[0] - p[1] * p[1]);
        interval q2 = 4 * p[0] * p[1];
        interval p1 = (p[0] * p[2] - p[1] * p[3]) / z_norm_sq;
        interval p2 = (p[0] * p[3] + p[1] * p[2]) / z_norm_sq;
        return q1_plus_mu * p2 - q2 * p1 - mu * q1_plus_mu;
    }

    /**
     * Compute symplectic form product of two given tangent vectors. 
     * We have ordered the coordinates as (z_1,z_2;w_1,w_2) with symplectic form
     * \omega= dw \wedge dz
    */
    interval omega(IVector &U, IVector &V) 
    {
        return U[2]*V[0] + U[3]*V[1] - U[0]*V[2] - U[1]*V[3]; // omega = dw wedge dz
    }

    /**
     * Compute frame vectors U, V of the symplectic trivialization of formula (2.5)
     * Note X_K = i DK, U = j DK, V = k DK, and ij=k and ki=j.
     * Hence
     * U = j (-i X_K)= k X_K  and V = k(-i) DK = -j X_K,
     * which we implement below.
     * In the paper, the standard quaternions are denoted by I, J, K.  
    */
    void frameVectors(IVector &X, IVector &U, IVector &V) 
    {
        IVector X_K = (*vf)(X);
        interval norm = X_K.euclNorm();

        U[0] = -X_K[3];
        U[1] = X_K[2];
        U[2] = -X_K[1];
        U[3] = X_K[0];

        V[0] = -X_K[1];
        V[1] = X_K[0];
        V[2] = X_K[3];
        V[3] = -X_K[2];

        U = U / norm;
        V = V / norm;
    }

    /**
     * Calculate 2 by 2 symplectic matrix corresponding to the linearized flow on the symplectic 
     * trivialization. This corresponds to formula (2.6).
    */
    void symplecticLinearization(IVector &p, IVector &q, IMatrix &Dphi, IMatrix &S) 
    {
        IVector Ui(4), Vi(4), Uf(4), Vf(4);
        frameVectors(p, Ui, Vi);
        frameVectors(q, Uf, Vf);

        IVector DphiU = Dphi * Ui;
        IVector DphiV = Dphi * Vi;
        S[0][0] = omega(DphiU, Vf);
        S[0][1] = omega(DphiV, Vf);
        S[1][0] = -omega(DphiU, Uf);
        S[1][1] = -omega(DphiV, Uf);
    }
};

#endif // SHOOTING_H