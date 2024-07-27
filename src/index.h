#ifndef INDEX_H
#define INDEX_H
#include <iostream>
#include <vector>
#include "shooting.h"
#include "capd/capdlib.h"
using namespace capd;
using namespace std;

/**
 * Compute retract rho: Sp(2) -> U(1) using formula (2.7).
 * 
 * 
*/
IMatrix retract(IMatrix A) 
{
    interval a = A[0][0];
    interval b = A[0][1];
    interval c = A[1][0];
    interval d = A[1][1];
    
    interval t1 = a + d;
    interval t2 = b - c;
    interval t3 = sqrt(t1 * t1 + t2 * t2);
    interval t4 = t1 / t3;
    interval t5 = t2 / t3;

    IMatrix U(2, 2);
    U[0][0] = t4;
    U[0][1] = t5;
    U[1][0] = -t5;
    U[1][1] = t4;
    return U;
}

/**
 * Get index of the cover which contains the given unitary (interval) matrix.
 *       U_1 = {x > 0}
 *       U_2 = {y > 0}
 *       U_3 = {x < 0}
 *       U_4 = {y < 0}
 * 
 *  
*/
int coveringIndex(IMatrix U) 
{
    interval X = U[0][0]; // Re(U)
    interval Y = U[1][0]; // Im(U)

    if (Y > 0)
        return 2;
    if (Y < 0)
        return 4;
    if (X < 0)
        return 3;
    if (X > 0)
        return 1;

    cout << "Proper covering not found" << endl;
    cout << "given complex number: " << X << " + " << Y << "i" << endl;
    throw runtime_error("Runtime error: proper covering not found.");
}


/**
 * Calculate winding number (up to error of 1/4) from a vector of covering indices.
 *   
 */
double windingNumber(std::vector<int> coveringIndices) 
{
    double wn = 0.0; // winding number
    int prevIdx = 1;

    for (int currIdx: coveringIndices) 
    {
        if ((currIdx - prevIdx + 4) % 4 == 1) 
            wn += 0.25; // positive quarter rotation (1->2->3->4->1)
        else if ((currIdx - prevIdx + 4) % 4 == 3)
            wn -= 0.25; // negative quarter rotation (1->4->3->2->1)

        prevIdx = currIdx;
    }
    return wn;
}


/**
 * Helper functions for evaluation of the Conley-Zehnder index using Proposition 2.12
 */
int parityType(IMatrix Dphi) 
{
    interval trace = Dphi[0][0] + Dphi[1][1];

    if (trace < 2) 
        return 1;
    if (trace > 2) 
        return 0;

    throw runtime_error("Parity type not determined.");
}

int parityCorrectedIndex(double wn, int parity) 
{
    int n = std::floor(wn);
    if (std::max(n, -n) % 2 == parity)
        return n;
    return n + 1;
}

/**
 * Get enclosure of linearized flow with respect the U, V frame from Formula (2.5)
 * Time steps are discretized in n_grid=128 steps; this is fine enough for the short orbits
 * considered in this paper.
 */
int computeCZIndex(interval mu, interval c, int axis, int sign, interval z, int k) 
{
    Shooting shooting(mu, c, axis, sign);
    shooting.setPoincareSection(1 - axis, Both); // shoot to original axis
    IOdeSolver *solver = shooting.solver;
    ITimeMap *tm = shooting.tm;
    IVector p = shooting.fixedPointLocus(z, sign);
    C0HORect2Set s0(p);
    interval return_time;
    interval prev_time, current_time;
    (*shooting.pm)(s0, return_time, k);
    cout << "Return time: " << return_time << endl;

    C1HORect2Set s1(p);
    std::vector<int> covering_indices;
    int parity;
    IMatrix S(2, 2);

    tm->stopAfterStep(true);
    do
    {
        (*tm)(return_time.rightBound(), s1);
        IOdeSolver::SolutionCurve curve = solver->getCurve();
        interval step = solver->getStep();
        interval time_domain = interval(0, 1) * step;

        int n_grid = 128;
        for(int i = 0; i < n_grid; i++)
        {
            interval t = interval(i, i+1) * step / n_grid;
            intersection(time_domain, t, t);
            current_time = prev_time + t;

            IVector q = curve(t); // enclosure of trajectory
            IMatrix Dphi = curve[t]; // enclosure of monodromy matrix
            shooting.symplecticLinearization(p, q, Dphi, S); // TODO: assert S is symplectic matrix
            IMatrix U = retract(S); // TODO: assert U is unitary
            int covering_index = coveringIndex(U*U);

            if (current_time.contains(return_time.leftBound()))
            {
                parity = parityType(S);
            }
            if (current_time.leftBound() > return_time.leftBound())
            {
                if (parity != parityType(S))
                    throw runtime_error("Index computation failed due to parity change during time enclosure of return time.");

                if (covering_index != covering_indices.back())
                    throw runtime_error("Index computation failed due to covering change during time enclosure of return time.");
            }

            covering_indices.push_back(covering_index);
            cout << "current time: " << current_time <<  endl;
            cout << "S: " << S <<  endl;
            cout << "U^2: " << (U*U)[0][0] << " + " << (U*U)[1][0] << "i" <<  endl;
            cout << "covering index: " << covering_index << endl << endl;
        }
        prev_time = tm->getCurrentTime();
    }
    while (!(*tm).completed());
    double wn = windingNumber(covering_indices);
    int index = parityCorrectedIndex(wn, parity);
    cout << "Winding number: " << wn << endl;
    cout << "Parity: " << parity << endl;
    cout << "Conley-Zehnder index: " << index << endl << endl;
    return index;
}
#endif // INDEX_H