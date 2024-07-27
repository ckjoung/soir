#ifndef VALIDATE_H
#define VALIDATE_H
#include <iostream>
#include <vector>
#include "capd/capdlib.h"
#include "shooting.h"
#include "utils.h"
using namespace std;
using namespace capd;

/**
 * Checks whether there exists a symmetric orbit in the given z parameter interval of the fixed point locus, by C0 shooting. 
 * 
 * @return 1 if existence is verified, 0 otherwise.
*/
int check_symmetric_existence_C0(interval mu, interval c, int axis, int sign, interval z)
{
    Shooting shooting(mu, c, axis, sign);
    interval zl = z.leftBound(), zr = z.rightBound(), wl, wr;
    IVector p = shooting.fixedPointLocus(z, sign);
    IVector pl = shooting.fixedPointLocus(zl, sign);
    IVector pr = shooting.fixedPointLocus(zr, sign);
    shooting.setPoincareSection(axis, Both); // shoot to perpendicular axis
    int coordw = slope_coordinate(axis); 

    try // C0 shooting at endpoints
    {
        wl = shooting.IC0Trace(pl)[coordw]; 
        wr = shooting.IC0Trace(pr)[coordw]; 
    }
    catch (std::exception& e) 
    {
        return 0; // C0 shooting at endpoints fail
    }

    if (wl * wr < 0) // sign of arriving slope different at endpoints
        return 1;
    return 0;
}

/**
 * Checks whether there is no symmetric orbit in the given z parameter interval of the fixed point locus, by C0 shooting. 
 * 
 * @return -1 if non-existence is verified, 0 otherwise.
*/
int check_symmetric_nonexistence_C0(interval mu, interval c, int axis, int sign, interval z)
{
    Shooting shooting(mu, c, axis, sign);
    IVector p = shooting.fixedPointLocus(z, sign);
    IVector q;

    int section = axis; // shoot to perpendicular axis
    shooting.setPoincareSection(section, Both); 
    int slope_coord = slope_coordinate(section); 
    try
    {
        q = shooting.IC0Trace(p); 
    }
    catch (std::exception& e) 
    {
        section = 1 - axis; // shoot to original axis
        shooting.setPoincareSection(section, Both); 
        slope_coord = slope_coordinate(section); 
        try
        {
            q = shooting.IC0Trace(p); 
        }
        catch (std::exception& e) 
        {
            return 0; // fail to shoot to both axes
        }
    }
    
    if (!q[slope_coord].contains(0.0)) 
        return -1; // non-existence verified
    else 
        return 0; 
}

/**
 * Checks whether there exists a symmetric orbit in the given z parameter interval of the fixed point locus, by C1 shooting.
 *  
 * @return -1 if non-existence is verified, 1 if local uniqueness verified, 0 otherwise.
*/
int check_symmetric_existence_C1(interval mu, interval c, int axis, int sign, interval z)
{
    Shooting shooting(mu, c, axis, sign);
    interval zl = z.leftBound(), zr = z.rightBound(), w, wl, wr;
    IVector p = shooting.fixedPointLocus(z, sign);
    IVector pl = shooting.fixedPointLocus(zl, sign);
    IVector pr = shooting.fixedPointLocus(zr, sign);
    IVector q;
    IMatrix Dphi, DP, DF;
    int section = axis; // shoot to perpendicular axis
    shooting.setPoincareSection(section, Both); 
    int slope_coord = slope_coordinate(section); 

    try 
    {
        wl = shooting.IC0Trace(pl)[slope_coord]; 
        wr = shooting.IC0Trace(pr)[slope_coord]; 
        shooting.IC1Trace(p, q, Dphi, DP);
        DF = shooting.fixedPointLocusDerivative(z, sign);
    }
    catch(const std::exception& e)
    {
        return 0; // shooting failed
    }

    interval dw = (DP * DF)[slope_coord][0];
    intersection(wl + (z - zl) * dw, wr + (z - zr) * dw, w);

    if (!w.contains(0.0))
        return -1; // non-existence verified
    else if (wl * wr < 0 && !dw.contains(0.0))
        return 1; // existence and local uniqueness verified
    else 
        return 0;
}

/**
 * Checks existence of a symmetric orbit within the given interval.
 * 
 * @return -1 if non-existence validated, 1 if existence and local uniqueness validated, 0 otherwise
*/
int check_symmetric_existence(interval mu, interval c, int axis, int sign, interval z)
{
    cout << "Try C0 validation..." << endl;
    int result0 = check_symmetric_nonexistence_C0(mu, c, axis, sign, z);
    if (result0 == -1) // non-existence validated
        return result0;

    cout << "Result indecisive. Try C1 validation..." << endl;
    int result1 = check_symmetric_existence_C1(mu, c, axis, sign, z);
    return result1;
}

/**
 * This function checks existence of symmetric orbit within a given interval below given action bound.
 * 
 * @return -1 if non-existence validated, 0 otherwise
*/
int check_symmetric_existence_below_action(interval mu, interval c, int axis, int sign, interval z, double action_bound)
{
    Shooting shooting(mu, c, axis, sign, true);
    IVector p = shooting.fixedPointLocus(z, sign);
    cout << "Initial point: " << p << endl;
    IVector q;
    int coord_action = 4;

    // initial shooting
    cout << "Try shooting to section (0: z1=0, 1: z2=0): " << 1 - axis << endl;
    try
    {
        shooting.setPoincareSection(1 - axis, Both); // shoot to original axis
        cout << "Shooting from: " << p << endl;
        q = shooting.IC0Trace(p);
        cout << "Arrived at: " << q << endl;
    }
    catch(std::exception& e)
    {
        cout << "Failed, instead shooting to (0: z1=0, 1: z2=0): " << axis << endl;
        try
        {
            shooting.setPoincareSection(axis, Both); // shoot to perpendicular axis
            cout << "Shooting from: " << p << endl;
            q = shooting.IC0Trace(p);
            cout << "Arrived at: " << q << endl;
        }
        catch(std::exception& e)
        {
            cout << "Both shooting failed." << endl;
            return 0; // both shooting fails
        }
    }

    C0HOTripletonSet s(p);
    q = (*shooting.pm)(s); // shoot to available section
    cout << "Current action: " << q[coord_action] << endl;
    cout << "Exceeds action bound (" << action_bound << ")?: " << (q[coord_action] > action_bound) << endl;

    // repeated shooting until action bound is reached
    shooting.setPoincareSection(1 - axis, Both); // shoot to original axis

    while (!(q[coord_action] > action_bound))
    {
        if ((q[0].contains(0.0) && q[3].contains(0.0)) || (q[1].contains(0.0) && q[2].contains(0.0)))
        {
            cout << "Arriving vector intersects fixed point locus, possible existence of symmetric orbit." << endl;
            return 0; // possible existence of symmetric orbit
        }

        // shoot to next intersection
        cout << "Try shooting to next intersection with section (0: z1=0, 1: z2=0): " << 1 - axis << endl;
        try 
        {
            cout << "Shooting from: " << q << endl;
            q = (*shooting.pm)(s);
            cout << "Arrived at: " << q << endl;
        }
        catch(std::exception& e)
        {
            cout << "Shooting failed." << endl;
            return 0; // shooting fails
        }

        cout << "Current action: " << q[coord_action] << endl;
        cout << "Exceeds action bound (" << action_bound << ")?: " << (q[coord_action] > action_bound) << endl;
    }
    return -1; // non-existence verified
}

/**
 * This function validates non-existence of symmetric orbit within a given interval below given action bound using a subdivision algorithm.
 * 
 * @return True if non-existence verified, false otherwise.
*/
bool validate_symmetric_nonexistence_below_action(interval mu, interval c, int axis, int sign, interval z, double action_bound, double delta=1e-3)
{
    queue<interval> Q;
    initialize_queue(Q, z.leftBound(), z.rightBound(), delta);
    int n_subdivision = 0, max_subdivisions = 1e3;
    while (!Q.empty() && n_subdivision < max_subdivisions)
    {
        cout << "Current queue size: " << Q.size() << endl;
        cout << "No. of subdivisions: " << n_subdivision << "/" << max_subdivisions << endl;
        interval z = Q.front();
        Q.pop();

        cout << "Try non-existence verification for interval: " << z << endl;
        int result = check_symmetric_existence_below_action(mu, c, axis, sign, z, action_bound);
        if (result == -1)
        {
            cout << "Non-existence verified" << endl << endl;
            continue;
        }

        cout << "Fail to verify non-existence, subdividing interval" << endl << endl;
        subdivide_and_push(Q, z);
        n_subdivision++;
    }
    bool result = Q.empty();
    return result;
}

/**
 * Validates retrograde or direct property of a periodic orbit, depending on given rotation sign.
 *  
 * @return True if the property is validated, False otherwise.
*/
bool validate_angle_monotonicity(interval mu, interval c, int axis, int sign, interval z, int rotation_sign)
{
    Shooting shooting(mu, c, axis, sign);
    shooting.setPoincareSection(axis, Both); // shoot to perpendicular axis
    IOdeSolver *solver = shooting.solver;
    ITimeMap *tm = shooting.tm;
    IVector p = shooting.fixedPointLocus(z, sign);
    C0HOTripletonSet s0(p);
    interval return_time;
    (*shooting.pm)(s0, return_time);
    cout << "Return time: " << return_time << endl;
    C0HOTripletonSet s(p);

    tm->stopAfterStep(true);
    do
    {
        (*tm)(return_time.rightBound(), s);
        IOdeSolver::SolutionCurve curve = solver->getCurve();
        interval step = solver->getStep();
        interval time_domain = interval(0, 1) * step;
        IVector q = curve(time_domain);
        interval phidot = shooting.angleDerivative(q);
        cout << "phidot: " << phidot << endl;
        if (!(phidot * rotation_sign > 0))
            return false;
    }
    while (!(*tm).completed());
    return true;
}

/**
 * Validates non-degeneracy of periodic orbit using Lemma 2.7 (trace condition)
 * @param k number of iterations of Poincare map for section z1=0 or z2=0
 
 * @return True if validated, False otherwise.
*/
bool validate_nondegeneracy(interval mu, interval c, int axis, int sign, interval z, int k=1) 
{
    Shooting shooting(mu, c, axis, sign);
    shooting.setPoincareSection(1 - axis, Both); // shoot to original axis
    IVector p = shooting.fixedPointLocus(z, sign);
    IVector q;
    IMatrix Dphi, DP;
    shooting.IC1Trace(p, q, Dphi, DP, k);

    IMatrix S(2, 2);
    shooting.symplecticLinearization(p, q, Dphi, S);
    interval trace = S[0][0] + S[1][1];
    cout << "S: " << S << endl;
    cout << "trace: " << trace << endl;
    if (trace.contains(2.0)) 
    {
        cout << "validation failed" << endl << endl;
        return false;
    }
    cout << "validation success" << endl << endl;
    return true;
}

#endif // VALIDATE_H