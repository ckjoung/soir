#include <iostream>
#include <vector>
#include <queue>
#include "capd/capdlib.h"
#include "shooting.h"
#include "validate.h"
#include "utils.h"
using namespace std;
using namespace capd;

/**
 * Validate that the retrograde orbit has minimal action among all symmetric orbits passing through given fixed point locus.
 * 
 * @return Pair containing validation result and validated location of orbit.
*/
pair<bool, interval> validate_retrograde_action_minimality(interval mu, interval c, int axis, double r)
{
    int sign = sign_retrograde(axis);
    Shooting shooting(mu, c, axis, sign, true);
    double x = shooting.searchSymmetricOrbit();
    interval retro_nbd = interval(x - r, x + r); // neighborhood of retrograde orbit

    // verify local uniqueness
    int local_uniqueness = check_symmetric_existence_C1(mu, c, axis, sign, retro_nbd);
    if (local_uniqueness != 1)
    {
        cout << "Failed to verify local uniqueness" << endl << endl;
        return make_pair(false, retro_nbd);
    }

    // compute action and verify minimal action locally
    cout << "Action computation and local verification of minimality" << endl;
    int coord_action = 4;
    int section = 1 - axis; // shoot to original axis 
    shooting.setPoincareSection(section, Both);

    IVector p = shooting.fixedPointLocus(retro_nbd, sign);
    C0HOTripletonSet s(p);
    cout << "Initial point: " << p << endl;

    IVector q = (*shooting.pm)(s);
    double action_bound = q[coord_action].rightBound();
    cout << "Shooting to original axis arrives at: " << q << endl;
    cout << "Action bound: " << action_bound << endl;

    IVector q2 = (*shooting.pm)(s);
    cout << "Shooting to original axis (again) arrives at: " << q2 << endl;
    if (!(q2[coord_action] > action_bound))
    {
        cout << "Failed to verify action minimality in a neighborhood of retrograde orbit" << endl << endl;
        return make_pair(false, retro_nbd);
    }
    cout << "Local action minimality verified" << endl << endl;
    
    // verify non-existence on complement of neighborhood
    cout << "Start verification of non-existnece of symmetric orbits on the complement" << endl;
    double b = shooting.searchHillBoundary();
    cout << "Boundary of Hill's region: " << b << endl << endl;

    bool result1 = validate_symmetric_nonexistence_below_action(mu, c, axis, sign, interval(0.0, retro_nbd.leftBound()), action_bound);
    bool result2 = validate_symmetric_nonexistence_below_action(mu, c, axis, sign, interval(retro_nbd.rightBound(), b), action_bound);
    bool result3 = validate_symmetric_nonexistence_below_action(mu, c, axis, -sign, interval(0.0, b), action_bound);
    bool success = (result1 && result2 && result3);

    return make_pair(success, retro_nbd);
}

int main(int argc, char *argv[]) 
{
    if (argc != 7)
    {
        cout << "Usage: retrograde_minimal_action <c_left> <c_right> <mu_left> <mu_right>\n"
             << "                                 <axis> <nbd_radius>\n";
        return 0;
    }

    const interval c(argv[1], argv[2]);
    const interval mu(argv[3], argv[4]);
    const int axis = stoi(argv[5]);
    const double r = stod(argv[6]); 

    cout.precision(16);
    cout << "Validating minimal action of retrograde orbit for\n"
         << "c interval: " << c << "\n"
         << "mu interval: " << mu << "\n"
         << "axis: " << axis << "\n";

    auto result = validate_retrograde_action_minimality(mu, c, axis, r);
    bool success = result.first;
    interval retro_nbd = result.second;

    cout << "Results\n"
         << "Retrograde orbit existence validated at: " << retro_nbd << "\n"
         << "Success: " << success << "\n";

    return 0;
}
