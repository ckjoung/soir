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
 * Validate that the retrograde and direct orbits have minimal action among all symmetric orbits passing through 
 * given fixed point locus.
 * 
 * @return Pair containing validation result and a pair of validated location of the two orbits.
*/
pair<bool, pair<Interval, Interval>> validate_symmetric_action_minimality(interval mu, interval c, int axis, double r)
{
    int sign_retro = sign_retrograde(axis);
    int sign_direct = -sign_retro;

    // verify direct orbit
    Shooting shooting(mu, c, axis, sign_direct, true);
    double x_direct = shooting.searchSymmetricOrbit();
    interval nbd_direct = interval(x_direct - r, x_direct + r); // neighborhood of direct orbit

    int local_uniqueness_direct = check_symmetric_existence_C1(mu, c, axis, sign_direct, nbd_direct);
    if (local_uniqueness_direct != 1)
    {
        cout << "Failed to verify direct orbit local uniqueness" << endl << endl;
        return make_pair(false, make_pair(interval(0.0), nbd_direct));
    }

    // compute action bound
    cout << "Action computation and local verification of minimality" << endl;
    int coord_action = 4;
    int section = 1 - axis; // shoot to original axis 
    shooting.setPoincareSection(section, Both);

    IVector p = shooting.fixedPointLocus(nbd_direct, sign_direct);
    C0HOTripletonSet s(p);
    cout << "Initial point: " << p << endl;

    IVector q = (*shooting.pm)(s);
    double action_bound = q[coord_action].rightBound();
    cout << "Shooting to original axis arrives at: " << q << endl;
    cout << "Action bound: " << action_bound << endl;

    // verify for neighborhood of direct orbit
    IVector q2 = (*shooting.pm)(s);
    cout << "Shooting to original axis (again) arrives at: " << q2 << endl;
    if (!(q2[coord_action] > action_bound))
    {
        cout << "Failed to verify action minimality in a neighborhood of direct orbit" << endl << endl;
        return make_pair(false, make_pair(interval(0.0), nbd_direct));
    }
    cout << "Local action minimality verified" << endl << endl;

    // verify retrograde orbit
    shooting.setFixedPointLocus(axis, sign_retro);
    double x_retro = shooting.searchSymmetricOrbit();
    interval nbd_retro = interval(x_retro - r, x_retro + r); // neighborhood of retrograde orbit

    int local_uniqueness_retro = check_symmetric_existence_C1(mu, c, axis, sign_retro, nbd_retro);
    if (local_uniqueness_retro != 1)
    {
        cout << "Failed to verify retrograde orbit local uniqueness" << endl << endl;
        return make_pair(false, make_pair(nbd_retro, nbd_direct));
    }

    // verify for neighborhood of retrograde orbit
    p = shooting.fixedPointLocus(nbd_retro, sign_retro);
    C0HOTripletonSet s1(p);
    q2 = (*shooting.pm)(s1, 2);
    if (!(q2[coord_action] > action_bound))
    {
        cout << "Failed to verify action minimality in a neighborhood of retrograde orbit" << endl << endl;
        return make_pair(false, make_pair(nbd_retro, nbd_direct));
    }

    // verify non-existence on complement of neighborhood
    double b = shooting.searchHillBoundary();
    bool result1 = validate_symmetric_nonexistence_below_action(mu, c, axis, sign_retro, interval(0.0, nbd_retro.leftBound()), action_bound);
    bool result2 = validate_symmetric_nonexistence_below_action(mu, c, axis, sign_retro, interval(nbd_retro.rightBound(), b), action_bound);
    bool result3 = validate_symmetric_nonexistence_below_action(mu, c, axis, sign_direct, interval(0.0, nbd_direct.leftBound()), action_bound);
    bool result4 = validate_symmetric_nonexistence_below_action(mu, c, axis, sign_direct, interval(nbd_direct.rightBound(), b), action_bound);
    bool success = (result1 && result2 && result3 && result4);

    return make_pair(success, make_pair(nbd_retro, nbd_direct));
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
    cout << "Validating minimal action of direct orbit for\n"
         << "c interval: " << c << "\n"
         << "mu interval: " << mu << "\n"
         << "axis: " << axis << "\n";

    auto result = validate_symmetric_action_minimality(mu, c, axis, r);
    bool success = result.first;
    interval retro_nbd = result.second.first;
    interval direct_nbd = result.second.second;

    cout << "Results\n"
         << "Retrograde orbit existence validated at: " << retro_nbd << "\n"
         << "Direct orbit existence validated at: " << direct_nbd << "\n"
         << "Success: " << success << "\n";
         
    return 0;
}
