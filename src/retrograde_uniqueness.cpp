/**
This code validates uniqueness of retrograde symmetric orbit with crossing number 2 for given mu and c intervals on the fixed point locus of the specified involution.
*/

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
 * This function checks whether the given interval corresponds to part of a retrograde orbit.
 * 
 * @return -1 if it is validated to be not retrograde, 0 otherwise
*/
int check_initial_retrograde(interval mu, interval c, int axis, int sign, interval z)
{
    double zr = z.rightBound();
    Shooting shooting(mu, c, axis, sign);
    IVector p = shooting.fixedPointLocus(zr, sign);
    interval w = p[w_coordinate(axis)];
    interval phidot = sign * (w - zr * mu); // w_2 - mu * z_1 or mu * z_2 - w_2
    if (phidot < 0)
        return -1; // non-retrograde initial condition validated
    return 0;
}

/**
 * This function validates non-existence of retrograde orbit within a given interval using a subdivision algorithm.
 * 
 * @return True if non-existence verified, false otherwise.
*/
bool validate_retrograde_nonexistence(interval mu, interval c, int axis, int sign, interval z)
{
    queue<interval> Q;
    double delta = 1e-3;
    int n_subdivision = 0, max_subdivisions = 1e3;
    initialize_queue(Q, z.leftBound(), z.rightBound(), delta);
    while (!Q.empty() && n_subdivision < max_subdivisions)
    {
        interval z = Q.front();
        Q.pop();
        cout << "Try non-existence verification for interval: " << z << endl;

        int initial_retro = check_initial_retrograde(mu, c, axis, sign, z);
        if (initial_retro == -1) 
        {
            cout << "Initial condition not retrograde" << endl;
            continue;
        }

        int result = check_symmetric_existence(mu, c, axis, sign, z);
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
 * Validates uniqueness of retrograde orbit up to symmetric orbits with crossing number two on given axis.
 * 
 * @return Pair containing validation result and validated location of orbit.
*/
pair<bool, interval> validate_retrograde_uniqueness(interval mu, interval c, int axis, double r=1e-3)
{
    int sign = sign_retrograde(axis);
    Shooting shooting(mu, c, axis, sign);
    double x = shooting.searchSymmetricOrbit();
    interval nbd = interval(x - r, x + r); // neighborhood of retrograde orbit
    cout << "Non-rigorous symmetric orbit search result: " << x << endl;

    // verify local uniqueness
    int local_uniqueness = check_symmetric_existence_C1(mu, c, axis, sign, nbd);
    if (local_uniqueness != 1)
    {
        cout << "Failed to verify local uniqueness" << endl << endl;
        return make_pair(false, nbd);
    }

    // verify non-existence on complement of neighborhood
    double b = shooting.searchHillBoundary();
    bool result1 = validate_retrograde_nonexistence(mu, c, axis, sign, interval(0.0, nbd.leftBound()));
    bool result2 = validate_retrograde_nonexistence(mu, c, axis, sign, interval(nbd.rightBound(), b));
    bool success = result1 && result2;
    return make_pair(success, nbd);
}

int main(int argc, char *argv[]) 
{
    if (argc != 7)
    {
        cout << "Usage: retrograde_uniqueness <c_left> <c_right> <mu_left> <mu_right>\n"
             << "                             <axis> <nbd_radius>\n";
        return 0;
    }

    const interval c(argv[1], argv[2]);
    const interval mu(argv[3], argv[4]);
    const int axis = stoi(argv[5]);
    const double r = stod(argv[6]); 

    cout.precision(16);
    cout << "Validating uniqueness of retrograde orbit for\n"
         << "c interval: " << c << "\n"
         << "mu interval: " << mu << "\n"
         << "axis: " << axis << "\n";

    auto result = validate_retrograde_uniqueness(mu, c, axis, r);
    bool success = result.first;
    interval retro_nbd = result.second;

    cout << "Results\n"
         << "Retrograde orbit existence validated at: " << retro_nbd << "\n"
         << "Uniqueness: " << success << "\n";
         
    return 0;
}
