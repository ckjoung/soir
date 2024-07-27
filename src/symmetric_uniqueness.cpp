/**
This code validates existence of exactly two symmetric orbits with crossing number 2 for given mu and c intervals on the fixed point locus of the specified involution.
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
 * This function validates non-existence of a symmetric orbit within a given interval using a subdivision algorithm.
 * 
 * @return True if non-existence verified, false otherwise.
*/
bool validate_symmetric_nonexistence(interval mu, interval c, int axis, int sign, interval z)
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
 * Validates uniqueness of symmetric orbit with crossing number two for given axis and sign.
 * 
 * @return Pair containing validation result and validated location of orbit.
*/
pair<bool, interval> validate_symmetric_uniqueness(interval mu, interval c, int axis, int sign, double r = 1e-3)
{
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
    bool result1 = validate_symmetric_nonexistence(mu, c, axis, sign, interval(0.0, nbd.leftBound()));
    bool result2 = validate_symmetric_nonexistence(mu, c, axis, sign, interval(nbd.rightBound(), b));
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
    cout << "Validating uniqueness of symmetric orbits for\n"
         << "c interval: " << c << "\n"
         << "mu interval: " << mu << "\n"
         << "axis: " << axis << "\n";

    int sign_retro = sign_retrograde(axis);
    int sign_direct = -sign_retro;
    auto result1 = validate_symmetric_uniqueness(mu, c, axis, sign_retro, r);
    auto result2 = validate_symmetric_uniqueness(mu, c, axis, sign_direct, r);
    bool success = result1.first && result2.first;

    cout << "Results\n"
         << "Retrograde orbit existence validated at: " << result1.second << "\n"
         << "Direct orbit existence validated at: " << result2.second << "\n"
         << "Retrograde orbit uniqueness: " << result1.first << "\n"
         << "Direct orbit uniqueness: " << result2.first << "\n"
         << "Success: " << success << "\n";
    
    return 0;
}
