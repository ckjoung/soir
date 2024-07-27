#include <iostream>
#include "capd/capdlib.h"
#include "shooting.h"
#include "validate.h"
#include "utils.h"
#include "index.h"
using namespace std;
using namespace capd;

int main(int argc, char *argv[]) 
{   
    if (argc != 11)
    {
        cout << "Usage: orbit_properties <c_left> <c_right> <mu_left> <mu_right>\n"
             << "                        <w_sign> <z1> <nbd_radius> <n_poincare_iterates>\n"
             << "                        <n_subdivision_angle> <n_subdivision_nondegeneracy>\n";
        return 0;
    }

    const interval c(argv[1], argv[2]);
    const interval mu(argv[3], argv[4]);
    const int sign = stoi(argv[5]);
    double z1 = stod(argv[6]);
    const double r = stod(argv[7]);
    const int k = stoi(argv[8]);
    const int Na = stoi(argv[9]);
    const int Nn = stoi(argv[10]);

    const int axis = 0; // z1-axis
    Shooting shooting(mu, c, axis, sign);

    if (z1 == 0.0) // perform line search if estimate not given
        z1 = shooting.searchSymmetricOrbit();
    interval nbd = interval(z1 - r, z1 + r);

    cout.precision(16);
    cout << "Validating orbit properties for\n"
         << "c interval: " << c << "\n"
         << "mu interval: " << mu << "\n"
         << "z1 interval: " << nbd << "\n";

    int existence = check_symmetric_existence_C0(mu, c, axis, sign, nbd);
    
    bool angle_monotonicity = true;
    for (int i = 0; i < Na; i++) 
    {
        interval zi = nbd.left() + diam(nbd)*interval(i, i+1) / Na;
        bool success = validate_angle_monotonicity(mu, c, axis, sign, zi, sign);  
        if (!success)
        {
            angle_monotonicity = false;
            break;
        }
    }

    bool nondegeneracy = true;
    for (int i = 0; i < Nn; i++) 
    {
        interval zi = nbd.left() + diam(nbd)*interval(i, i+1) / Nn;
        bool success = validate_nondegeneracy(mu, c, axis, sign, zi, k);    
        if (!success)
        {
            nondegeneracy = false;
            break;
        }
    }

    cout << "Results\n"
         << "Validation interval: " << nbd << "\n"
         << "Existence: " << existence << "\n"
         << "Angle monotonicity: " << angle_monotonicity << "\n"
         << "Nondegeneracy: " << nondegeneracy << "\n";

    return 0;
}
