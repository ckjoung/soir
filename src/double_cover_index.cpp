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
    if (argc != 4)
    {
        cout << "Usage: orbit_properties <c> <z1> <nbd_radius>\n";
        return 0;
    }
    
    const interval c(argv[1], argv[1]);
    const interval mu = interval(99)/100; // 0.99
    double z1 = stod(argv[2]);
    const double r = stod(argv[3]);

    const int axis = 0; // z1-axis
    const int sign = -1;
    Shooting shooting(mu, c, axis, sign);
    interval nbd = interval(z1 - r, z1 + r);

    cout.precision(16);
    cout << "Calculating Conley-Zehnder index for\n"
         << "c interval: " << c << "\n"
         << "mu interval: " << mu << "\n"
         << "z1 interval: " << nbd << "\n";

    int existence = check_symmetric_existence_C0(mu, c, axis, sign, nbd);
    int index1, index2;
    try
    {
        index1 = computeCZIndex(mu, c, axis, sign, nbd, 1);
    }
    catch (std::exception& e)
    {
        index1 = -1;
    }
    try
    {
        index2 = computeCZIndex(mu, c, axis, sign, nbd, 2);
    }
    catch (std::exception& e)
    {
        index2 = -1;
    }
    bool nondegeneracy2 = validate_nondegeneracy(mu, c, axis, sign, nbd, 2);

    cout << "Results\n"
         << "Validation interval: " << nbd << "\n"
         << "Existence: " << existence << "\n"
         << "Index of half orbit: " << index1 << "\n"
         << "Index of full orbit: " << index2 << "\n"
         << "Non-degeneracy of full orbit: " << nondegeneracy2 << "\n";
         
    return 0;
}
