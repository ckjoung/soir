#ifndef UTILS_H
#define UTILS_H
#include <iostream>
#include <queue>
#include "capd/capdlib.h"
using namespace std;
using namespace capd;

// initialize queue by subdividing interval with given endpoints
void initialize_queue(queue<interval> &Q, double a, double b, double delta)
{
    double left = a;
    while (left < b) 
    {
        double right = std::min(left + delta, b);
        Q.push(interval(left, right));

        left += delta;
    }
}

// subdivide interval and push to queue
void subdivide_and_push(queue<interval> &Q, interval x, int n_subdivide=10)
{
    for (int i = 0; i < n_subdivide; i++) 
        Q.push(x.left() + diam(x) * interval(i, i+1) / n_subdivide);  
}

// sign corresponding to retrograde orbit
int sign_retrograde(int axis)
{
    if (axis == 0) // z1-axis
        return 1; 
    else // z2-axis
        return -1;
}

// coordinate of arriving slope
int slope_coordinate(int section)
{
    if (section == 0) // z1 = 0
        return 3; // w2 coordinate
    else // z2 = 0
        return 2; // w1 coordinate
}

// w coordinate for fixed point locus
int w_coordinate(int axis)
{
    if (axis == 0) // z1-axis
        return 3; // w2 coordinate
    else // z2-axis
        return 2; // w1 coordinate
}


#endif // UTILS_H