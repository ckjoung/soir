#ifndef RTBP_H
#define RTBP_H
#include <iostream>
using namespace std;
using namespace capd;
using namespace capd::autodiff;

/**
 * Derivatives of the Levi-Civita Hamiltonian $K_{\mu,c}$ defined in formula (2.2).
 * The action component is obtained as $\lambda(X_K)$, where we have 
 * $$ \lambda = \sum_{i=1}^2 (w_i dz_i -z_i dw_i)$$.
 *  
 */
std::string LCRTBPStr(bool includeAction=false) 
{
    string z1Eq = "w1-(2*z1^2+2*z2^2)*z2-mu*z2";
    string z2Eq = "w2+(2*z1^2+2*z2^2)*z1-mu*z1";
    string w1Eq = "-2*c*z1-4*z1*(-w1*z2+w2*z1)-(2*z1^2+2*z2^2)*w2+mu*w2+2*mu*z1/sqrt((2*z1^2-2*z2^2-1)^2+16*z1^2*z2^2)-(1/2)*mu*(z1^2+z2^2)*((8*(2*z1^2-2*z2^2-1))*z1+32*z1*z2^2)/((2*z1^2-2*z2^2-1)^2+16*z1^2*z2^2)^(3/2)";
    string w2Eq = "-2*c*z2-4*z2*(-w1*z2+w2*z1)+(2*z1^2+2*z2^2)*w1+mu*w1+2*mu*z2/sqrt((2*z1^2-2*z2^2-1)^2+16*z1^2*z2^2)-(1/2)*mu*(z1^2+z2^2)*(-(8*(2*z1^2-2*z2^2-1))*z2+32*z1^2*z2)/((2*z1^2-2*z2^2-1)^2+16*z1^2*z2^2)^(3/2)";
    if (!includeAction)
        return "par:mu,c;time:t;var:z1,z2,w1,w2;fun:" + z1Eq + "," + z2Eq + "," + w1Eq + "," + w2Eq + ";";
    
    string lEq = "w1*(w1-(2*z1^2+2*z2^2)*z2-mu*z2)+w2*(w2+(2*z1^2+2*z2^2)*z1-mu*z1)-z1*(-2*c*z1-4*z1*(-w1*z2+w2*z1)-(2*z1^2+2*z2^2)*w2+mu*w2+2*mu*z1/sqrt((2*z1^2-2*z2^2-1)^2+16*z1^2*z2^2)-(1/2)*mu*(z1^2+z2^2)*((8*(2*z1^2-2*z2^2-1))*z1+32*z1*z2^2)/((2*z1^2-2*z2^2-1)^2+16*z1^2*z2^2)^(3/2))-z2*(-2*c*z2-4*z2*(-w1*z2+w2*z1)+(2*z1^2+2*z2^2)*w1+mu*w1+2*mu*z2/sqrt((2*z1^2-2*z2^2-1)^2+16*z1^2*z2^2)-(1/2)*mu*(z1^2+z2^2)*(-(8*(2*z1^2-2*z2^2-1))*z2+32*z1^2*z2)/((2*z1^2-2*z2^2-1)^2+16*z1^2*z2^2)^(3/2))";
    return "par:mu,c;time:t;var:z1,z2,w1,w2,l;fun:" + z1Eq + "," + z2Eq + "," + w1Eq + "," + w2Eq + "," + lEq + ";";
}

#endif // RTBP_H