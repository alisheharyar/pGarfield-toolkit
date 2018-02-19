#ifndef LORGAMMA_H
#define LORGAMMA_H

/*
Functions for safe manipulations with beta and gamma kinematical
parameters without loss of precision at very extremal values,
such as beta very close to zero or to unity.

The main representation of gamma is gamma - 1, to assure that
gamma will be meaningful for small beta.

Author I. B. Smirnov, 1999 - 2002.
*/

namespace Heed {

double lorgamma_1(double beta);  // gamma - 1
double lorbeta(const double gamma_1);
double lorbeta2(const double gamma_1);
double lorbeta(const double momentum, const double mass);

}

#endif
