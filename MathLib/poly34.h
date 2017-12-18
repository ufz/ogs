// poly34.h : solution of cubic and quartic equation
// (c) Khashin S.I. http://math.ivanovo.ac.ru/dalgebra/Khashin/index.html
// khash2 (at) gmail.com

// x - array of size 3
// return 3: 3 real roots x[0], x[1], x[2]
// return 1: 1 real root x[0] and pair of complex roots: x[1]張*x[2]
int   SolveP3(double *x, double a, double b, double c);			// solve cubic equation x^3 + a*x^2 + b*x + c = 0
