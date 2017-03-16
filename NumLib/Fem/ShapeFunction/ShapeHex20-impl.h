/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

namespace
{

inline double ShapeFunctionHexHQ_Corner(const double r, const double s, const double t)
{
    return 0.125 * (r - 1) * (s - 1) * (t - 1) * (r + s + t + 2.0);
}

inline double ShapeFunctionHexHQ_Middle(const double r, const double s, const double t)
{
    return 0.25 * (1 - r * r) * (s - 1) * (t - 1);
}

inline double dShapeFunctionHexHQ_Corner(const double r, const double s, const double t, const int ty)
{
    switch(ty)
    {
    case 0:
        return 0.125 * (s - 1) * (t - 1) * (2.0 * r + s + t + 1.0);
    case 1:
        return 0.125 * (t - 1) * (r - 1) * (2.0 * s + r + t + 1.0);
    case 2:
        return 0.125 * (r - 1) * (s - 1) * (2.0 * t + s + r + 1.0);
    }
    return 0.0;
}

inline double dShapeFunctionHexHQ_Middle(const double r, const double s, const double t, const int ty)
{
    switch(ty)
    {
    case 0:
        return -0.5 * r * (s - 1) * (t - 1);
    case 1:
        return 0.25 * (1 - r * r) * (t - 1);
    case 2:
        return 0.25 * (1 - r * r) * (s - 1);
    }
    return 0.0;
}

}

namespace NumLib
{

template <class T_X, class T_N>
void ShapeHex20::computeShapeFunction(const T_X &rst, T_N &N)
{
    const double r = rst[0];
    const double s = rst[1];
    const double t = rst[2];

    N[0] = ShapeFunctionHexHQ_Corner(r,s,t);
    N[1] = ShapeFunctionHexHQ_Corner(-r,s,t);
    N[2] = ShapeFunctionHexHQ_Corner(-r,-s,t);
    N[3] = ShapeFunctionHexHQ_Corner(r,-s,t);
    N[4] = ShapeFunctionHexHQ_Corner(r,s,-t);
    N[5] = ShapeFunctionHexHQ_Corner(-r,s,-t);
    N[6] = ShapeFunctionHexHQ_Corner(-r,-s,-t);
    N[7] = ShapeFunctionHexHQ_Corner(r,-s,-t);

    N[8] = ShapeFunctionHexHQ_Middle(r,s,t);
    N[10] = ShapeFunctionHexHQ_Middle(r,-s,t);
    N[14] = ShapeFunctionHexHQ_Middle(r,-s,-t);
    N[12] = ShapeFunctionHexHQ_Middle(r,s,-t);

    N[11] = ShapeFunctionHexHQ_Middle(s,t,r);
    N[15] = ShapeFunctionHexHQ_Middle(s,-t,r);
    N[13] = ShapeFunctionHexHQ_Middle(s,-t,-r);
    N[9]  =  ShapeFunctionHexHQ_Middle(s,t,-r);

    N[16] = ShapeFunctionHexHQ_Middle(t,r,s);
    N[17] = ShapeFunctionHexHQ_Middle(t,-r,s);
    N[18] = ShapeFunctionHexHQ_Middle(t,-r,-s);
    N[19] = ShapeFunctionHexHQ_Middle(t,r,-s);
}

template <class T_X, class T_N>
void ShapeHex20::computeGradShapeFunction(const T_X &rst, T_N &dNdr)
{
    int co;
    const double r = rst[0];
    const double s = rst[1];
    const double t = rst[2];
    const static double sign1[] = {-1.0, 1.0,1.0};
    const static double sign2[] = { 1.0,-1.0,1.0};
    const static double sign3[] = { 1.0, 1.0,-1.0};
    for(int i = 0; i < 3; i++)
    {
        dNdr[20 * i + 0] = dShapeFunctionHexHQ_Corner(r,s,t,i);
        dNdr[20 * i + 1] = sign1[i] * dShapeFunctionHexHQ_Corner(-r,s,t,i);
        dNdr[20 * i + 2] = sign1[i] * sign2[i] * dShapeFunctionHexHQ_Corner(-r,-s,t,i);
        dNdr[20 * i + 3] = sign2[i] * dShapeFunctionHexHQ_Corner(r,-s,t,i);
        dNdr[20 * i + 4] = sign3[i] * dShapeFunctionHexHQ_Corner(r,s,-t,i);
        dNdr[20 * i + 5] = sign1[i] * sign3[i] * dShapeFunctionHexHQ_Corner(-r,s,-t,i);
        dNdr[20 * i + 6] = sign1[i] * sign2[i] * sign3[i] * dShapeFunctionHexHQ_Corner(-r,-s,-t,i);
        dNdr[20 * i + 7] = sign2[i] * sign3[i] * dShapeFunctionHexHQ_Corner(r,-s,-t,i);

        dNdr[20 * i + 8] =  dShapeFunctionHexHQ_Middle(r,s,t,i);
        dNdr[20 * i + 10] = sign2[i] * dShapeFunctionHexHQ_Middle(r,-s,t,i);
        dNdr[20 * i + 14] = sign2[i] * sign3[i] * dShapeFunctionHexHQ_Middle(r,-s,-t,i);
        dNdr[20 * i + 12] = sign3[i] * dShapeFunctionHexHQ_Middle(r,s,-t,i);

        co = (i + 2) % 3;
        dNdr[20 * i + 11] = dShapeFunctionHexHQ_Middle(s,t,r,co);
        dNdr[20 * i + 15] = sign3[i] * dShapeFunctionHexHQ_Middle(s,-t,r,co);
        dNdr[20 * i + 13] = sign1[i] * sign3[i] * dShapeFunctionHexHQ_Middle(s,-t,-r,co);
        dNdr[20 * i + 9] =  sign1[i] * dShapeFunctionHexHQ_Middle(s,t,-r,co);

        co = (i + 1) % 3;
        dNdr[20 * i + 16] = dShapeFunctionHexHQ_Middle(t,r,s,co);
        dNdr[20 * i + 17] = sign1[i] * dShapeFunctionHexHQ_Middle(t,-r,s,co);
        dNdr[20 * i + 18] = sign1[i] * sign2[i] * dShapeFunctionHexHQ_Middle(t,-r,-s,co);
        dNdr[20 * i + 19] = sign2[i] * dShapeFunctionHexHQ_Middle(t,r,-s,co);
    }
}

}

