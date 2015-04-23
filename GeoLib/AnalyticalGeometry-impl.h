/**
 * \copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

namespace GeoLib
{

template <class T_POINT>
void getNewellPlane (const std::vector<T_POINT*>& pnts,
                     MathLib::Vector3 &plane_normal,
                     double& d)
{
    d = 0;
    MathLib::Vector3 centroid;
    std::size_t n_pnts(pnts.size());
    for (size_t i(n_pnts - 1), j(0); j < n_pnts; i = j, j++) {
        plane_normal[0] += ((*(pnts[i]))[1] - (*(pnts[j]))[1])
                           * ((*(pnts[i]))[2] + (*(pnts[j]))[2]); // projection on yz
        plane_normal[1] += ((*(pnts[i]))[2] - (*(pnts[j]))[2])
                           * ((*(pnts[i]))[0] + (*(pnts[j]))[0]); // projection on xz
        plane_normal[2] += ((*(pnts[i]))[0] - (*(pnts[j]))[0])
                           * ((*(pnts[i]))[1] + (*(pnts[j]))[1]); // projection on xy

        centroid += *(pnts[j]);
    }

    plane_normal.normalize();
    d = MathLib::scalarProduct(centroid, plane_normal) / n_pnts;
}

template<class T_MATRIX>
void compute2DRotationMatrixToX(MathLib::Vector3 const& v,
        T_MATRIX & rot_mat)
{
    const double cos_theta = v[0];
    const double sin_theta = v[1];
    rot_mat(0,0) = rot_mat(1,1) = cos_theta;
    rot_mat(0,1) = sin_theta;
    rot_mat(1,0) = -sin_theta;
}

template <class T_MATRIX>
void compute3DRotationMatrixToX(MathLib::Vector3  const& v,
        T_MATRIX & rot_mat)
{
    // a vector on the plane
    MathLib::Vector3 yy(0, 0, 0);
    if (fabs(v[0]) > 0.0 && fabs(v[1]) + fabs(v[2]) < std::numeric_limits<double>::epsilon()) {
        yy[2] = 1.0;
    } else if (fabs(v[1]) > 0.0 && fabs(v[0]) + fabs(v[2]) < std::numeric_limits<double>::epsilon()) {
        yy[0] = 1.0;
    } else if (fabs(v[2]) > 0.0 && fabs(v[0]) + fabs(v[1]) < std::numeric_limits<double>::epsilon()) {
        yy[1] = 1.0;
    } else {
        for (unsigned i = 0; i < 3; i++) {
            if (fabs(v[i]) > 0.0) {
                yy[i] = -v[i];
                break;
            }
        }
    }
    // z"_vec
    MathLib::Vector3 zz(MathLib::crossProduct(v, yy));
    zz.normalize();
    // y"_vec
    yy = MathLib::crossProduct(zz, v);
    yy.normalize();

    for (unsigned i=0; i<3; ++i) {
        rot_mat(0, i) = v[i];
        rot_mat(1, i) = yy[i];
        rot_mat(2, i) = zz[i];
    }
}

template <class T_MATRIX>
void computeRotationMatrixToXY2(const std::vector<MathLib::Point3d*> &vec_pt,
        T_MATRIX & rot_mat)
{
    assert(vec_pt.size()>2);
    // x"_vec
    MathLib::Vector3 xx((*vec_pt[0]), (*vec_pt[1]));
    xx.normalize();
    // a vector on the plane
    MathLib::Vector3 yy((*vec_pt[1]), (*vec_pt[2]));
    // z"_vec
    MathLib::Vector3 zz(MathLib::crossProduct(xx, yy));
    zz.normalize();
    // y"_vec
    yy = MathLib::crossProduct(zz, xx);
    yy.normalize();

    for (size_t i=0; i<3; ++i) {
        rot_mat(0, i) = xx[i];
        rot_mat(1, i) = yy[i];
        rot_mat(2, i) = zz[i];
    }
}

} // end namespace GeoLib

