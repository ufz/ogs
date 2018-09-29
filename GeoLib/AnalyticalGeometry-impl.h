/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

namespace GeoLib
{

template <typename InputIterator>
void getNewellPlane (InputIterator pnts_begin, InputIterator pnts_end,
                     MathLib::Vector3 &plane_normal,
                     double& d)
{
    MathLib::Vector3 centroid;
    for (auto i=std::prev(pnts_end), j=pnts_begin; j!=pnts_end; i = j, ++j) {
        auto &pt_i = *(*i);
        auto &pt_j = *(*j);
        plane_normal[0] += (pt_i[1] - pt_j[1])
                           * (pt_i[2] + pt_j[2]); // projection on yz
        plane_normal[1] += (pt_i[2] - pt_j[2])
                           * (pt_i[0] + pt_j[0]); // projection on xz
        plane_normal[2] += (pt_i[0] - pt_j[0])
                           * (pt_i[1] + pt_j[1]); // projection on xy

        centroid += pt_j;
    }

    plane_normal.normalize();
    auto n_pnts(std::distance(pnts_begin, pnts_end));
    assert(n_pnts > 2);
    if (n_pnts > 0)
        d = MathLib::scalarProduct(centroid, plane_normal) / n_pnts;
}

template <class T_POINT>
void getNewellPlane (const std::vector<T_POINT*>& pnts,
                     MathLib::Vector3 &plane_normal,
                     double& d)
{
    getNewellPlane(pnts.begin(), pnts.end(), plane_normal, d);
}

template <class T_POINT>
std::pair<MathLib::Vector3, double> getNewellPlane (const std::vector<T_POINT>& pnts)
{
    MathLib::Vector3 plane_normal;
    MathLib::Vector3 centroid;
    std::size_t n_pnts(pnts.size());
    for (std::size_t i = n_pnts - 1, j = 0; j < n_pnts; i = j, j++) {
        plane_normal[0] += (pnts[i][1] - pnts[j][1])
                           * (pnts[i][2] + pnts[j][2]); // projection on yz
        plane_normal[1] += (pnts[i][2] - pnts[j][2])
                           * (pnts[i][0] + pnts[j][0]); // projection on xz
        plane_normal[2] += (pnts[i][0] - pnts[j][0])
                           * (pnts[i][1] + pnts[j][1]); // projection on xy

        centroid += pnts[j];
    }

    plane_normal.normalize();
    double d = MathLib::scalarProduct(centroid, plane_normal) / n_pnts;
    return std::make_pair(plane_normal, d);
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
    rot_mat(0,2) = rot_mat(1,2) = rot_mat(2,0) = rot_mat(2,1) = 0.0;
    rot_mat(2,2) = 1.0;
}

template <class T_MATRIX>
void compute3DRotationMatrixToX(MathLib::Vector3  const& v,
        T_MATRIX & rot_mat)
{
    // a vector on the plane
    MathLib::Vector3 yy(0, 0, 0);
    if (std::abs(v[0]) > 0.0 && std::abs(v[1]) + std::abs(v[2]) < std::numeric_limits<double>::epsilon()) {
        yy[2] = 1.0;
    } else if (std::abs(v[1]) > 0.0 && std::abs(v[0]) + std::abs(v[2]) < std::numeric_limits<double>::epsilon()) {
        yy[0] = 1.0;
    } else if (std::abs(v[2]) > 0.0 && std::abs(v[0]) + std::abs(v[1]) < std::numeric_limits<double>::epsilon()) {
        yy[1] = 1.0;
    } else {
        for (unsigned i = 0; i < 3; i++) {
            if (std::abs(v[i]) > 0.0) {
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
void computeRotationMatrixToXY(MathLib::Vector3 const& n,
        T_MATRIX & rot_mat)
{
    // check if normal points already in the right direction
    if (n[0] == 0 && n[1] == 0) {
        rot_mat(0,1) = 0.0;
        rot_mat(0,2) = 0.0;
        rot_mat(1,0) = 0.0;
        rot_mat(1,1) = 1.0;
        rot_mat(1,2) = 0.0;
        rot_mat(2,0) = 0.0;
        rot_mat(2,1) = 0.0;

        if (n[2] > 0) {
            // identity matrix
            rot_mat(0,0) = 1.0;
            rot_mat(2,2) = 1.0;
        } else {
            // rotate by pi about the y-axis
            rot_mat(0,0) = -1.0;
            rot_mat(2,2) = -1.0;
        }

        return;
    }

    // sqrt (n_1^2 + n_3^2)
    double const h0(sqrt(n[0]*n[0]+n[2]*n[2]));

    // In case the x and z components of the normal are both zero the rotation
    // to the x-z-plane is not required, i.e. only the rotation in the z-axis is
    // required. The angle is either pi/2 or 3/2*pi. Thus the components of
    // rot_mat are as follows.
    if (h0 < std::numeric_limits<double>::epsilon()) {
        rot_mat(0,0) = 1.0;
        rot_mat(0,1) = 0.0;
        rot_mat(0,2) = 0.0;
        rot_mat(1,0) = 0.0;
        rot_mat(1,1) = 0.0;
        if (n[1] > 0)
            rot_mat(1,2) = -1.0;
        else
            rot_mat(1,2) = 1.0;
        rot_mat(2,0) = 0.0;
        if (n[1] > 0)
            rot_mat(2,1) = 1.0;
        else
            rot_mat(2,1) = -1.0;
        rot_mat(2,2) = 0.0;
        return;
    }

    double h1(1 / n.getLength());

    // general case: calculate entries of rotation matrix
    rot_mat(0, 0) = n[2] / h0;
    rot_mat(0, 1) = 0;
    rot_mat(0, 2) = - n[0] / h0;
    rot_mat(1, 0) = - n[1]*n[0]/h0 * h1;
    rot_mat(1, 1) = h0 * h1;
    rot_mat(1, 2) = - n[1]*n[2]/h0 * h1;
    rot_mat(2, 0) = n[0] * h1;
    rot_mat(2, 1) = n[1] * h1;
    rot_mat(2, 2) = n[2] * h1;
}

template <typename InputIterator>
void rotatePoints(
        MathLib::DenseMatrix<double> const& rot_mat,
        InputIterator pnts_begin, InputIterator pnts_end)
{
    double* tmp (nullptr);
    for (auto it=pnts_begin; it!=pnts_end; ++it) {
        tmp = rot_mat * (*it)->getCoords();
        for (std::size_t j(0); j < 3; j++)
            (*(*it))[j] = tmp[j];
        delete [] tmp;
    }
}

template <typename InputIterator1, typename InputIterator2>
MathLib::DenseMatrix<double> rotatePointsToXY(InputIterator1 p_pnts_begin,
                                              InputIterator1 p_pnts_end,
                                              InputIterator2 r_pnts_begin,
                                              InputIterator2 r_pnts_end)
{
    assert(std::distance(p_pnts_begin, p_pnts_end) > 2);

    // calculate supporting plane
    MathLib::Vector3 plane_normal;
    double d;
    // compute the plane normal
    GeoLib::getNewellPlane(p_pnts_begin, p_pnts_end, plane_normal, d);

    // rotate points into x-y-plane
    MathLib::DenseMatrix<double> rot_mat(3, 3);
    computeRotationMatrixToXY(plane_normal, rot_mat);
    rotatePoints(rot_mat, r_pnts_begin, r_pnts_end);

    for (auto it=r_pnts_begin; it!=r_pnts_end; ++it)
        (*(*it))[2] = 0.0; // should be -= d but there are numerical errors

    return rot_mat;
}

template <typename P>
void rotatePoints(MathLib::DenseMatrix<double> const& rot_mat,
    std::vector<P*> const& pnts)
{
    rotatePoints(rot_mat, pnts.begin(), pnts.end());
}

} // end namespace GeoLib

