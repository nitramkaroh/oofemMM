/*
 *
 *                 #####    #####   ######  ######  ###   ###
 *               ##   ##  ##   ##  ##      ##      ## ### ##
 *              ##   ##  ##   ##  ####    ####    ##  #  ##
 *             ##   ##  ##   ##  ##      ##      ##     ##
 *            ##   ##  ##   ##  ##      ##      ##     ##
 *            #####    #####   ##      ######  ##     ##
 *
 *
 *             OOFEM : Object Oriented Finite Element Code
 *
 *               Copyright (C) 1993 - 2013   Borek Patzak
 *
 *
 *
 *       Czech Technical University, Faculty of Civil Engineering,
 *   Department of Structural Mechanics, 166 29 Prague, Czech Republic
 *
 *  This library is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU Lesser General Public
 *  License as published by the Free Software Foundation; either
 *  version 2.1 of the License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *  Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public
 *  License along with this library; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 */

#ifndef dofiditemh
#define dofiditemh

#include "enumitem.h"

#include <string>

namespace oofem {
#define DofIDItem_DEF \
    ENUM_ITEM_WITH_VALUE(Undef, 0) /* Error value */ \
    ENUM_ITEM_WITH_VALUE(D_u, 1) /* u-displacement (in direction of x-axis) */ \
    ENUM_ITEM_WITH_VALUE(D_v, 2) /* v-displacement (in direction of y-axis) */ \
    ENUM_ITEM_WITH_VALUE(D_w, 3) /* w-displacement (in direction of z-axis) */ \
    ENUM_ITEM_WITH_VALUE(R_u, 4) /* Rotation around x-axis (right hand rule assumed) */ \
    ENUM_ITEM_WITH_VALUE(R_v, 5) /* Rotation around y-axis */ \
    ENUM_ITEM_WITH_VALUE(R_w, 6) /* Rotation around z-axis */ \
  \
    ENUM_ITEM_WITH_VALUE(V_u, 7) /* u-velocity (in direction of x-axis) */ \
    ENUM_ITEM_WITH_VALUE(V_v, 8) /* v-velocity (in direction of y-axis) */ \
    ENUM_ITEM_WITH_VALUE(V_w, 9) /* w-velocity (in direction of z-axis) */ \
  \
    ENUM_ITEM_WITH_VALUE(T_f, 10) /* Temperature field */ \
    ENUM_ITEM_WITH_VALUE(P_f, 11) /* Pressure field */ \
    ENUM_ITEM_WITH_VALUE(G_0, 12) /* DOF for gradient formulation no. 0 */ \
    ENUM_ITEM_WITH_VALUE(G_1, 13) /* DOF for gradient formulation no. 1 */ \
    ENUM_ITEM_WITH_VALUE(C_1, 14) /* Mass concentration of the first constituent */ \
    ENUM_ITEM_WITH_VALUE(W_u, 15) /* u-component of change in director field (in direction of x-axis) */ \
    ENUM_ITEM_WITH_VALUE(W_v, 16) /* v-component of change in director field (in direction of y-axis) */ \
    ENUM_ITEM_WITH_VALUE(W_w, 17) /* w-component of change in director field (in direction of z-axis) */ \
    ENUM_ITEM_WITH_VALUE(Gamma, 18) /* inhomogenous thickness strain in direction of the directorfield m */ \
      /*Micromorphic degrees of freedom, i.e. components of full micromorphic displacement gradient tensor */ \
      ENUM_ITEM_WITH_VALUE(M_X11, 19)					\
      ENUM_ITEM_WITH_VALUE(M_X22, 20)					\
      ENUM_ITEM_WITH_VALUE(M_X33, 21)					\
      ENUM_ITEM_WITH_VALUE(M_X23, 22)					\
      ENUM_ITEM_WITH_VALUE(M_X13, 23)					\
      ENUM_ITEM_WITH_VALUE(M_X12, 24)					\
      ENUM_ITEM_WITH_VALUE(M_X32, 25)					\
      ENUM_ITEM_WITH_VALUE(M_X31, 26)					\
      ENUM_ITEM_WITH_VALUE(M_X21, 27)					\
      /*Microstretch degrees of freedom, i.e. components of symetric part of  micromorphic displacement gradient tensor */ \
      ENUM_ITEM_WITH_VALUE(M_E11, 28)					\
      ENUM_ITEM_WITH_VALUE(M_E22, 29)					\
      ENUM_ITEM_WITH_VALUE(M_E33, 30)					\
      ENUM_ITEM_WITH_VALUE(M_E23, 31)					\
      ENUM_ITEM_WITH_VALUE(M_E13, 32)					\
      ENUM_ITEM_WITH_VALUE(M_E12, 33)					\
      /*Micropolar(Cosserat) degrees of freedom, i.e. components of skew part of  micromorphic displacement gradient tensor */ \
      ENUM_ITEM_WITH_VALUE(M_W1, 34)					\
      ENUM_ITEM_WITH_VALUE(M_W2, 35)					\
      ENUM_ITEM_WITH_VALUE(M_W3, 36)					\
      /*Microstrain norm degree of freedom, i.e. trace of micromorphic displacement gradient(or other equivalent micromorphic "strain" norm) */ \
      ENUM_ITEM_WITH_VALUE(M_D, 37)					\
      /*Microstrain penalty degree of freedom, within this formulation, the srain gradient or coupled stress theory can be implemented with C0 elements*/ \
      ENUM_ITEM_WITH_VALUE(M_P, 38)					\
      /*Micromorphic cumulative strain degree of freedom, within this formulation, the Aifantis gradient plasticity can be implemented with C0 elements*/ \
      ENUM_ITEM_WITH_VALUE(M_MP, 39)					\
      /*Warping dofs*/							\
      ENUM_ITEM_WITH_VALUE(Warp_PsiTheta, 40) /* Relative twist times deplanation function, used by Trwarp element*/ \
      ENUM_ITEM_WITH_VALUE(Warp_Theta, 41) /* Relative twist, used by Trwarp element*/ \
      /*Hybrid-Trefftz dofs - old ready to remove*/			\
      ENUM_ITEM_WITH_VALUE(D_u_edge_const, 42) /* Constant part of boundary u-displacement used by Trefftz element*/ \
      ENUM_ITEM_WITH_VALUE(D_u_edge_lin, 43) /* Linear part of boundary u-displacement used by Trefftz element*/ \
      ENUM_ITEM_WITH_VALUE(D_v_edge_const, 44) /* Constant part of boundary v-displacement used by Trefftz element*/ \
      ENUM_ITEM_WITH_VALUE(D_v_edge_lin, 45) /* Linear part of boundary v-displacement used by Trefftz element*/ \


/**
 * Type representing particular dof type. Values of this type describe the physical meaning of
 * available DOFs.
 * @note{The implementation of Node::computeGNTransformation rely on D_u, D_v and D_w (R_u, R_v, R_w) order.
 * Do not change their order and do not insert any values between these values.}
 */
enum DofIDItem {
    DofIDItem_DEF
    MaxDofID = 500
};

#undef ENUM_ITEM
#undef ENUM_ITEM_WITH_VALUE
#undef enumitem_h

std :: string __DofIDItemToString(DofIDItem _value);
} // end namespace oofem
#endif // dofiditem_h
