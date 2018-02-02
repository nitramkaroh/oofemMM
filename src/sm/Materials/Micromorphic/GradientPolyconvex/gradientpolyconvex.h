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
 *               Copyright (C) 1993 - 2015   Borek Patzak
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

#ifndef gradientpolyconvex_h
#define gradientpolyconvex_h


#include "../sm/Materials/Micromorphic/micromorphicmaterialextensioninterface.h"
#include "../sm/Materials/Micromorphic/micromorphicms.h"
#include "../sm/Materials/isolinearelasticmaterial.h"
#include "../sm/Materials/structuralmaterial.h"
#include "cltypes.h"


///@name Input fields for MicromorphLEmat
//@{
#define _IFT_GradientPolyconvexMaterial_Name "gradientpolyconvexmaterial"
//#define _IFT_MicromorphicMaterial_Hk "hk"
//#define _IFT_MicromorphicMaterial_Ak "ak"
//@}

namespace oofem {


/**
 * MicromorphicLinearElasticMaterial
 */
class GradientPolyconvexMaterial : public IsotropicLinearElasticMaterial, MicromorphicMaterialExtensionInterface
{
protected:

public:
    GradientPolyconvexMaterial(int n, Domain * d);
    virtual ~GradientPolyconvexMaterial();

    // definition
    virtual const char *giveInputRecordName() const { return _IFT_GradientPolyconvexMaterial_Name; }
    virtual const char *giveClassName() const { return "GradientPolyconvexMaterial"; }

    virtual IRResultType initializeFrom(InputRecord *ir);

    virtual Interface *giveInterface(InterfaceType t) {
        if ( t == MicromorphicMaterialExtensionInterfaceType ) {
            return static_cast< MicromorphicMaterialExtensionInterface * >(this);
        } else {
            return NULL;
        }
    }

    virtual void giveStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep);
    virtual void giveMicromorphicMatrix_dSigdUgrad(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep);
    virtual void giveMicromorphicMatrix_dSigdPhi(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep);
    virtual void giveMicromorphicMatrix_dSdUgrad(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep);
    virtual void giveMicromorphicMatrix_dSdPhi(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep);
    virtual void giveMicromorphicMatrix_dMdPhiGrad(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep);

    virtual void giveGeneralizedStressVectors (FloatArray &sigma, FloatArray &s, FloatArray &S, GaussPoint *gp, const FloatArray &totalStrain, const FloatArray &micromorphicVar, const FloatArray micromorphicVarGrad, TimeStep *tStep);


    virtual int giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep);

protected:
  virtual MaterialStatus *CreateStatus(GaussPoint *gp) const { return new MicromorphicMaterialStatus(1, domain, gp); }
                                                                            
};


} // end namespace oofem
#endif
