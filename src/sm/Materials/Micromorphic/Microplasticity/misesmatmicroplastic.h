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

#ifndef MisesMatMicroplastic_h

#include "../sm/Materials/misesmat.h"
#include "../sm/Materials/Micromorphic/micromorphicmaterialextensioninterface.h"
#include "../sm/Materials/Micromorphic/micromorphicms.h"
#include "../sm/Materials/linearelasticmaterial.h"

#include "cltypes.h"

///@name Input fields for MisesMatMicroplastic
//@{
#define _IFT_MisesMatMicroplastic_Name "misesmatmicroplastic"
//@}

namespace oofem {
/**
 * Gradient Mises material.
 */
class MisesMatMicroplastic : public MisesMat, MicromorphicMaterialExtensionInterface
{
protected:
 /// Reference to the basic elastic material.
// LinearElasticMaterial *linearElasticMaterial;


public:
    MisesMatMicroplastic(int n, Domain * d);
    virtual ~MisesMatMicroplastic();

    // definition
    virtual const char *giveInputRecordName() const { return _IFT_MisesMatMicroplastic_Name; }
    virtual const char *giveClassName() const { return "MisesMatMicroplastic"; }

    virtual IRResultType initializeFrom(InputRecord *ir);
    // virtual int hasMaterialModeCapability(MaterialMode mode);

    virtual Interface *giveInterface(InterfaceType t) {
        if ( t == MicromorphicMaterialExtensionInterfaceType ) {
            return static_cast< MicromorphicMaterialExtensionInterface * >(this);
        } else {
            return NULL;
        }
    }

    MaterialStatus* CreateStatus(GaussPoint *gp) const;

    virtual void giveStiffnessMatrix(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep);
    virtual void giveMicromorphicMatrix_dSigdUgrad(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep);
    virtual void giveMicromorphicMatrix_dSigdPhi(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep);
    virtual void giveMicromorphicMatrix_dSdUgrad(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep);
    virtual void giveMicromorphicMatrix_dSdPhi(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep);
    virtual void giveMicromorphicMatrix_dMdPhiGrad(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep);

    virtual void giveGeneralizedStressVectors (FloatArray &sigma, FloatArray &s, FloatArray &S, GaussPoint *gp, const FloatArray &totalStrain, const FloatArray &micromorphicVar, const FloatArray micromorphicVarGrad, TimeStep *tStep);
    virtual void giveFiniteStrainGeneralizedStressVectors (FloatArray &sigma, FloatArray &s, FloatArray &M, GaussPoint *gp, const FloatArray &displacementGradient, const FloatArray &micromorphicVar, const FloatArray micromorphicVarGrad, TimeStep *tStep){;}
    
    void performPlasticityReturn(GaussPoint *gp, const FloatArray &totalStrain, const FloatArray &micromorphicVar, const FloatArray &micromorphicVarGrad);

    LinearElasticMaterial *giveLinearElasticMaterial() { return linearElasticMaterial; }

    virtual bool isStressTensorSymmetric(){return true;}
    

protected:
    //   virtual MaterialStatus *CreateStatus(GaussPoint *gp) const { return new MisesMatMicroplasticStatus(1, MicromorphicMaterialStatus :: domain, gp); }
};


/**
 * Gradient Mises maaterial status.
 */
 class MisesMatMicroplasticStatus : public MicromorphicMaterialStatus
{

    /// Plastic strain (initial).
    FloatArray plasticStrain;

    /// Plastic strain (final).
    FloatArray tempPlasticStrain;

    /// Deviatoric trial stress - needed for tangent stiffness.
    FloatArray trialStressDev;

    /// Cumulative plastic strain (initial).
    double kappa;

    /// Cumulative plastic strain (final).
    double tempKappa;

public:
    MisesMatMicroplasticStatus(int n, Domain * d, GaussPoint * g);
    virtual ~MisesMatMicroplasticStatus();


    void letTrialStressDevBe(const FloatArray &values) { trialStressDev = values; }
    const FloatArray &giveTrialStressDev() { return trialStressDev; }


    void setTempCumulativePlasticStrain(double value) { tempKappa = value; }

    double giveCumulativePlasticStrain() { return kappa; }
    double giveTempCumulativePlasticStrain() { return tempKappa; }


   void letTempPlasticStrainBe(const FloatArray &values) { tempPlasticStrain = values; }
    
    const FloatArray &givePlasticStrain() { return plasticStrain; }
    const FloatArray &getTempPlasticStrain() const { return tempPlasticStrain; }

    virtual void initTempStatus();

    virtual void updateYourself(TimeStep *tStep);

    //   virtual void printOutputAt(FILE *file, TimeStep *tStep);

    virtual const char *giveClassName() const { return "MisesMatMicroplasticStatus"; }


};




} // end namespace oofem
#define MisesMatMicroplastic_h
#endif
