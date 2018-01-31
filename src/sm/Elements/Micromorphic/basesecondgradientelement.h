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
 *
 *               Copyright (C) 1993 - 2013   Borek Patzak
 *
 *
 *
 *       Czech Technical University, Faculty of Civil Engineering,
 *   Department of Structural Mechanics, 166 29 Prague, Czech Republic
 *
 *  This library is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU Lesser Base Public
 *  License as published by the Free Software Foundation; either
 *  version 2.1 of the License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *  Lesser Base Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser Base Public
 *  License along with this library; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 */

#ifndef micromorphicelement_h
#define micromorphicelement_h

#include "../sm/Elements/structuralelement.h"
#include "../sm/Elements/nlstructuralelement.h"

namespace oofem {
/**
 * Base class for second gradient continua.
 * This class contains constrained microstrain, microstretch, Couple-stress, or microdilatational continua
 * Regarding the specific formulation, new micro degrees of freedom are introduced
 * Moreover, lagrangian multiplier is introduce to constrain the micromorphic dofs to obtain the second gradient continue
 * @author Martin Horak A simplified second gradient model for dilatant materials: Theory and numerical implementation
 */
class BaseSecondGradientElement
{
protected:
    IntArray displacementDofsOrdering, micromorphicDofsOrdering, lagrangianMultiplierDofsOrdering;
    bool isStressTensorSymmetric;
    IntArray locationArray_u, locationArray_m, locationArray_l;
    

public:
    BaseSecondGradientElement();
    virtual ~BaseSecondGradientElement() { }

    virtual IRResultType initializeFrom(InputRecord *ir);

protected:
  
    /// Pure virtual functions
    virtual NLStructuralElement *giveElement() = 0;
    void computeB_uMatrixAt(GaussPoint *gp, FloatMatrix &B, NLStructuralElement *element, bool isStressTensorSymmetric);

    virtual void computeMicromorphicNMatrixAt(GaussPoint *, FloatMatrix &) = 0;
    virtual void computeMicromorphicBMatrixAt(GaussPoint *, FloatMatrix &) = 0;

    virtual void computeLagrangianMultiplierNMatrixAt(GaussPoint *, FloatMatrix &) = 0;

    virtual int giveNumberOfMicromorphicDofs() = 0;
    virtual int giveNumberOfDisplacementDofs() = 0;
    virtual int giveNumberOfLagrangianMultipliersDofs() = 0;
    virtual int giveNumberOfDofs() = 0;

    virtual void giveDofManDofIDMask_u(IntArray &answer) = 0;
    virtual void giveDofManDofIDMask_m(IntArray &answer) = 0;
    virtual void giveDofManDofIDMask_l(IntArray &answer) = 0;
    /// End of pure virtual functions


    bool symmetricFormulation = false;
  

   virtual void computeStiffnessMatrix(FloatMatrix &, MatResponseMode, TimeStep *);

    void computeStiffnessMatrix_uu(FloatMatrix &, MatResponseMode, TimeStep *);
    void computeStiffnessMatrix_ul(FloatMatrix &, MatResponseMode, TimeStep *);

    void computeStiffnessMatrix_mm(FloatMatrix &, MatResponseMode, TimeStep *);
    void computeStiffnessMatrix_ml(FloatMatrix &, MatResponseMode, TimeStep *);

    void computeStiffnessMatrix_lu(FloatMatrix &, MatResponseMode, TimeStep *);
    void computeStiffnessMatrix_lm(FloatMatrix &, MatResponseMode, TimeStep *);




    void computeMassMatrix(FloatMatrix &answer, TimeStep *tStep);

    void computeMassMatrix_uu(FloatMatrix &answer, TimeStep *tStep);
    void computeMassMatrix_mm(FloatMatrix &answer, TimeStep *tStep);
    void computeMassMatrix_ll(FloatMatrix &answer, TimeStep *tStep);




    void computeGeneralizedStressVectors(FloatArray &sigma, FloatArray &s, FloatArray &M, FloatArray &relativeStrain, GaussPoint *gp, TimeStep *tStep);
    void computeDisplacementGradient(FloatArray &answer, GaussPoint *gp, TimeStep *tStep, bool isStressTensorSymmetric);
    void computeMicromorphicVars(FloatArray &micromorphVar, FloatArray &micromorphVarGrad, IntArray IdMask_m, GaussPoint *gp, TimeStep *tStep);   
    void computeLagrangianMultipliers(FloatArray &lagrangianMulriplier, IntArray IdMask_l, GaussPoint *gp, TimeStep *tStep);   


    void giveInternalForcesVector(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord);
    void giveStandardInternalForcesVector(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord);
    void giveMicromorphicInternalForcesVector(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord);
    void giveLagrangianMultipliersInternalForcesVector(FloatArray &answer, TimeStep *tStep, int useUpdatedGpRecord);



    void computeForceLoadVector(FloatArray &answer, TimeStep *tStep, ValueModeType mode);
    void computeLocForceLoadVector(FloatArray &answer, TimeStep *tStep, ValueModeType mode);

    virtual IntArray &giveDisplacementDofsOrdering() {return displacementDofsOrdering;}
    virtual IntArray &giveMicromorphicDofsOrdering() {return micromorphicDofsOrdering;}
    virtual IntArray &giveLagrangianMulripliersDofsOrdering() {return micromorphicDofsOrdering;}
    // void giveLocationArrayOfDofIDs( IntArray &answer, const UnknownNumberingScheme &s, const IntArray &dofIdArray );
    void giveLocationArrayOfDofIDs(IntArray &locationArray_u, IntArray &locationArray_m, IntArray &locationArray_l, const UnknownNumberingScheme &s, const IntArray &dofIdArray_u,const IntArray &dofIdArray_m, const IntArray &dofIdArray_l );
    virtual void postInitialize();
    virtual void updateInternalState(TimeStep *tStep);




};
} // end namespace oofem

#endif
