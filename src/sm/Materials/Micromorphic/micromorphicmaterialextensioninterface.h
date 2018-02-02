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

#ifndef micromorphicmaterialextensioninterface_h
#define micromorphicmaterialextensioninterface_h

#include "interface.h"
#include "matresponsemode.h"
#include "domain.h"

///@name micromorphicmaterialextensioninterface
//@{
#define _IFT_MicromorphicMaterialExtensionInterface_Hk "hk"
#define _IFT_MicromorphicMaterialExtensionInterface_Ak "ak"
//@}

namespace oofem {
class FloatMatrix;
class FloatArray;
class GaussPoint;
class TimeStep;



/**
 * Material interface for gradient material models.
 */
class MicromorphicMaterialExtensionInterface : public Interface
{
protected:
    Domain *dom;

    /**
     * Micromorphic parameter
     * Choosen as a penalty parameter in several cases
     */
    double Hk;


    /**
     * Micromorphic paramater related to the square of the internal length
     */
    double Ak;

public:
    /**
     * Constructor. Creates material with given number, belonging to given domain.
     * @param d Domain to which new material will belong.
     */
    MicromorphicMaterialExtensionInterface(Domain *d);
    /// Destructor.
    virtual ~MicromorphicMaterialExtensionInterface() { }

    /// Left upper block
    virtual void giveMicromorphicMatrix_dSigdUgrad(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) = 0;

    virtual void giveMicromorphicMatrix_dSigdPhi(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) = 0;

    virtual void giveMicromorphicMatrix_dSdUgrad(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) = 0;


    virtual void giveMicromorphicMatrix_dSdPhi(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) = 0;

    virtual void giveMicromorphicMatrix_dMdPhiGrad(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *tStep) = 0;


    /// micromorhpic stresses
    virtual void giveGeneralizedStressVectors (FloatArray &sigma, FloatArray &s, FloatArray &M, GaussPoint *gp, const FloatArray &displacementGradient, const FloatArray &micromorphicVar, const FloatArray micromorphicVarGrad, TimeStep *tStep) = 0;
    virtual void giveFiniteStrainGeneralizedStressVectors (FloatArray &sigma, FloatArray &s, FloatArray &M, GaussPoint *gp, const FloatArray &displacementGradient, const FloatArray &micromorphicVar, const FloatArray micromorphicVarGrad, TimeStep *tStep) = 0;

    
    virtual bool isStressTensorSymmetric(){return false;}


  
};

}
#endif
