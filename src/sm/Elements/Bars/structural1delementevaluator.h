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

#ifndef structural1delementevaluator_h
#define structural1delementevaluator_h

#include "../sm/Elements/structuralelementevaluator.h"

namespace oofem {
class GaussPoint;
class FloatMatrix;
class FloatArray;
class IntArray;

/**
 * Base class for 1D elements.
 *
 * 
 */
class Structural1DElementEvaluator :  public StructuralElementEvaluator
{

protected:
  /** 
   * To facilitate the transformation of 2d elements into 3d, the complexity of transformation from 3d to
   *  local 2d system can be efficiently hidden in custom FEICellGeometry wrapper, that performs the transformation
   * into local system. This way, the existing 2d intrpolation classes can be used.
   * The element maintain its FEICellGeometry, which is accesible through the giveCellGeometryWrapper. 
   * Generalization to 3d then would require only substitution of the geometry warpper and definition of 
   * element transformation matrix.
   */
  FEICellGeometry* cellGeometryWrapper;

public:
    /**
     * Constructor. Creates element with given number, belonging to given domain.
     * @param n Element number.
     * @param d Domain to which new material will belong.
     */
    Structural1DElementEvaluator() : StructuralElementEvaluator() { }
    /// Destructor.
    virtual ~Structural1DElementEvaluator(){ }
  
    
protected:
/**
 * Assemble interpolation matrix at given IP
 * In case of IGAElements, N is assumed to contain only nonzero interpolation functions
 */
    virtual void computeNMatrixAt(FloatMatrix &answer, GaussPoint *gp);
    /**
     * Assembles the strain-displacement matrix of the receiver at given integration point
     * In case of IGAElements, B is assumed to contain only contribution from nonzero interpolation functions
     */
    virtual void computeBMatrixAt(FloatMatrix &answer, GaussPoint *gp);
    virtual double computeVolumeAround(GaussPoint *gp);
    virtual void computeStressVector(FloatArray &answer, const FloatArray &strain, GaussPoint *gp, TimeStep *tStep);
    virtual void computeConstitutiveMatrixAt(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep);
    void giveDofManDofIDMask(int inode, IntArray &answer) const {
      answer = {D_u};
    }
    
};




} // end namespace oofem
#endif // structural1delement_h
