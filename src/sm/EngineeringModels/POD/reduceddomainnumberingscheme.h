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

#ifndef reduceddomainnumberingscheme_h
#define reduceddomainnumberingscheme_h

#include "unknownnumberingscheme.h"
#include "dof.h"
#include "domain.h"
#include "dofmanager.h"

namespace oofem {
/**
 * Numbering scheme that takes into account only list of selected nodes
 *
 * @author Martin Horak
 */
class ReducedDomainNumberingScheme : public UnknownNumberingScheme
{
protected:
    Domain *domain;
    /// Container storing particular equation numbers for each node
    IntArray nodalEquationNumbers;
    /// Selected nodes
    IntArray selectedNodes;
    /// Last given number of equation
    int neq;
    /// Last given number of prescribed equation
    int pres_neq;
    /// Flag controlling wether the numbering has been initialized or not
    bool isInitialized;
    /// map form dofid to equation number
    std::map<int, std::map<int,int>> *equationMap;


public:
    /// Constructor
    ReducedDomainNumberingScheme(IntArray selectedNodes);
    /// Destructor
    virtual ~ReducedDomainNumberingScheme();

    /**
     * Initializes the receiver
     */
    virtual void init(Domain *domain, TimeStep *tStep);
    virtual bool isDefault() const { return true; }
    virtual int giveDofEquationNumber(Dof *dof) const;
    virtual int giveRequiredNumberOfDomainEquation() const;

    /// Returns total number of equations
    virtual int giveTotalNumberOfEquations() const;

    /// Returns total number of prescribed equations
    virtual int giveTotalNumberOfPrescribedEquations() const;

    /// Resets the numbering in order to start numbering again from 1
    virtual void reset();
};

} // end namespace oofem
#endif // reduceddomainnumberingschemes_h
