/*
 * File:   enrichmentitem.h
 * Author: chamrova
 *
 * Created on October 26, 2008, 11:42 AM
 */

#ifndef _ENRICHMENTITEM_H
#define	_ENRICHMENTITEM_H

#include "femcmpnn.h"
#include "domain.h"
#include "geometry.h"
#include "flotmtrx.h"
#include "enrichmentfunction.h"

class Geometry;


class EnrichmentItem : public FEMComponent
{
    public:
         /// Constructor
         EnrichmentItem (int n,Domain* aDomain);
         /// initializes EnrichmentItem from InputRecord
         IRResultType initializeFrom (InputRecord* ir);
         /// returns class name
         const char* giveClassName () const { return "EnrichmentItem" ; }
         /// Accessor
         BasicGeometry* giveGeometry();
         /// Checks whether EnrichmentItem interacts element
         bool interacts(Element* element);
         /// Computes intersection points with Element
         void computeIntersectionPoints(AList<FloatArray>* intersectionPoints, Element *element);
         /// Computes number intersection points with Element
         double computeNumberOfIntersectionPoints(Element *element);
         /// Accessor
         EnrichmentFunction* giveEnrichmentFunction() {return ef;}
         /// Sets EnrichmentFunction
         void setEnrichmentFunction(EnrichmentFunction *ef);
         /// Gives number of dofs
         int giveNumberOfDofs() { return ef->giveNumberOfDofs(); }
         /** Sets DofId Array of an Enrichment Item */
         void setDofIdArray(IntArray & dofId) { this->dofsId = dofId; }
         /// Accessor
         IntArray* getDofIdArray(){ return &dofsId; }
         /// Finds out whether a DofManager is enriched
         bool isDofManEnriched(int nodeNumber);

    protected:
         /// Geometry associated with EnrichmentItem
         BasicGeometry* geometry;
         /// EnrichmentFunction associated with the EnrichmentItem
         EnrichmentFunction *ef;
         /// Additional dofIds from Enrichment
         IntArray dofsId;
};

/** Concrete representation of EnrichmentItem */
class CrackTip : public EnrichmentItem {
    public:
        CrackTip (int n,Domain* aDomain) : EnrichmentItem(n, aDomain) {}
};

/** Concrete representation of EnrichmentItem */
class CrackInterior : public EnrichmentItem {
    public:
        CrackInterior (int n,Domain* aDomain) : EnrichmentItem(n, aDomain) {}

};
#endif	/* _ENRICHMENTITEM_H */