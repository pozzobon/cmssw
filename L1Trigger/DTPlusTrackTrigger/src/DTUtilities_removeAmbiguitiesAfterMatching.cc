/*! \class DTUtilities
 *  \author Pierluigi Zotto
 *  \author Nicola Pozzobon
 *  \brief Utilities of L1 DT + Track Trigger for the HL-LHC
 *  \date 2014, Jul 9
 */

#include "L1Trigger/DTPlusTrackTrigger/interface/DTUtilities.h"

/// Method to get rid of ambiguities (setting a rejection flag)
void DTUtilities::removeAmbiguitiesAfterMatching()
{
  /// Combine together the two vector of pointers
  /// Then you can change the objects as you are working with pointers
  std::vector< DTMatch* > tempVector = theDTMatchContainer->at(1);
  tempVector.insert( tempVector.end(),
                     theDTMatchContainer->at(2).begin(), theDTMatchContainer->at(2).end() );

  /// Check the size is correct
  unsigned int sumSizes = theDTMatchContainer->at(1).size() + theDTMatchContainer->at(2).size();
  unsigned int mergedSize = tempVector.size();

  if ( sumSizes != mergedSize )
  {
    exit(0);
  }

  /// Find DTMatches with same L1 Track attached
  for ( unsigned int iDt = 0; iDt < mergedSize; iDt++ )
  {
    if ( tempVector.at(iDt)->getRejectionFlag() )
    {
      continue;
    }

    if ( tempVector.at(iDt)->getPtMatchedTrackPtr().isNull() )
    {
      continue;
    }

    /// Record the pointer to the TTTrack
    const TTTrack< Ref_PixelDigi_ >* tkPtr1 = tempVector.at(iDt)->getPtMatchedTrackPtr().get(); /// Exploit uniqueness of pointers from edm::Ptr 

    for ( unsigned int jDt = (iDt + 1); jDt < mergedSize; jDt++ )
    {
      if ( tempVector.at(jDt)->getRejectionFlag() )
      {
        continue;
      }

      if ( tempVector.at(jDt)->getPtMatchedTrackPtr().isNull() )
      {
        continue;
      }

      /// Get the pointer to the TTTrack
      const TTTrack< Ref_PixelDigi_ >* tkPtr2 = tempVector.at(jDt)->getPtMatchedTrackPtr().get(); /// Exploit uniqueness of pointers from edm::Ptr 

      if ( tkPtr1 == tkPtr2 )
      {
        DTMatch* dtMatch1 = tempVector.at(iDt);
        DTMatch* dtMatch2 = tempVector.at(jDt);

#ifdef npDEBUG
       std::cerr << "Different DTMatches matched to same track!" << std::endl;
#endif

        if ( dtMatch1->getDTCode() <= dtMatch2->getDTCode() )
        {
          dtMatch1->setRejectionFlag(true);
        }
        else
        { 
          dtMatch2->setRejectionFlag(true);
        }
      } /// End "matched to the very same track
    }
  } /// End L II track rejection

  return;
}

