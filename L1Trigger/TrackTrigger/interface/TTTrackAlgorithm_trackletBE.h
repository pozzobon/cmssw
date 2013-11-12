/*! \class   TTTrackAlgorithm_trackletBE
 *  \brief   Tracklet-based algorithm for the BE layout.
 *  \details After moving from SimDataFormats to DataFormats,
 *           the template structure of the class was maintained
 *           in order to accomodate any types other than PixelDigis
 *           in case there is such a need in the future.
 *
 *  \author Anders Ryd
 *  \author Emmanuele Salvati
 *  \author Nicola Pozzobon
 *  \date   2013, Jul 18
 *
 */

#ifndef L1_TRACK_TRIGGER_TRACK_ALGO_TRACKLETBE_H
#define L1_TRACK_TRIGGER_TRACK_ALGO_TRACKLETBE_H

#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/ModuleFactory.h"
#include "FWCore/Framework/interface/ESProducer.h"

#include "L1Trigger/TrackTrigger/interface/TTTrackAlgorithm.h"
#include "L1Trigger/TrackTrigger/interface/TTTrackAlgorithmRecord.h"
#include "L1Trigger/TrackTrigger/interface/TKSector.h"
#include "L1Trigger/TrackTrigger/interface/TKBarrel.h"
#include "L1Trigger/TrackTrigger/interface/TKEndcap.h"

#include "DataFormats/GeometryCommonDetAlgo/interface/MeasurementPoint.h"
#include "Geometry/CommonTopologies/interface/Topology.h"
#include "Geometry/TrackerGeometryBuilder/interface/PixelGeomDetUnit.h"
#include "Geometry/CommonTopologies/interface/PixelTopology.h"

#include <boost/shared_ptr.hpp>
#include <memory>
#include <string>
#include <map>

template< typename T >
class TTTrackAlgorithm_trackletBE : public TTTrackAlgorithm< T >
{
  private :
    /// Data members
    double       mMagneticField;
    unsigned int nSectors;
    unsigned int nWedges;

    std::vector< std::vector< double > > tableRPhiBB;
    std::vector< std::vector< double > > tableZBB;
    std::vector< std::vector< double > > tableRPhiBE;
    std::vector< std::vector< double > > tableZBE;
    std::vector< std::vector< double > > tableRPhiEB;
    std::vector< std::vector< double > > tableZEB;
    std::vector< std::vector< double > > tableRPhiEE;
    std::vector< std::vector< double > > tableZEE;

    std::vector< std::vector< double > > tableRPhiBE_PS;
    std::vector< std::vector< double > > tableZBE_PS;
    std::vector< std::vector< double > > tableRPhiEB_PS;
    std::vector< std::vector< double > > tableZEB_PS;
    std::vector< std::vector< double > > tableRPhiEE_PS;
    std::vector< std::vector< double > > tableZEE_PS;

  public:
    /// Constructors
    TTTrackAlgorithm_trackletBE( const StackedTrackerGeometry *aStackedGeom,
                                 double aMagneticField,
                                 unsigned int aSectors,
                                 unsigned int aWedges,
                                 std::vector< std::vector< double > > aTableRPhiBB,
                                 std::vector< std::vector< double > > aTableZBB,
                                 std::vector< std::vector< double > > aTableRPhiBE,
                                 std::vector< std::vector< double > > aTableZBE,
                                 std::vector< std::vector< double > > aTableRPhiBE_PS, 
                                 std::vector< std::vector< double > > aTableZBE_PS, 
                                 std::vector< std::vector< double > > aTableRPhiEB,
                                 std::vector< std::vector< double > > aTableZEB,
                                 std::vector< std::vector< double > > aTableRPhiEB_PS,
                                 std::vector< std::vector< double > > aTableZEB_PS,
                                 std::vector< std::vector< double > > aTableRPhiEE,
                                 std::vector< std::vector< double > > aTableZEE,
                                 std::vector< std::vector< double > > aTableRPhiEE_PS,
                                 std::vector< std::vector< double > > aTableZEE_PS )
      : TTTrackAlgorithm< T > ( aStackedGeom, __func__ )
    {
      mMagneticField = aMagneticField;
      nSectors = aSectors;
      nWedges = aWedges;

      tableRPhiBB = aTableRPhiBB;
      tableZBB = aTableZBB;
      tableRPhiBE = aTableRPhiBE;
      tableZBE = aTableZBE;
      tableRPhiEB = aTableRPhiEB;
      tableZEB = aTableZEB;
      tableRPhiEE = aTableRPhiEE;
      tableZEE = aTableZEE;

      tableRPhiBE_PS = aTableRPhiBE_PS;
      tableZBE_PS = aTableZBE_PS;
      tableRPhiEB_PS = aTableRPhiEB_PS;
      tableZEB_PS = aTableZEB_PS;
      tableRPhiEE_PS = aTableRPhiEE_PS;
      tableZEE_PS = aTableZEE_PS;
    }

    /// Destructor
    ~TTTrackAlgorithm_trackletBE(){}

    /// Fill sectors with stubs and structures
    void FillSectors( std::map< std::pair< unsigned int, unsigned int >, TKSector< T > > *outputSectorMap,
                      edm::Handle< std::vector< TTStub< T > > > &input ) const;

    /// Find the Seeds
    void FindSeeds( std::vector< TTTrack< T > > *output,
                    std::map< std::pair< unsigned int, unsigned int >, TKSector< T > > *inputSectorMap ) const;

    /// Propagate the Seeds
    void FindMatches( std::vector< TTTrack< T > > *output,
                      std::map< std::pair< unsigned int, unsigned int >, TKSector< T > > *inputSectorMap ) const;

    /// Return the number of Sectors
    unsigned int ReturnNumberOfSectors() const { return nSectors; } /// Phi
    unsigned int ReturnNumberOfWedges() const  { return nWedges; } /// Eta

    /// Return the value of the magnetic field
    double ReturnMagneticField() const { return mMagneticField; }

    /// Fit the Track
    /// Take it from the parent class

}; /// Close class

/*! \brief   Implementation of methods
 *  \details Here, in the header file, the methods which do not depend
 *           on the specific type <T> that can fit the template.
 *           Other methods, with type-specific features, are implemented
 *           in the source file.
 */

/// Fill sectors with stubs
template < typename T >
void TTTrackAlgorithm_trackletBE< T >::FillSectors( std::map< std::pair< unsigned int, unsigned int >, TKSector< T > > *outputSectorMap,
                                                    edm::Handle< std::vector< TTStub< T > > > &input ) const
{
  /// Prepare output
  outputSectorMap->clear();

  /// Loop over input stubs
  typename std::vector< TTStub< T > >::const_iterator inputIter;
  unsigned int j = 0; /// Counter needed to build the edm::Ptr to the TTStub
  for ( inputIter = input->begin();
        inputIter != input->end();
        ++inputIter )
  {
    /// Make the pointer to be put in the map and, later on, in the Track
    edm::Ptr< TTStub< T > > tempStubPtr( input, j++ );

    /// Calculate Sector
    /// From 0 to nSectors-1
    /// Sector 0 centered on Phi = 0 and symmetric around it
    GlobalPoint stubPosition = TTTrackAlgorithm< T >::theStackedTracker->findGlobalPosition( tempStubPtr.get() );
    double stubPhi = stubPosition.phi();
    stubPhi += M_PI/nSectors;
    if ( stubPhi < 0 )
    {
      stubPhi += 2*M_PI;
    }
    unsigned int thisSector = floor( 0.5*stubPhi*nSectors/M_PI );

    /// Calculate Wedge
    /// From 0 to nWedges-1
    double stubEta = stubPosition.eta();
    stubEta += 2.5; /// Bring eta = -2.5 to 0
    if ( stubEta < 0.0 || stubEta > 5.0 )
    {
      /// Accept only stubs within -2.5, 2.5 range
      continue;
    }
    unsigned int thisWedge = floor( stubEta*nWedges/5.0 );

    /// Build the key to the map (by Sector / Wedge)
    std::pair< unsigned int, unsigned int > mapKey = std::make_pair( thisSector, thisWedge );

    /// Get the DetId of the Pt-module of this stub
    StackedTrackerDetId detIdStub( inputIter->getDetId() );

    /// Check if this stub is already present in some sector
    /// in the output map
    typename std::map< std::pair< unsigned int, unsigned int >, TKSector< T > >::iterator sectorMapIter;
    sectorMapIter = outputSectorMap->find( mapKey );

    if ( sectorMapIter == outputSectorMap->end() )
    {
      /// This stub belongs to a new sector, not stored yet
      /// Create the sector object
      TKSector< T > newSector( TTTrackAlgorithm< T >::theStackedTracker );

      /// Different actions for Barrel and Endcap Stubs
      if ( detIdStub.isBarrel() )
      {
        unsigned int iLayer = detIdStub.iLayer();

        /// First time the Sector is filled
        /// Create the Barrel container from scratch
        TKBarrel< T > newBarrel( TTTrackAlgorithm< T >::theStackedTracker, mMagneticField );
        newBarrel.addStubPtr( tempStubPtr );
        newSector.addBarrel( iLayer, newBarrel );
      }
      else if ( detIdStub.isEndcap() )
      {
        unsigned int iSide = detIdStub.iSide();
        unsigned int iDisk = detIdStub.iDisk();

        /// First time the Sector is filled
        /// Create the Endcap container from scratch
        TKEndcap< T > newEndcap( TTTrackAlgorithm< T >::theStackedTracker, mMagneticField );
        newEndcap.addStubPtr( tempStubPtr );
        newSector.addEndcap( std::make_pair( iSide, iDisk ), newEndcap );
      }

      outputSectorMap->insert( std::make_pair( mapKey, newSector ) );

    } /// End of new Sector case
    else
    {
      /// The Sector already exists
      /// so the Stub must be properly attached to it
      /// Different actions for Barrel and Endcap Stubs
      if ( detIdStub.isBarrel() )
      {
        unsigned int iLayer = detIdStub.iLayer();

        /// Check if the Barrel Layer already exists
        typename std::map< unsigned int, TKBarrel< T > >::iterator barrelMapIter;
        std::map< unsigned int, TKBarrel< T > > tempBarrelMap = (*outputSectorMap)[mapKey].getBarrelLayerMap();
        barrelMapIter = tempBarrelMap.find( iLayer );

        if ( barrelMapIter == tempBarrelMap.end() )
        {
          /// New Barrel
          TKBarrel< T > aNewBarrel( TTTrackAlgorithm< T >::theStackedTracker, mMagneticField );
          aNewBarrel.addStubPtr( tempStubPtr );
          (*outputSectorMap)[mapKey].addBarrel( iLayer, aNewBarrel );
        }
        else
        {
          /// Already existing Barrel
          (*outputSectorMap)[mapKey].addStubPtrToBarrel( iLayer, tempStubPtr );
        }
      }
      else if ( detIdStub.isEndcap() )
      {
        unsigned int iSide = detIdStub.iSide();
        unsigned int iDisk = detIdStub.iDisk();
        std::pair< unsigned int, unsigned int > tempKey = std::make_pair( iSide, iDisk );

        /// Check if the Endcap Disk already exists
        typename std::map< std::pair< unsigned int, unsigned int >, TKEndcap< T > >::iterator endcapMapIter;
        std::map< std::pair< unsigned int, unsigned int >, TKEndcap< T > > tempEndcapMap = (*outputSectorMap)[mapKey].getEndcapDiskMap();
        endcapMapIter = tempEndcapMap.find( tempKey );

        if ( endcapMapIter == tempEndcapMap.end() )
        {
          /// New Endcap
          TKEndcap< T > aNewEndcap( TTTrackAlgorithm< T >::theStackedTracker, mMagneticField );
          aNewEndcap.addStubPtr( tempStubPtr );
          (*outputSectorMap)[mapKey].addEndcap( tempKey, aNewEndcap );
        }
        else
        {
          /// Already existing Endcap
          (*outputSectorMap)[mapKey].addStubPtrToEndcap( tempKey, tempStubPtr );
        }
      }
    } /// End of already existing Sector case
  } /// End of loop over input Stubs
}

/// Find the Seeds
template< typename T >
void TTTrackAlgorithm_trackletBE< T >::FindSeeds( std::vector< TTTrack< T > > *output,
                                                  std::map< std::pair< unsigned int, unsigned int >, TKSector< T > > *inputSectorMap ) const
{
  /// Prepare output
  output->clear();

  /// Loop over the input map
  typename std::map< std::pair< unsigned int, unsigned int >, TKSector< T > >::iterator mapIter;
  typename std::map< std::pair< unsigned int, unsigned int >, TKSector< T > >::iterator anotherMapIter;

  for ( mapIter = inputSectorMap->begin();
        mapIter != inputSectorMap->end();
        mapIter++ )
  {
    /// Get the current Sector container and the corresponding
    /// Barrel and Endcap containers
    TKSector< T > thisSector0 = mapIter->second;
    std::map< unsigned int, TKBarrel< T > > mapBarrelLayer0 = thisSector0.getBarrelLayerMap();
    std::map< std::pair< unsigned int, unsigned int >, TKEndcap< T > > mapEndcapDisk0 = thisSector0.getEndcapDiskMap();

    /// Get the Sector/Wedge of the present list
    unsigned int curSector0 = mapIter->first.first + this->ReturnNumberOfSectors(); /// This is to use the %nSectors later
    unsigned int curWedge0 = mapIter->first.second;

    /// Loop over the sector and its two neighbors
    for ( unsigned int iSector = 0; iSector < 2; iSector++ )
    {
      for ( unsigned int iWedge = 0; iWedge < 2; iWedge++)
      {
        /// Find the correct sector index
        unsigned int curSector = ( curSector0 + iSector -1 )%(this->ReturnNumberOfSectors());
        int curWedge = curWedge0 + iWedge - 1;
        if ( curWedge < 0 || curWedge >= (int)(this->ReturnNumberOfWedges()) )
          continue;

        /// Now we are in the correct Sector/Wedge pair
        anotherMapIter = inputSectorMap->find( std::make_pair( curSector, curWedge ) );

        /// Skip if the Sector does not exist or is empty
        if ( anotherMapIter == inputSectorMap->end() )
          continue;

        /// Get the Sector of the nested loop, to get the second stub
        /// in order to make Tracklets
        TKSector< T > thisSector1 = anotherMapIter->second;
        std::map< unsigned int, TKBarrel< T > > mapBarrelLayer1 = thisSector1.getBarrelLayerMap();
        std::map< std::pair< unsigned int, unsigned int >, TKEndcap< T > > mapEndcapDisk1 = thisSector1.getEndcapDiskMap();

        /// Here we have all the containers of Stubs in this Sector
        /// already mapped by Layer and Disk ...

        /// xxxx0 is meant to store the inner Stub of the candidate Tracklet
        /// xxxx1 is meant to store the outer Stub of the candidate Tracklet
        /// This is safe as for each Sector we are also
        /// looking at all its neighbors

        /// Loop over all the Barrel containers for inner Stubs
        typename std::map< unsigned int, TKBarrel< T > >::iterator mapBarrelIter;
        for ( mapBarrelIter = mapBarrelLayer0.begin();
              mapBarrelIter != mapBarrelLayer0.end();
              ++mapBarrelIter )
        {
          /// If there is the corresponding next Layer in the Sector,
          /// then you can look for "traditional" Tracklets
          if ( mapBarrelLayer1.find( mapBarrelIter->first + 1 ) != mapBarrelLayer1.end() )
          {
            mapBarrelLayer0[mapBarrelIter->first].findSeeds( mapBarrelLayer1[mapBarrelIter->first + 1], mapIter->first.first, mapIter->first.second );
          }

          /// If there is also the first Disk in the Sector,
          /// then you can look for "mixed" Tracklets too
          for ( unsigned int aSide = 1; aSide <= 2; aSide++ )
          {
            if ( mapEndcapDisk1.find( std::make_pair( aSide, 1 ) ) != mapEndcapDisk1.end() )
            {
              mapBarrelLayer0[mapBarrelIter->first].findSeeds( mapEndcapDisk1[std::make_pair( aSide, 1 )], mapIter->first.first, mapIter->first.second );
            }
          }

          std::vector< TTTrack< T > > tempTracks = mapBarrelLayer0[mapBarrelIter->first].getTracks();
          output->insert( output->end(), tempTracks.begin(), tempTracks.end() );

        } /// End of loop over the barrel map

        /// Loop over all the Endcap containers for inner Stubs
        typename std::map< std::pair< unsigned int, unsigned int>, TKEndcap< T > >::iterator mapEndcapIter;
        for ( mapEndcapIter = mapEndcapDisk0.begin();
              mapEndcapIter != mapEndcapDisk0.end();
              ++mapEndcapIter )
        {
          /// If there is the corresponding next Disk in the Sector,
          /// then you can look for "traditional" Tracklets
          if ( mapEndcapDisk1.find( std::make_pair( mapEndcapIter->first.first, mapEndcapIter->first.second + 1 ) ) != mapEndcapDisk1.end() )
          {
              mapEndcapDisk0[mapEndcapIter->first].findSeeds( mapEndcapDisk1[std::make_pair( mapEndcapIter->first.first, mapEndcapIter->first.second + 1 )],
                                                              mapIter->first.first, mapIter->first.second );
          }

          std::vector< TTTrack< T > > tempTracks = mapEndcapDisk0[mapEndcapIter->first].getTracks();
          output->insert( output->end(), tempTracks.begin(), tempTracks.end() );

        } /// End of loop over the endcap map
      }
    } /// End of loop over the neighbors
  } /// End of loop over the input sector map
}

/// Propagate the Seeds
template< typename T >
void TTTrackAlgorithm_trackletBE< T >::FindMatches( std::vector< TTTrack< T > > *output,
                                                    std::map< std::pair< unsigned int, unsigned int >, TKSector< T > > *inputSectorMap ) const
{
  /// Prepare output
  output->clear();

  /// Loop over the input map
  typename std::map< std::pair< unsigned int, unsigned int >, TKSector< T > >::iterator mapIter;
  typename std::map< std::pair< unsigned int, unsigned int >, TKSector< T > >::iterator anotherMapIter;
  for ( mapIter = inputSectorMap->begin();
        mapIter != inputSectorMap->end();
        mapIter++ )
  {
    /// Get the current Sector container
    TKSector< T > thisSector0 = mapIter->second;
    std::map< unsigned int, TKBarrel< T > > mapBarrelLayer0 = thisSector0.getBarrelLayerMap();
    std::map< std::pair< unsigned int, unsigned int >, TKEndcap< T > > mapEndcapDisk0 = thisSector0.getEndcapDiskMap();

    /// Get the Sector/Wedge of the present list
    unsigned int curSector0 = mapIter->first.first + this->ReturnNumberOfSectors(); /// This is to use the %nSectors later
    unsigned int curWedge0 = mapIter->first.second;

    /// Loop over the Sector and its two neighbors
    for ( unsigned int iSector = 0; iSector < 2; iSector++ )
    {
      for ( unsigned int iWedge = 0; iWedge < 2; iWedge++)
      {
        /// Find the correct Sector index
        unsigned int curSector = ( curSector0 + iSector -1 )%(this->ReturnNumberOfSectors());
        int curWedge = curWedge0 + iWedge - 1;
        if ( curWedge < 0 || curWedge >= (int)(this->ReturnNumberOfWedges()) )
          continue;

        /// Now we are in the correct Sector/Wedge pair
        anotherMapIter = inputSectorMap->find( std::make_pair( curSector, curWedge ) );

        /// Skip to the next one if the Sector does not exist or is empty
        if ( anotherMapIter == inputSectorMap->end() )
          continue;

        /// Get the Sector to search for Stubs to be matched to the Seed
        TKSector< T > thisSector1 = anotherMapIter->second;
        std::map< unsigned int, TKBarrel< T > > mapBarrelLayer1 = thisSector1.getBarrelLayerMap();
        std::map< std::pair< unsigned int, unsigned int >, TKEndcap< T > > mapEndcapDisk1 = thisSector1.getEndcapDiskMap();

        /// Here we have all the containers of Stubs in this Sector
        /// already mapped by Layer and Disk ...

        /// The containers also do contain the Seeds corresponding to them,
        /// i.e., a Barrel contains the Seeds with inner Stub which is
        /// present onto that specific Barrel, etc.

        /// xxxx0 contains the Seeds to be propagated
        /// xxxx1 contains the Stubs to be attached to the Seed

        /// Now propagate Seeds towards all other containers within
        /// the Sector and its neighbors
        typename std::map< unsigned int, TKBarrel< T > >::iterator mapBarrelIter;
        typename std::map< unsigned int, TKBarrel< T > >::iterator anotherMapBarrelIter;
        typename std::map< std::pair< unsigned int, unsigned int >, TKEndcap< T > >::iterator mapEndcapIter;
        typename std::map< std::pair< unsigned int, unsigned int >, TKEndcap< T > >::iterator anotherMapEndcapIter;

        /// Loop over Barrel-Barrel Seeds and Barrel-Endcap Seeds
        for ( mapBarrelIter = mapBarrelLayer0.begin();
              mapBarrelIter != mapBarrelLayer0.end();
              ++mapBarrelIter )
        {

          /// Propagate to other Barrels
          for ( anotherMapBarrelIter = mapBarrelLayer1.begin();
                anotherMapBarrelIter != mapBarrelLayer1.end();
                ++anotherMapBarrelIter )
          {
            /// In case the candidate Stub is in the same Layer as
            /// the ones of the Seed, skip to the next one!
            if ( mapBarrelIter->first == anotherMapBarrelIter->first ||
                 mapBarrelIter->first + 1 == anotherMapBarrelIter->first )
              continue;

            /// Redundancies and additional cross checks in "findMatches"
            mapBarrelLayer0[mapBarrelIter->first].findMatches( mapBarrelLayer1[anotherMapBarrelIter->first],
                                                               tableRPhiBB.at(mapBarrelIter->first).at(anotherMapBarrelIter->first),
                                                               tableZBB.at(mapBarrelIter->first).at(anotherMapBarrelIter->first) );
          }

          /// Propagate to other Endcaps
          for ( anotherMapEndcapIter = mapEndcapDisk1.begin();
                anotherMapEndcapIter != mapEndcapDisk1.end();
                ++anotherMapEndcapIter )
          {
            /// Prepare the matching windows
            std::pair< double, double > pairRPhi = std::make_pair( tableRPhiBE_PS.at(mapBarrelIter->first).at(anotherMapEndcapIter->first.second),
                                                                   tableRPhiBE.at(mapBarrelIter->first).at(anotherMapEndcapIter->first.second) );
            std::pair< double, double > pairR = std::make_pair( tableZBE_PS.at(mapBarrelIter->first).at(anotherMapEndcapIter->first.second),
                                                                tableZBE.at(mapBarrelIter->first).at(anotherMapEndcapIter->first.second) );

            /// Redundancies and additional cross checks in "findMatches"
            mapBarrelLayer0[mapBarrelIter->first].findMatches( mapEndcapDisk1[anotherMapEndcapIter->first], pairRPhi, pairR );
          }

          std::vector< TTTrack< T > > tempTracks = mapBarrelLayer0[mapBarrelIter->first].getTracks();
          output->insert( output->end(), tempTracks.begin(), tempTracks.end() );

        } /// End of loop over the Barrel map

        /// Loop over the Endcap-Endcap Seeds
        for ( mapEndcapIter = mapEndcapDisk0.begin();
              mapEndcapIter != mapEndcapDisk0.end();
              ++mapEndcapIter )
        {
          /// Propagate to other Barrels
          for ( anotherMapBarrelIter = mapBarrelLayer1.begin();
                anotherMapBarrelIter != mapBarrelLayer1.end();
                ++anotherMapBarrelIter )
          {
            /// Prepare the matching windows
            std::pair< double, double > pairRPhi = std::make_pair( tableRPhiEB_PS.at(mapEndcapIter->first.second).at(anotherMapBarrelIter->first),
                                                                   tableRPhiEB.at(mapEndcapIter->first.second).at(anotherMapBarrelIter->first) );
            std::pair< double, double > pairZ = std::make_pair( tableZEB_PS.at(mapEndcapIter->first.second).at(anotherMapBarrelIter->first),
                                                                tableZEB.at(mapEndcapIter->first.second).at(anotherMapBarrelIter->first) );

            /// Redundancies and additional cross checks in "findMatches"
            mapEndcapDisk0[mapEndcapIter->first].findMatches( mapBarrelLayer1[anotherMapBarrelIter->first], pairRPhi, pairZ );
          }

          /// Propagate to other Endcaps
          for ( anotherMapEndcapIter = mapEndcapDisk1.begin();
                anotherMapEndcapIter != mapEndcapDisk1.end();
                ++anotherMapEndcapIter )
          {
            /// In case the candidate Stub is in the same Disks as
            /// the ones of the Seed, Skip to the next one!
            if ( mapEndcapIter->first == anotherMapEndcapIter->first ||
                 std::make_pair( mapEndcapIter->first.first, mapEndcapIter->first.second + 1 ) == anotherMapEndcapIter->first )
              continue;

            /// Prepare the matching windows
            std::pair< double, double > pairRPhi = std::make_pair( tableRPhiEE_PS.at(mapEndcapIter->first.second).at(anotherMapEndcapIter->first.second),
                                                                   tableRPhiEE.at(mapEndcapIter->first.second).at(anotherMapEndcapIter->first.second) );
            std::pair< double, double > pairR = std::make_pair( tableZEE_PS.at(mapEndcapIter->first.second).at(anotherMapEndcapIter->first.second),
                                                                tableZEE.at(mapEndcapIter->first.second).at(anotherMapEndcapIter->first.second) );

            /// Redundancies and additional cross checks in "findMatches"
            mapEndcapDisk0[mapEndcapIter->first].findMatches( mapEndcapDisk1[anotherMapEndcapIter->first], pairRPhi, pairR );
          }

          std::vector< TTTrack< T > > tempTracks = mapEndcapDisk0[mapEndcapIter->first].getTracks();
          output->insert( output->end(), tempTracks.begin(), tempTracks.end() );

        } /// End of loop over the Endcap map
      }
    } /// End of loop over the neighbors
  } /// End of loop over the input sector map
}






/*! \class   ES_TTTrackAlgorithm_trackletBE
 *  \brief   Class to declare the algorithm to the framework
 *
 *  \author Nicola Pozzobon
 *  \date   2013, Jul 18
 *
 */

template< typename T >
class ES_TTTrackAlgorithm_trackletBE : public edm::ESProducer
{
  private:
    /// Data members
    boost::shared_ptr< TTTrackAlgorithm< T > > _theAlgo;

    /// Number of Sectors
    unsigned int mSectors;
    unsigned int mWedges;

    /// Projection windows
    std::vector< std::vector< double > > setRhoPhiWinBB;
    std::vector< std::vector< double > > setZWinBB;
    std::vector< std::vector< double > > setRhoPhiWinBE;
    std::vector< std::vector< double > > setZWinBE;
    std::vector< std::vector< double > > setRhoPhiWinEB;
    std::vector< std::vector< double > > setZWinEB;
    std::vector< std::vector< double > > setRhoPhiWinEE;
    std::vector< std::vector< double > > setZWinEE;

    /// PS Modules variants
    /// NOTE these are not needed for the Barrel-Barrel case
    std::vector< std::vector< double > > setRhoPhiWinBE_PS;
    std::vector< std::vector< double > > setZWinBE_PS;
    std::vector< std::vector< double > > setRhoPhiWinEB_PS;
    std::vector< std::vector< double > > setZWinEB_PS;
    std::vector< std::vector< double > > setRhoPhiWinEE_PS;
    std::vector< std::vector< double > > setZWinEE_PS;

  public:
    /// Constructor
    ES_TTTrackAlgorithm_trackletBE( const edm::ParameterSet & p )
      : mSectors( p.getParameter< int >("NumSectors") ),
        mWedges( p.getParameter< int >("NumWedges") )
    {
      std::vector< edm::ParameterSet > vPSet;
      std::vector< edm::ParameterSet >::const_iterator iPSet;

      vPSet = p.getParameter< std::vector< edm::ParameterSet > >("ProjectionWindowsBarrelBarrel");
      for ( iPSet = vPSet.begin(); iPSet != vPSet.end(); iPSet++ )
      {
        setRhoPhiWinBB.push_back( iPSet->getParameter< std::vector< double > >("RhoPhiWin") );
        setZWinBB.push_back( iPSet->getParameter< std::vector< double > >("ZWin") );
      }

      vPSet = p.getParameter< std::vector< edm::ParameterSet > >("ProjectionWindowsBarrelEndcap");
      for ( iPSet = vPSet.begin(); iPSet != vPSet.end(); iPSet++ )
      {
        setRhoPhiWinBE.push_back( iPSet->getParameter< std::vector< double > >("RhoPhiWin") );
        setZWinBE.push_back( iPSet->getParameter< std::vector< double > >("ZWin") );
        setRhoPhiWinBE_PS.push_back( iPSet->getParameter< std::vector< double > >("RhoPhiWinPS") );
        setZWinBE_PS.push_back( iPSet->getParameter< std::vector< double > >("ZWinPS") );
      }

      vPSet = p.getParameter< std::vector< edm::ParameterSet > >("ProjectionWindowsEndcapBarrel");
      for ( iPSet = vPSet.begin(); iPSet != vPSet.end(); iPSet++ )
      {
        setRhoPhiWinEB.push_back( iPSet->getParameter< std::vector< double > >("RhoPhiWin") );
        setZWinEB.push_back( iPSet->getParameter< std::vector< double > >("ZWin") );
        setRhoPhiWinEB_PS.push_back( iPSet->getParameter< std::vector< double > >("RhoPhiWinPS") );
        setZWinEB_PS.push_back( iPSet->getParameter< std::vector< double > >("ZWinPS") );
      }

      vPSet = p.getParameter< std::vector< edm::ParameterSet > >("ProjectionWindowsEndcapEndcap");
      for ( iPSet = vPSet.begin(); iPSet != vPSet.end(); iPSet++ )
      {
        setRhoPhiWinEE.push_back( iPSet->getParameter< std::vector< double > >("RhoPhiWin") );
        setZWinEE.push_back( iPSet->getParameter< std::vector< double > >("ZWin") );
        setRhoPhiWinEE_PS.push_back( iPSet->getParameter< std::vector< double > >("RhoPhiWinPS") );
        setZWinEE_PS.push_back( iPSet->getParameter< std::vector< double > >("ZWinPS") );
      }

      setWhatProduced( this );
    }

    /// Destructor
    virtual ~ES_TTTrackAlgorithm_trackletBE() {}

    /// Implement the producer
    boost::shared_ptr< TTTrackAlgorithm< T > > produce( const TTTrackAlgorithmRecord & record )
    {
      /// Get magnetic field
      edm::ESHandle< MagneticField > magnet;
      record.getRecord< IdealMagneticFieldRecord >().get(magnet);
      double mMagneticFieldStrength = magnet->inTesla( GlobalPoint(0,0,0) ).z();
      double mMagneticFieldRounded = (floor(mMagneticFieldStrength*10.0 + 0.5))/10.0;

      edm::ESHandle< StackedTrackerGeometry > StackedTrackerGeomHandle;
      record.getRecord< StackedTrackerGeometryRecord >().get( StackedTrackerGeomHandle );

      TTTrackAlgorithm< T >* TTTrackAlgo =
        new TTTrackAlgorithm_trackletBE< T >( &(*StackedTrackerGeomHandle),
                                              mMagneticFieldRounded,
                                              mSectors,
                                              mWedges,
                                              setRhoPhiWinBB,
                                              setZWinBB,
                                              setRhoPhiWinBE,
                                              setZWinBE,
                                              setRhoPhiWinBE_PS,
                                              setZWinBE_PS,
                                              setRhoPhiWinEB,
                                              setZWinEB,
                                              setRhoPhiWinEB_PS,
                                              setZWinEB_PS,
                                              setRhoPhiWinEE,
                                              setZWinEE,
                                              setRhoPhiWinEE_PS,
                                              setZWinEE_PS );

      _theAlgo = boost::shared_ptr< TTTrackAlgorithm< T > >( TTTrackAlgo );
      return _theAlgo;
    }

};

#endif

