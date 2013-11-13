/*! \class   TTTrackAlgorithm_trackletLB
 *  \brief   Tracklet-based algorithm for the LB layout.
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

#ifndef L1_TRACK_TRIGGER_TRACK_ALGO_TRACKLETLB_H
#define L1_TRACK_TRIGGER_TRACK_ALGO_TRACKLETLB_H

#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/ModuleFactory.h"
#include "FWCore/Framework/interface/ESProducer.h"

#include "L1Trigger/TrackTrigger/interface/TTTrackAlgorithm.h"
#include "L1Trigger/TrackTrigger/interface/TTTrackAlgorithmRecord.h"
#include "L1Trigger/TrackTrigger/interface/TKSector.h"
#include "L1Trigger/TrackTrigger/interface/TKBarrel.h"

#include "DataFormats/GeometryCommonDetAlgo/interface/MeasurementPoint.h"
#include "Geometry/CommonTopologies/interface/Topology.h"
#include "Geometry/TrackerGeometryBuilder/interface/PixelGeomDetUnit.h"
#include "Geometry/CommonTopologies/interface/PixelTopology.h"

#include <boost/shared_ptr.hpp>
#include <memory>
#include <string>
#include <map>

template< typename T >
class TTTrackAlgorithm_trackletLB : public TTTrackAlgorithm< T >
{
  private :
    /// Data members
    double       mMagneticField;
    unsigned int nSectors;
    unsigned int nWedges;

    std::vector< std::vector< double > > tableRPhi;
    std::vector< std::vector< double > > tableZ;

  public:
    /// Constructors
    TTTrackAlgorithm_trackletLB( const StackedTrackerGeometry *aStackedGeom,
                                 double aMagneticField,
                                 unsigned int aSectors,
                                 unsigned int aWedges,
                                 std::vector< std::vector< double > > aTableRPhi,
                                 std::vector< std::vector< double > > aTableZ )
      : TTTrackAlgorithm< T > ( aStackedGeom, __func__ )
    {
      mMagneticField = aMagneticField;
      nSectors = aSectors;
      nWedges = aWedges;

      tableRPhi = aTableRPhi;
      tableZ = aTableZ;
    }

    /// Destructor
    ~TTTrackAlgorithm_trackletLB(){}

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
void TTTrackAlgorithm_trackletLB< T >::FillSectors( std::map< std::pair< unsigned int, unsigned int >, TKSector< T > > *outputSectorMap,
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
    } /// End of already existing Sector case
  } /// End of loop over input Stubs
}

/// Find the Seeds
template< typename T >
void TTTrackAlgorithm_trackletLB< T >::FindSeeds( std::vector< TTTrack< T > > *output,
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
    /// Barrel containers
    TKSector< T > thisSector0 = mapIter->second;
    std::map< unsigned int, TKBarrel< T > > mapBarrelLayer0 = thisSector0.getBarrelLayerMap();

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
            mapBarrelLayer0[mapBarrelIter->first].findSeeds( mapBarrelLayer1[mapBarrelIter->first + 1], mapIter->first.first, mapIter->first.second, true );
          }

          std::vector< TTTrack< T > > tempTracks = mapBarrelLayer0[mapBarrelIter->first].getTracks();
          output->insert( output->end(), tempTracks.begin(), tempTracks.end() );

        } /// End of loop over the barrel map
      }
    } /// End of loop over the neighbors
  } /// End of loop over the input sector map
}

/// Propagate the Seeds
template< typename T >
void TTTrackAlgorithm_trackletLB< T >::FindMatches( std::vector< TTTrack< T > > *output,
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

        /// Loop over Barrel-Barrel Seeds
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
            /// NOTE: differently from the BE case, here the tables are
            /// encoded by super-layer, SL 0 does not exist, SL 1 is Barrels 1 + 2
            /// etc, which means that SL = (L+1)/2 (integer division)
            mapBarrelLayer0[mapBarrelIter->first].findMatches( mapBarrelLayer1[anotherMapBarrelIter->first],
                                                               tableRPhi.at((mapBarrelIter->first+1)/2).at((anotherMapBarrelIter->first+1)/2),
                                                               tableZ.at((mapBarrelIter->first+1)/2).at((anotherMapBarrelIter->first+1)/2) );
          }

          std::vector< TTTrack< T > > tempTracks = mapBarrelLayer0[mapBarrelIter->first].getTracks();
          output->insert( output->end(), tempTracks.begin(), tempTracks.end() );

        } /// End of loop over the Barrel map
      }
    } /// End of loop over the neighbors
  } /// End of loop over the input sector map
}






/*! \class   ES_TTTrackAlgorithm_trackletLB
 *  \brief   Class to declare the algorithm to the framework
 *
 *  \author Nicola Pozzobon
 *  \date   2013, Jul 18
 *
 */

template< typename T >
class ES_TTTrackAlgorithm_trackletLB : public edm::ESProducer
{
  private:
    /// Data members
    boost::shared_ptr< TTTrackAlgorithm< T > > _theAlgo;

    /// Number of Sectors
    unsigned int  mSectors;
    unsigned int  mWedges;

    /// projection windows
    std::vector< std::vector< double > > setRhoPhiWin;
    std::vector< std::vector< double > > setZWin;

  public:
    /// Constructor
    ES_TTTrackAlgorithm_trackletLB( const edm::ParameterSet & p )
      : mSectors( p.getParameter< int >("NumSectors") ), mWedges( p.getParameter< int >("NumWedges") )
    {
      std::vector< edm::ParameterSet > vPSet = p.getParameter< std::vector< edm::ParameterSet > >("ProjectionWindows");
      std::vector< edm::ParameterSet >::const_iterator iPSet;
      for ( iPSet = vPSet.begin(); iPSet != vPSet.end(); iPSet++ )
      {
        setRhoPhiWin.push_back( iPSet->getParameter< std::vector< double > >("RhoPhiWin") );
        setZWin.push_back( iPSet->getParameter< std::vector< double > >("ZWin") );
      }

      setWhatProduced( this );
    }

    /// Destructor
    virtual ~ES_TTTrackAlgorithm_trackletLB() {}

    /// Implement the producer
    boost::shared_ptr< TTTrackAlgorithm< T > > produce( const TTTrackAlgorithmRecord & record )
    {
      /// Get magnetic field
      edm::ESHandle< MagneticField > magnet;
      record.getRecord< IdealMagneticFieldRecord >().get(magnet);
      double mMagneticFieldStrength = magnet->inTesla(GlobalPoint(0,0,0)).z();
      double mMagneticFieldRounded = (floor(mMagneticFieldStrength*10.0 + 0.5))/10.0;

      edm::ESHandle< StackedTrackerGeometry > StackedTrackerGeomHandle;
      record.getRecord< StackedTrackerGeometryRecord >().get( StackedTrackerGeomHandle );

      TTTrackAlgorithm< T >* TTTrackAlgo =
        new TTTrackAlgorithm_trackletLB< T >( &(*StackedTrackerGeomHandle),
                                              mMagneticFieldRounded,
                                              mSectors,
                                              mWedges,
                                              setRhoPhiWin,
                                              setZWin );

      _theAlgo = boost::shared_ptr< TTTrackAlgorithm< T > >( TTTrackAlgo );
      return _theAlgo;
    }

};

#endif

