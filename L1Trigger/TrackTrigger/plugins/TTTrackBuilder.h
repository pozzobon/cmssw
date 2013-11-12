/*! \class   TTTrackBuilder
 *  \brief   Plugin to load the Track finding algorithm and produce the
 *           collection of Tracks that goes in the event content.
 *  \details After moving from SimDataFormats to DataFormats,
 *           the template structure of the class was maintained
 *           in order to accomodate any types other than PixelDigis
 *           in case there is such a need in the future.
 *
 *  \author Nicola Pozzobon
 *  \date   2013, Jul 18
 *
 */

#ifndef L1TK_TRACK_BUILDER_H
#define L1TK_TRACK_BUILDER_H

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "L1Trigger/TrackTrigger/interface/TTTrackAlgorithm.h"
#include "L1Trigger/TrackTrigger/interface/TTTrackAlgorithmRecord.h"

#include <memory>
#include <map>
#include <vector>

template<  typename T  >
class TTTrackBuilder : public edm::EDProducer
{
  public:
    /// Constructor
    explicit TTTrackBuilder( const edm::ParameterSet& iConfig );

    /// Destructor;
    ~TTTrackBuilder();

  private:
    /// Tracking algorithm
    const StackedTrackerGeometry           *theStackedTracker;
    edm::ESHandle< TTTrackAlgorithm< T > > theTrackFindingAlgoHandle;
    edm::InputTag                          TTStubsInputTag;

    /// Other stuff
    bool enterAssociativeMemoriesWorkflow;

    /// Mandatory methods
    virtual void beginRun( const edm::Run& run, const edm::EventSetup& iSetup );
    virtual void endRun( const edm::Run& run, const edm::EventSetup& iSetup );
    virtual void produce( edm::Event& iEvent, const edm::EventSetup& iSetup );

}; /// Close class

/*! \brief   Implementation of methods
 *  \details Here, in the header file, the methods which do not depend
 *           on the specific type <T> that can fit the template.
 *           Other methods, with type-specific features, are implemented
 *           in the source file.
 */

/// Constructors
template< typename T >
TTTrackBuilder< T >::TTTrackBuilder( const edm::ParameterSet& iConfig )
{
  produces< std::vector< TTTrack< T > > >( "Seeds" );
  produces< std::vector< TTTrack< T > > >( "NoDup" );
  TTStubsInputTag = iConfig.getParameter< edm::InputTag >( "TTStubsBricks" );
  enterAssociativeMemoriesWorkflow = iConfig.getParameter< bool >( "AssociativeMemories" );
}

/// Destructor
template< typename T >
TTTrackBuilder< T >::~TTTrackBuilder() {}

/// Begin run
template< typename T >
void TTTrackBuilder< T >::beginRun( const edm::Run& run, const edm::EventSetup& iSetup )
{
  /// Get the geometry references
  edm::ESHandle< StackedTrackerGeometry > StackedTrackerGeomHandle;
  iSetup.get< StackedTrackerGeometryRecord >().get( StackedTrackerGeomHandle );
  theStackedTracker = StackedTrackerGeomHandle.product();

  /// Get the tracking algorithm 
  iSetup.get< TTTrackAlgorithmRecord >().get( theTrackFindingAlgoHandle );
  /// Print some information when loaded
  std::cout  << std::endl;
  std::cout  << "TTTrackBuilder< " << templateNameFinder< T >() << " > loaded modules:"
             << "\n\tTTTrackAlgorithm:\t" << theTrackFindingAlgoHandle->AlgorithmName()
             << std::endl;
  std::cout  << std::endl;
}

/// End run
template< typename T >
void TTTrackBuilder< T >::endRun( const edm::Run& run, const edm::EventSetup& iSetup ) {}

/// Implement the producer
template< typename T >
void TTTrackBuilder< T >::produce( edm::Event& iEvent, const edm::EventSetup& iSetup )
{
  /// Prepare output
  /// The temporary collection is used to store tracks
  /// before removal of duplicates
  std::auto_ptr< std::vector< TTTrack< T > > > TTTracksSeedsForOutput( new std::vector< TTTrack< T > > );
  std::auto_ptr< std::vector< TTTrack< T > > > TTTracksForOutput( new std::vector< TTTrack< T > > );
  std::auto_ptr< std::vector< TTTrack< T > > > TTTracksForOutputPurged( new std::vector< TTTrack< T > > );
  std::vector< TTTrack< T > > *tempTrackCollection = new std::vector< TTTrack< T > >();
  std::vector< TTTrack< T > > *tempPurgedCollection = new std::vector< TTTrack< T > >();

  /// Get the Stubs already stored away
  edm::Handle< std::vector< TTStub< T > > > TTStubHandle;
  iEvent.getByLabel( TTStubsInputTag, TTStubHandle );

  if ( enterAssociativeMemoriesWorkflow )
  {
    /// Enter AM 
    std::cerr << "TEST: AM workflow" << std::endl;
    theTrackFindingAlgoHandle->PatternFinding();
    theTrackFindingAlgoHandle->PatternRecognition();

  } /// End AM workflow
  else
  {
    /// Tracklet-based approach

    /// Map the Stubs per Sector/Wedge
    std::map< std::pair< unsigned int, unsigned int >, TKSector< T > > *trackerSectorMap;
    trackerSectorMap = new std::map< std::pair< unsigned int, unsigned int >, TKSector< T > >();
    theTrackFindingAlgoHandle->FillSectors( trackerSectorMap, TTStubHandle );

    /// Here trackerSectorMap contains all the Sectors of the event
    /// Each sector contains a map of Barrels and Endcaps

    /// Build the Seeds
    theTrackFindingAlgoHandle->FindSeeds( tempTrackCollection, trackerSectorMap );

    /// Put the seeds into the output
    TTTracksSeedsForOutput->assign( tempTrackCollection->begin(), tempTrackCollection->end() );
    tempTrackCollection->clear();

/*
    /// Get the number of sectors
    unsigned int nSectors = theTrackFindingAlgoHandle->ReturnNumberOfSectors();
    unsigned int nWedges = theTrackFindingAlgoHandle->ReturnNumberOfWedges();
*/

    /// Here trackerSectorMap contains all the Sectors of the event
    /// Each sector contains a map of Barrels and Endcaps
    /// TTTracksSeedsForOutput contains all the seeds built so far
    /// but there is also a copy of such seeds within the sectors
    /// which is the one which is actually propagated while
    /// the copies stored in TTTracksSeedsForOutput are safe
    theTrackFindingAlgoHandle->FindMatches( tempTrackCollection, trackerSectorMap );

    /// Fit the tracks
    typename std::vector< TTTrack< T > >::iterator trackIter;
    for ( trackIter = tempTrackCollection->begin();
          trackIter != tempTrackCollection->end();
          ++trackIter )
    {
      /// Here the seed is completed with all its matched stubs
      /// The seed is now a track and it is time to fit it
      theTrackFindingAlgoHandle->FitTrack( *trackIter );

      /// Refit tracks if needed
      if ( trackIter->getLargestResIdx() > -1 )
      {
        if ( trackIter->getStubPtrs().size() > 3 && trackIter->getChi2() > 100.0 )
        {
          trackIter->removeStubPtr( trackIter->getLargestResIdx() );
          trackIter->setLargestResIdx( -1 );
          theTrackFindingAlgoHandle->FitTrack( *trackIter );
        }
      } /// End of refit
    } /// End of loop over tracks

    /// Put the cleaned tracks into the output
    TTTracksForOutput->assign( tempTrackCollection->begin(), tempTrackCollection->end() );

  } /// End of non-AM

  theTrackFindingAlgoHandle->RemoveDuplicates( tempPurgedCollection, tempTrackCollection );

  /// Put the cleaned tracks into the output
  TTTracksForOutputPurged->assign( tempPurgedCollection->begin(), tempPurgedCollection->end() );

  /// Put in the event content
  iEvent.put( TTTracksSeedsForOutput, "Seeds" );
  //iEvent.put( TTTracksForOutput, "All" );
  iEvent.put( TTTracksForOutputPurged, "NoDup" );
}

#endif

