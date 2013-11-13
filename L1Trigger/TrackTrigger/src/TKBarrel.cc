/*! \class   TKBarrel
 *  \brief   Implementation of the specific methods of TKBarrel
 *  \details Here, in the source file, the methods which do depend
 *           on the specific type <T> that can fit the template.
 *
 *  \author Nicola Pozzobon
 *  \author Anders Ryd
 *  \date   2013, Nov 8
 *
 */

#include "L1Trigger/TrackTrigger/interface/TKBarrel.h"
#include "L1Trigger/TrackTrigger/interface/TKEndcap.h"

template< >
void TKBarrel< Ref_PixelDigi_ >::findSeeds( TKBarrel< Ref_PixelDigi_ > aBarrel, unsigned int aSector, unsigned int aWedge, bool wantHermetic ) const
{
  /// The Stubs in this Barrel are stored in theStubPtrs

  /// Get the Stubs from the second Barrel
  std::vector< edm::Ptr< TTStub< Ref_PixelDigi_ > > > otherStubPtrs = aBarrel.getStubPtrs();

  /// Double loop over all pairs of stubs
  for ( unsigned int i = 0; i < theStubPtrs->size(); i++ )
  {
    edm::Ptr< TTStub< Ref_PixelDigi_ > > thisStubPtr = theStubPtrs->at(i);

    StackedTrackerDetId detId1( thisStubPtr->getDetId() );
    GlobalPoint pos1 = theStackedTracker->findGlobalPosition( thisStubPtr.get() );
    double rho1 = pos1.perp();
    double phi1 = pos1.phi();
    double z1 = pos1.z();

    /// Layer-Disk-Ring constraint
    if ( wantHermetic )
    {
      if ( detId1.iLayer() % 2 == 0 )
        continue;
    }
    else
    {
      if ( detId1.iLayer() > 4 )
        continue;
    }

    bool barrelSeed2S = !( theStackedTracker->isPSModule( detId1 ) );

    /// Find the index of the first element of the nested loop
    unsigned int startIndex = 0;
    if ( this == &aBarrel )
    {
      /// If they are in the same Sector/Wedge the loop can be simplified
      /// This means theStubPtrs is the same as otherStubPtrs
      startIndex = i + 1;
    }

    /// Loop over other Barrel Stubs
    for ( unsigned int k = startIndex; k < otherStubPtrs.size(); k++ )
    {
      edm::Ptr< TTStub< Ref_PixelDigi_ > > otherStubPtr = otherStubPtrs.at(k);

      StackedTrackerDetId detId2( otherStubPtr->getDetId() );
      GlobalPoint pos2 = theStackedTracker->findAverageGlobalPosition( otherStubPtr->getClusterPtr(0).get() );
      double rho2 = pos2.perp();
      double phi2 = pos2.phi();
      double z2 = pos2.z();

      /// Skip same layer pairs
      /// Skip pairs with distance larger than 1 layer
      if ( detId2.iLayer() != detId1.iLayer() + 1 )
        continue;

      if ( wantHermetic )
      {
        if ( detId2.iPhi() != detId1.iPhi() )
          continue;
      }

      /// Safety cross-check
      if ( rho1 > rho2 )
      {
        continue;
      }

      /// Apply cosine theorem to find 1/(2R) = rInvOver2
#include "L1Trigger/TrackTrigger/src/TTTrackAlgorithm_trackletBE_SeedTrigonometry.icc"

      /// Perform cut on Pt
      if ( fabs(rInvOver2) > mMagneticField*0.0015*0.5 )
        continue;

      /// Calculate tracklet parameters with helix model
      /// roughPt, z0, cotTheta0, phi0
#include "L1Trigger/TrackTrigger/src/TTTrackAlgorithm_trackletBE_SeedParameters.icc"

      /// Correct for seeds in 2S Barrel layers
      if ( barrelSeed2S )
      {
        if ( fabs( z1 - z2 ) < 10 )
        {
          z0 = 0;
          cotTheta0 = z1 / rhoPsi1;
        }
      }

      /// Perform projected vertex cut
      if ( fabs(z0) > 30.0 )
        continue;

      /// Create the Seed in the form of a Track and store it in the output
      std::vector< edm::Ptr< TTStub< Ref_PixelDigi_ > > > tempVec;
      tempVec.push_back( thisStubPtr );
      tempVec.push_back( otherStubPtr );
      TTTrack< Ref_PixelDigi_ > tempTrack( tempVec );
      tempTrack.setSector( aSector );
      tempTrack.setWedge( aWedge );
      tempTrack.setRInv( 2*rInvOver2 );
      tempTrack.setMomentum( GlobalVector( roughPt*cos(phi0),
                                           roughPt*sin(phi0),
                                           roughPt*cotTheta0 ) );
      tempTrack.setVertex( GlobalPoint( 0, 0, z0 ) );
      theTracks->push_back( tempTrack );
    } /// End of loop over other Barrel Stubs
  } /// End of loop over these Stubs
}

template< >
void TKBarrel< Ref_PixelDigi_ >::findSeeds( TKEndcap< Ref_PixelDigi_ > anEndcap, unsigned int aSector, unsigned int aWedge ) const
{
  /// The Stubs in this Barrel are stored in theStubPtrs

  /// Get the Stubs from the Endcap
  std::vector< edm::Ptr< TTStub< Ref_PixelDigi_ > > > otherStubPtrs = anEndcap.getStubPtrs();

  /// Double loop over all pairs of stubs
  for ( unsigned int i = 0; i < theStubPtrs->size(); i++ )
  {
    edm::Ptr< TTStub< Ref_PixelDigi_ > > thisStubPtr = theStubPtrs->at(i);

    StackedTrackerDetId detId1( thisStubPtr->getDetId() );
    GlobalPoint pos1 = theStackedTracker->findAverageGlobalPosition( thisStubPtr->getClusterPtr(0).get() );
    double rho1 = pos1.perp();
    double phi1 = pos1.phi();
    double z1 = pos1.z();

    /// Layer-Disk-Ring constraint
    if ( detId1.iLayer() > 4 )
      continue;

    /// Loop over Endcap Stubs
    for ( unsigned int k = 0; k < otherStubPtrs.size(); k++ )
    {
      edm::Ptr< TTStub< Ref_PixelDigi_ > > otherStubPtr = otherStubPtrs.at(k);

      StackedTrackerDetId detId2( otherStubPtr->getDetId() );
      GlobalPoint pos2 = theStackedTracker->findAverageGlobalPosition( otherStubPtr->getClusterPtr(0).get() );
      double rho2 = pos2.perp();
      double phi2 = pos2.phi();
      double z2 = pos2.z();

      /// Layer-Disk-Ring constraint
      if ( detId1.iLayer() > 4 )
        continue;

      if ( detId2.iDisk() > 1 || detId2.iRing() > 11 )
        continue;

      if ( detId1.iLayer() == 1 && detId2.iRing() < 1 )
      {
        continue;
      }
      else if ( detId1.iLayer() == 2 && detId2.iRing() < 4 )
      {
        continue;
      }
      else if ( detId1.iLayer() == 3 && detId2.iRing() < 8 )
      {
        continue;
      }

      /// Safety cross-check
      if ( rho1 > rho2 )
      {
        continue;
      }

      /// Apply cosine theorem to find 1/(2R) = rInvOver2
#include "L1Trigger/TrackTrigger/src/TTTrackAlgorithm_trackletBE_SeedTrigonometry.icc"

      /// Perform cut on Pt
      if ( fabs(rInvOver2) > mMagneticField*0.0015*0.5 )
        continue;

      /// Calculate tracklet parameters with helix model
      /// roughPt, z0, cotTheta0, phi0
#include "L1Trigger/TrackTrigger/src/TTTrackAlgorithm_trackletBE_SeedParameters.icc"

      /// Correct for Endcap 2S in the seed
      bool endcapSeedPS = theStackedTracker->isPSModule( detId2 );

      if ( !endcapSeedPS )
      {
        if ( fabs( z1 - z2 ) < 10 )
        {
          z0 = 0;
          cotTheta0 = z1 / rhoPsi1; /// Use barrel coordinates!
        }
      }

      /// Perform projected vertex cut
      if ( fabs(z0) > 30.0 )
        continue;

      /// Create the Seed in the form of a Track and store it in the output
      std::vector< edm::Ptr< TTStub< Ref_PixelDigi_ > > > tempVec;
      tempVec.push_back( thisStubPtr );
      tempVec.push_back( otherStubPtr );
      TTTrack< Ref_PixelDigi_ > tempTrack( tempVec );
      tempTrack.setSector( aSector );
      tempTrack.setWedge( aWedge );
      tempTrack.setRInv( 2*rInvOver2 );
      tempTrack.setMomentum( GlobalVector( roughPt*cos(phi0),
                                           roughPt*sin(phi0),
                                           roughPt*cotTheta0 ) );
      tempTrack.setVertex( GlobalPoint( 0, 0, z0 ) );
      theTracks->push_back( tempTrack );

    } /// End of loop over other Barrel Stubs
  } /// End of loop over these Stubs
}

template< >
void TKBarrel< Ref_PixelDigi_ >::findMatches( TKBarrel< Ref_PixelDigi_ > aBarrel, double rPhiWindow, double zWindow ) const
{
  /// The Stubs in this Barrel are stored in theStubPtrs
  /// The Seeds in this Barrel are stored in theTracks
  /// Moreover, by construction, the first two Stubs in each Track, are the seed ones!

  /// Get the Stubs from the second Barrel
  std::vector< edm::Ptr< TTStub< Ref_PixelDigi_ > > > otherStubPtrs = aBarrel.getStubPtrs();

  /// Loop over all the Seeds in the Container and over all the Stubs in the second Barrel
  for ( unsigned int it = 0; it < theTracks->size(); it++ )
  {
    /// Get the Seed Stubs
    std::vector< edm::Ptr< TTStub< Ref_PixelDigi_ > > > theSeedStubPtrs = theTracks->at(it).getStubPtrs();

    /// Get the Seed Momentum and Vtx
    GlobalVector curMomentum = theTracks->at(it).getMomentum();
    GlobalPoint curVertex = theTracks->at(it).getVertex();
    double curRInv = theTracks->at(it).getRInv();
    double curPhi = curMomentum.phi();
    double curTheta = curMomentum.theta();
    double curZVtx = curVertex.z();

    /// Loop over the Stubs to be attached
    for ( unsigned int is = 0; is < otherStubPtrs.size(); is++ )
    {
      edm::Ptr< TTStub< Ref_PixelDigi_ > > candidate = otherStubPtrs.at(is);

      /// Check that the candidate is NOT the one under examination
      if ( theSeedStubPtrs.at(0) == candidate ||
           theSeedStubPtrs.at(1) == candidate )
        continue;

      /// Skip if the Stub is in one of the Seed Layers/Disks  
      StackedTrackerDetId stDetId0( theSeedStubPtrs.at(0)->getDetId() );
      StackedTrackerDetId stDetId1( theSeedStubPtrs.at(1)->getDetId() );
      StackedTrackerDetId stDetIdCand( candidate->getDetId() );

      /// Case 1: Barrel-Barrel Seed towards Barrel Stub
      if ( stDetId0.isBarrel() && stDetId1.isBarrel() )
      {
        if ( stDetId0.iLayer() == stDetIdCand.iLayer() ||
             stDetId1.iLayer() == stDetIdCand.iLayer() )
          continue;
      }

      /// Case 2: Barrel-Endcap Seed towards Barrel Stub
      if ( stDetId0.isBarrel() && stDetId1.isEndcap() )
      {
        if ( stDetId0.iLayer() == stDetIdCand.iLayer() )
          continue;
      }

      /// Get the candidate Stub position
      GlobalPoint candPos = theStackedTracker->findGlobalPosition( candidate.get() );
      double rhoCand = candPos.perp();
      double phiCand = candPos.phi();
      double zCand = candPos.z();

      /// Calculate deltaRPhi and deltaZ
#include "L1Trigger/TrackTrigger/src/TTTrackAlgorithm_trackletBE_BarrelSeedPropagation.icc"

      if ( deltaRPhi < rPhiWindow &&
           deltaZ < zWindow )
      {
        theTracks->at(it).addStubPtr( candidate );
      }
    } /// End of loop over the Stubs

    /// Now the Seed is longer than earlier, but its first two Stubs are still
    /// the ones from the Tracklet finding

  } /// End of loop over the Seeds
}

template< >
void TKBarrel< Ref_PixelDigi_ >::findMatches( TKEndcap< Ref_PixelDigi_ > anEndcap,
                                              std::pair< double, double > rPhiWindow, std::pair< double, double > rWindow ) const
{
  /// The Stubs in this Barrel are stored in theStubPtrs
  /// The Seeds in this Barrel are stored in theTracks
  /// Moreover, by construction, the first two Stubs in each Track, are the seed ones!

  /// Get the Stubs from the Endcap
  std::vector< edm::Ptr< TTStub< Ref_PixelDigi_ > > > otherStubPtrs = anEndcap.getStubPtrs();

  /// Loop over all the Seeds in the Container and over all the Stubs in the Endcap
  for ( unsigned int it = 0; it < theTracks->size(); it++ )
  {
    /// Get the Seed Stubs
    std::vector< edm::Ptr< TTStub< Ref_PixelDigi_ > > > theSeedStubPtrs = theTracks->at(it).getStubPtrs();

    /// Get the Seed Momentum and Vtx
    GlobalVector curMomentum = theTracks->at(it).getMomentum();
    GlobalPoint curVertex = theTracks->at(it).getVertex();
    double curRInv = theTracks->at(it).getRInv();
    double curPhi = curMomentum.phi();
    double curTheta = curMomentum.theta();
    double curZVtx = curVertex.z();

    /// Loop over the Stubs to be attached
    for ( unsigned int is = 0; is < otherStubPtrs.size(); is++ )
    {
      edm::Ptr< TTStub< Ref_PixelDigi_ > > candidate = otherStubPtrs.at(is);

      /// Check that the candidate is NOT the one under examination
      /// By construction, this is possible only for the second Stub in
      /// the Seed, and only for Barrel-Endcap seeds
      if ( theSeedStubPtrs.at(1) == candidate )
        continue;

      /// Skip if the Stub is in one of the Seed Layers/Disks  
      StackedTrackerDetId stDetId0( theSeedStubPtrs.at(0)->getDetId() );
      StackedTrackerDetId stDetId1( theSeedStubPtrs.at(1)->getDetId() );
      StackedTrackerDetId stDetIdCand( candidate->getDetId() );

      /// Case 1: Barrel-Endcap Seed towards Endcap Stub
      if ( stDetId0.isBarrel() && stDetId1.isEndcap() )
      {
        if ( stDetId1.iSide() == stDetIdCand.iSide() )
        {
          if ( stDetId1.iDisk() == stDetIdCand.iDisk() )
            continue;
        }
      }

      bool endcapCandPS = stDetIdCand.isEndcap() && theStackedTracker->isPSModule( stDetIdCand );

      /// Get the candidate Stub position
      GlobalPoint candPos = theStackedTracker->findGlobalPosition( candidate.get() );
      double rhoCand = candPos.perp();
      double phiCand = candPos.phi();
      double zCand = candPos.z();

      /// Calculate a correction for non-pointing-strips in square modules
      /// Relevant angle is the one between hit and module center, with
      /// vertex at (0, 0). Take snippet from HitMatchingAlgorithm_window201*
      /// POSITION IN TERMS OF PITCH MULTIPLES:
      ///       0 1 2 3 4 5 6 7 8 9 ...
      /// COORD: 0 1 2 3 4 5 6 7 8 9 ...
      /// OUT   | | | | | |x| | | | | | | | | |
      ///
      /// IN    | | | |x|x| | | | | | | | | | |
      ///             THIS is 3.5 (COORD) and 4.0 (POS)
      /// The center of the module is at NROWS/2 (position) and NROWS-0.5 (coordinates)
      StackedTrackerDetId stDetId( candidate->getClusterPtr(0)->getDetId() );
      const GeomDetUnit* det0 = theStackedTracker->idToDetUnit( stDetId, 0 );
      const PixelGeomDetUnit* pix0 = dynamic_cast< const PixelGeomDetUnit* >( det0 );
      const PixelTopology* top0 = dynamic_cast< const PixelTopology* >( &(pix0->specificTopology()) );
      std::pair< float, float > pitch0 = top0->pitch();
      MeasurementPoint stubCoord = candidate->getClusterPtr(0)->findAverageLocalCoordinates();
      double stubTransvDispl = pitch0.first * ( stubCoord.x() - (top0->nrows()/2 - 0.5) ); /// Difference in coordinates is the same as difference in position

      if ( zCand > 0 )
      {
        stubTransvDispl = - stubTransvDispl;
      }

      /// Calculate deltaRPhi and deltaRho
#include "L1Trigger/TrackTrigger/src/TTTrackAlgorithm_trackletBE_EndcapSeedPropagation.icc"

      if ( endcapCandPS )
      {
        if ( deltaRPhi < rPhiWindow.first &&
             deltaR < rWindow.first )
        {
          theTracks->at(it).addStubPtr( candidate );
        }
      }
      else
      {
        if ( deltaRPhi < rPhiWindow.second &&
             deltaR < rWindow.second )
        {
          theTracks->at(it).addStubPtr( candidate );
        }
      }
    } /// End of loop over the Stubs

    /// Now the Seed is longer than earlier, but its first two Stubs are still
    /// the ones from the Tracklet finding

  } /// End of loop over the Seeds
}

