/*! \class   TKBarrel
 *  \brief   Dedicated Barrel class for Tracklet-based
 *           geometry containers
 *
 *  \author Nicola Pozzobon
 *  \author Anders Ryd
 *  \date   2013, Nov 8
 *
 */

#ifndef L1_TRACK_TRIGGER_TKBARREL_H
#define L1_TRACK_TRIGGER_TKBARREL_H

#include "L1Trigger/TrackTrigger/interface/TKGeomBase.h"

template< typename T >
class TKEndcap;

template< typename T >
class TKBarrel : public TKGeomBase< T >
{
  protected:
    /// Data members
    const StackedTrackerGeometry *theStackedTracker;
    double                       mMagneticField;

  private:
    /// Data members
    std::vector< edm::Ptr< TTStub< T > > > *theStubPtrs;
    std::vector< TTTrack< T > >            *theTracks;

  public:
    /// Constructors
    TKBarrel(){}

    TKBarrel( const StackedTrackerGeometry *aStackedGeom, double aMagneticField )
      : theStackedTracker( aStackedGeom ),
        mMagneticField( aMagneticField )
    {
      theStubPtrs = new std::vector< edm::Ptr< TTStub< T > > >();
      theTracks = new std::vector< TTTrack< T > >();
    }

    /// Destructor
    ~TKBarrel(){}

    /// Stubs
    std::vector< edm::Ptr< TTStub< T > > > getStubPtrs() const { return *theStubPtrs; }
    void addStubPtr( edm::Ptr< TTStub< T > > aStubPtr )        { theStubPtrs->push_back( aStubPtr ); }

    /// Tracks
    std::vector< TTTrack< T > > getTracks() const { return *theTracks; }
    void addTrack( TTTrack< T > aTrack )          { theTracks->push_back( aTrack ); }

    /// Find Tracklets
    void findSeeds( TKBarrel< T > aBarrel, unsigned int aSector, unsigned int aWedge, bool wantHermetic = false ) const {}
    void findSeeds( TKEndcap< T > anEndcap, unsigned int aSector, unsigned int aWedge ) const {}


    /// Find Matches
    void findMatches( TKBarrel< T > aBarrel, double rPhiWindow, double zWindow ) const {}
    void findMatches( TKEndcap< T > anEndcap, std::pair< double, double > rPhiWindow, std::pair< double, double > rWindow  ) const {}


}; /// Close class

template< >
void TKBarrel< Ref_PixelDigi_ >::findSeeds( TKBarrel< Ref_PixelDigi_ > aBarrel, unsigned int aSector, unsigned int aWedge, bool wantHermetic ) const;

template< >
void TKBarrel< Ref_PixelDigi_ >::findSeeds( TKEndcap< Ref_PixelDigi_ > anEndcap, unsigned int aSector, unsigned int aWedge  ) const;

template< >
void TKBarrel< Ref_PixelDigi_ >::findMatches( TKBarrel< Ref_PixelDigi_ > aBarrel, double rPhiWindow, double zWindow ) const;

template< >
void TKBarrel< Ref_PixelDigi_ >::findMatches( TKEndcap< Ref_PixelDigi_ > anEndcap, std::pair< double, double > rPhiWindow, std::pair< double, double > rWindow  ) const;

#endif

