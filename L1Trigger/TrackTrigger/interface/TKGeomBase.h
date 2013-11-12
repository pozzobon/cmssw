/*! \class   TKGeomBase
 *  \brief   Base class Tracklet-related geometry containers
 *
 *  \author Nicola Pozzobon
 *  \author Anders Ryd
 *  \date   2013, Nov 8
 *
 */

#ifndef L1_TRACK_TRIGGER_TKGEOM_BASE_H
#define L1_TRACK_TRIGGER_TKGEOM_BASE_H

#include "DataFormats/L1TrackTrigger/interface/TTTypes.h"
#include "Geometry/TrackerGeometryBuilder/interface/StackedTrackerGeometry.h"

template< typename T >
class TKGeomBase
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
    TKGeomBase(){}

    TKGeomBase( const StackedTrackerGeometry *aStackedGeom, double aMagneticField )
      : theStackedTracker( aStackedGeom ),
        mMagneticField( aMagneticField )
    {
      theStubPtrs = new std::vector< edm::Ptr< TTStub< T > > >();
      theTracks = new std::vector< TTTrack< T > >();
    }

    /// Destructor
    virtual ~TKGeomBase(){}

    /// Stubs
    std::vector< edm::Ptr< TTStub< T > > > getStubPtrs() const { return *theStubPtrs; }

    /// Tracks
    std::vector< TTTrack< T > > getTracks() const { return *theTracks; }

    /// Find Tracklets
    virtual void findSeeds( TKGeomBase< T > ) const {}
    virtual void findMatches( TKGeomBase< T > ) const {}

    /// Helper methods
    double DeltaPhi( double phi1, double phi2 ) const;
    double CosineTheorem( double a, double b, double phi ) const;
    double FindRInvOver2( double rho1, double rho2, double phi1, double phi2 ) const;

}; /// Close class

/// Helper methods
template< typename T >
double TKGeomBase< T >::DeltaPhi( double phi1, double phi2 ) const
{
  double deltaPhi = phi1 - phi2;
  if ( fabs(deltaPhi) >= M_PI )
  {
    if ( deltaPhi>0 )
      deltaPhi = deltaPhi - 2*M_PI;
    else
      deltaPhi = 2*M_PI + deltaPhi;
  }
  return deltaPhi;
}

template< typename T >
double TKGeomBase< T >::CosineTheorem( double a, double b, double phi ) const
{
  return sqrt( a*a + b*b - 2*a*b*cos(phi) );
}

template< typename T >
double TKGeomBase< T >::FindRInvOver2( double rho1, double rho2, double phi1, double phi2 ) const
{
  /// Calculate angle between the vectors
  double deltaPhi = this->DeltaPhi( phi1, phi2 );

  /// Apply cosine theorem to find the distance between the vectors
  double distance = this->CosineTheorem( rho1, rho2, deltaPhi );

  /// Apply sine theorem to find 1/(2R)
  return sin(deltaPhi)/distance; /// Sign is maintained to keep track of the charge
}

#endif

