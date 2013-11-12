/*! \class   TTTrackAlgorithm
 *  \brief   Base class for any algorithm to be used
 *           in TTTrackBuilder
 *  \details After moving from SimDataFormats to DataFormats,
 *           the template structure of the class was maintained
 *           in order to accomodate any types other than PixelDigis
 *           in case there is such a need in the future.
 *
 *  \author Nicola Pozzobon
 *  \date   2013, Jul 18
 *
 */

#ifndef L1_TRACK_TRIGGER_TRACK_ALGO_BASE_H
#define L1_TRACK_TRIGGER_TRACK_ALGO_BASE_H

#include "DataFormats/L1TrackTrigger/interface/TTTypes.h"
#include "DataFormats/L1TrackTrigger/interface/TTTrack.h"
#include "DataFormats/GeometryCommonDetAlgo/interface/MeasurementPoint.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "Geometry/CommonDetUnit/interface/GeomDetUnit.h"
#include "Geometry/CommonTopologies/interface/Topology.h"
#include "Geometry/CommonTopologies/interface/PixelTopology.h"
#include "Geometry/TrackerGeometryBuilder/interface/PixelGeomDetUnit.h"
#include "Geometry/TrackerGeometryBuilder/interface/StackedTrackerGeometry.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "L1Trigger/TrackTrigger/interface/TKSector.h"

#include <sstream>
#include <map>
#include <string>
#include "classNameFinder.h"

template< typename T >
class TTTrackAlgorithm
{
  protected:
    /// Data members
    const StackedTrackerGeometry *theStackedTracker;
    std::string                  className_;

  public:
    /// Constructors
    TTTrackAlgorithm( const StackedTrackerGeometry *aStackedGeom, std::string fName )
      : theStackedTracker( aStackedGeom )
    {
      className_ = classNameFinder< T >(fName);
    }

    /// Destructor
    virtual ~TTTrackAlgorithm(){}

    /// Fill sectors with stubs and structures
    virtual void FillSectors( std::map< std::pair< unsigned int, unsigned int >, TKSector< T > > *outputSectorMap,
                              edm::Handle< std::vector< TTStub< Ref_PixelDigi_ > > > &input ) const
    {
      outputSectorMap->clear();
    }

    /// Find the Seeds
    virtual void FindSeeds( std::vector< TTTrack< T > > *output,
                            std::map< std::pair< unsigned int, unsigned int >, TKSector< T > > *inputSectorMap ) const
    {
      output->clear();
    }

    /// Propagate the Seeds
    virtual void FindMatches( std::vector< TTTrack< T > > *output,
                              std::map< std::pair< unsigned int, unsigned int >, TKSector< T > > *inputSectorMap ) const
    {
      output->clear();
    }

    /// Remove the duplicates
    virtual void RemoveDuplicates( std::vector< TTTrack< T > > *output, std::vector< TTTrack< T > > *input ) const;

    /// AM Pattern Finding
    virtual void PatternFinding() const
    {}

    /// AM Pattern Recognition
    virtual void PatternRecognition() const
    {}

    virtual unsigned int ReturnNumberOfSectors() const { return 1; } 
    virtual unsigned int ReturnNumberOfWedges() const  { return 1; }
    virtual double ReturnMagneticField() const         { return 1.0; }

    /// Fit the Track
    virtual void FitTrack( TTTrack< T > &seed ) const;

    /// Algorithm name
    virtual std::string AlgorithmName() const { return className_; }

    /// Helper methods
    double DeltaPhi( double phi1, double phi2 ) const;
    double CosineTheorem( double a, double b, double phi ) const;
    double FindRInvOver2( double rho1, double rho2, double phi1, double phi2 ) const;

}; /// Close class





/*! \brief   Implementation of methods
 *  \details Here, in the header file, the methods which do not depend
 *           on the specific type <T> that can fit the template.
 *           Other methods, with type-specific features, are implemented
 *           in the source file.
 */

/// Remove the duplicates
template< typename T >
void TTTrackAlgorithm< T >::RemoveDuplicates( std::vector< TTTrack< T > > *output, std::vector< TTTrack< T > > *input ) const
{
  /// Prepare output
  output->clear();

  /// Prepare the vector of booleans to matk the tracks to be deleted
  std::vector< bool > toBeDeleted;
  toBeDeleted.assign( input->size(), false );

  for ( unsigned int i = 0; i < input->size(); i++ )
  {
    /// This check is necessary as the bool may be reset in a previous iteration
    if ( toBeDeleted.at(i) )
      continue;

    /// Check if the track has min 3 stubs
    if ( input->at(i).getStubPtrs().size() < 3 )
      continue;

    /// Count the number of PS stubs
    unsigned int nPSi = 0;
    for ( unsigned int is = 0; is < input->at(i).getStubPtrs().size(); is++ )
    {
      StackedTrackerDetId stDetId( input->at(i).getStubPtrs().at(is)->getDetId() );
      bool isPS = theStackedTracker->isPSModule( stDetId );
      if ( isPS )
        nPSi++;
    }

    bool hasBL1i = input->at(i).hasStubInBarrel(1);

    /// Nested loop to compare tracks with each other
    for ( unsigned int j = i+1 ; j < input->size(); j++ )
    {
      /// This check is necessary as the bool may be reset in a previous iteration
      if ( toBeDeleted.at(j) )
        continue;

      /// Check if they are the same track
      if ( input->at(i).isTheSameAs( input->at(j) ) )
      {
        /// Check if the track has min 3 stubs
        if ( input->at(j).getStubPtrs().size() < 3 )
          continue;

        /// Count the number of PS stubs
        unsigned int nPSj = 0;
        for ( unsigned int js = 0; js < input->at(j).getStubPtrs().size(); js++ )
        {
          StackedTrackerDetId stDetId( input->at(j).getStubPtrs().at(js)->getDetId() );
          bool isPS = theStackedTracker->isPSModule( stDetId );
          if ( isPS )
            nPSj++;
        }

        /// Choose the one with the largest number of PS stubs
        if ( nPSi > nPSj )
        {
          toBeDeleted[j] = true;
          continue;
        }
        else if ( nPSi < nPSj )
        {
          toBeDeleted[i] = true;
          continue;
        }
        /// Here we are if the two tracks have the same number of PS stubs

        /// Check which one has a stub in Barrel L1
        bool hasBL1j = input->at(j).hasStubInBarrel(1);

        if ( hasBL1i || hasBL1j )
        {
          if ( !hasBL1i )
          {
            toBeDeleted[i] = true;
            continue;
          }
          if ( !hasBL1j )
          {
            toBeDeleted[j] = true;
            continue;
          }
        }
        /// We get here only if both have BL1 or both have not

        /// Compare Chi2
        if ( input->at(i).getChi2Red() > input->at(j).getChi2Red() )
        {
          toBeDeleted[i] = true;
        }
        else
        {
          toBeDeleted[j] = true;
        }
        continue;
      }
    }
  } /// End of loop to assign the boolean flags

  /// Store only the non-deleted tracks
  /// (already tried using vector::erase, but something
  /// odd had happened ... that's why it is this way)
  for ( unsigned int i = 0; i < input->size(); i++ )
  {
    if ( toBeDeleted.at(i) )
      continue;

    output->push_back( input->at(i) );
  }

}

/// Fit the track
template< >
void TTTrackAlgorithm< Ref_PixelDigi_ >::FitTrack( TTTrack< Ref_PixelDigi_ > &seed ) const;

/// Helper methods
template< typename T >
double TTTrackAlgorithm< T >::DeltaPhi( double phi1, double phi2 ) const
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
double TTTrackAlgorithm< T >::CosineTheorem( double a, double b, double phi ) const
{
  return sqrt( a*a + b*b - 2*a*b*cos(phi) );
}

template< typename T >
double TTTrackAlgorithm< T >::FindRInvOver2( double rho1, double rho2, double phi1, double phi2 ) const
{
  /// Calculate angle between the vectors
  double deltaPhi = this->DeltaPhi( phi1, phi2 );

  /// Apply cosine theorem to find the distance between the vectors
  double distance = this->CosineTheorem( rho1, rho2, deltaPhi );

  /// Apply sine theorem to find 1/(2R)
  return sin(deltaPhi)/distance; /// Sign is maintained to keep track of the charge
}

#endif

