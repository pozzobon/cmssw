/*! \class   TKSector
 *  \brief   Base class Tracklet-related geometry containers
 *
 *  \author Nicola Pozzobon
 *  \author Anders Ryd
 *  \date   2013, Nov 8
 *
 */

#ifndef L1_TRACK_TRIGGER_TKSECTOR_H
#define L1_TRACK_TRIGGER_TKSECTOR_H

#include "L1Trigger/TrackTrigger/interface/TKBarrel.h"
#include "L1Trigger/TrackTrigger/interface/TKEndcap.h"

template< typename T >
class TKSector
{
  protected:
    /// Data members
    const StackedTrackerGeometry *theStackedTracker;

  private:
    /// Data members
    std::map< unsigned int, TKBarrel< T > >                            theBarrelLayerMap;
    std::map< std::pair< unsigned int, unsigned int >, TKEndcap< T > > theEndcapDiskMap;

  public:
    /// Constructors
    TKSector(){}

    TKSector( const StackedTrackerGeometry *aStackedGeom )
    {
      theStackedTracker = aStackedGeom;
      theBarrelLayerMap.clear();
      theEndcapDiskMap.clear();
    }

    /// Destructor
    ~TKSector(){}

    /// Access to data members
    std::map< unsigned int, TKBarrel< T > > getBarrelLayerMap() const                           { return theBarrelLayerMap; }
    std::map< std::pair< unsigned int, unsigned int >, TKEndcap< T > > getEndcapDiskMap() const { return theEndcapDiskMap; }

    /// Add elements
    void addBarrel( unsigned int aLayer, TKBarrel< T > aBarrel )                            { theBarrelLayerMap.insert( std::make_pair( aLayer, aBarrel ) ); }
    void addEndcap( std::pair< unsigned int, unsigned int > aDisk, TKEndcap< T > anEndcap ) { theEndcapDiskMap.insert( std::make_pair( aDisk, anEndcap ) ); }

    void addStubPtrToBarrel( unsigned int aLayer, edm::Ptr< TTStub< T > > aStubPtr )                           { theBarrelLayerMap[aLayer].addStubPtr( aStubPtr ); }
    void addStubPtrToEndcap( std::pair< unsigned int, unsigned int > aDisk, edm::Ptr< TTStub< T > > aStubPtr ) { theEndcapDiskMap[aDisk].addStubPtr( aStubPtr ); }

}; /// Close class

#endif

