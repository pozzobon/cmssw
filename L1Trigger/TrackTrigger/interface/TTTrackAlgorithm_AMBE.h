/*! \class   TTTrackAlgorithm_AMBE
 *  \brief   Skeleton for AM-based track finder simulation.
 *  \details After moving from SimDataFormats to DataFormats,
 *           the template structure of the class was maintained
 *           in order to accomodate any types other than PixelDigis
 *           in case there is such a need in the future.
 *
 *  \author Nicola Pozzobon
 *  \author Sebastien Viret
 *  \date   2013, Jul 18
 *
 */

#ifndef L1_TRACK_TRIGGER_TRACK_ALGO_ASSOBE_H
#define L1_TRACK_TRIGGER_TRACK_ALGO_ASSOBE_H

#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/ModuleFactory.h"
#include "FWCore/Framework/interface/ESProducer.h"

#include "L1Trigger/TrackTrigger/interface/TTStubAlgorithm.h"
#include "L1Trigger/TrackTrigger/interface/TTTrackAlgorithm.h"
#include "L1Trigger/TrackTrigger/interface/TTTrackAlgorithmRecord.h"

#include "DataFormats/GeometryCommonDetAlgo/interface/MeasurementPoint.h"
#include "Geometry/CommonTopologies/interface/Topology.h"
#include "Geometry/TrackerGeometryBuilder/interface/PixelGeomDetUnit.h"
#include "Geometry/CommonTopologies/interface/PixelTopology.h"

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/shared_ptr.hpp>
#include <memory>
#include <string>
#include <map>
#include <vector>
#include <iostream>
#include <fstream>

#include "L1Trigger/TrackTrigger/interface/CMSPatternLayer.h"
#include "L1Trigger/TrackTrigger/interface/PatternFinder.h"
#include "L1Trigger/TrackTrigger/interface/SectorTree.h"
#include "L1Trigger/TrackTrigger/interface/Hit.h"

#ifndef __APPLE__
BOOST_CLASS_EXPORT_IMPLEMENT(CMSPatternLayer) 
#endif

template< typename T >
class TTTrackAlgorithm_AMBE : public TTTrackAlgorithm< T >
{
  private :
    /// Data members
    double        mMagneticField;
    unsigned int  nSectors;
    unsigned int  nWedges;
    std::string   nBKName;
    int           nThresh;
    SectorTree    m_st; //NP7 cambiare nome
    PatternFinder *m_pf; //NP7 cambiare nome

    //    std::vector<Hit*> m_hits;

  public:
    /// Constructors
    TTTrackAlgorithm_AMBE( const StackedTrackerGeometry *aStackedGeom,
                           double aMagneticField,
                           unsigned int aSectors,
                           unsigned int aWedges,
                           std::string aBKName,
                           int aThresh )
      : TTTrackAlgorithm< T > ( aStackedGeom, __func__ )
    {
      mMagneticField = aMagneticField;
      nSectors = aSectors;
      nWedges  = aWedges;
      nBKName  = aBKName;
      nThresh  = aThresh;

      std::cout << "Loading pattern bank file : " << std::endl;
      std::cout << nBKName << std::endl;

      std::ifstream ifs(nBKName.c_str());
      boost::archive::text_iarchive ia(ifs);

      ia >> m_st;
      m_pf = new PatternFinder( m_st.getSuperStripSize(), nThresh, &m_st, "", "" );
    }

    /// Destructor
    ~TTTrackAlgorithm_AMBE(){}

    /// Bank Loading 
    void LoadBank() const;

    /// Pattern Finding
    void PatternFinding( std::vector< TTTrack< T > > &output,                      
                         edm::Handle< std::vector< TTStub< T > > > &input ) const;

    /// Return the number of Sectors
    unsigned int ReturnNumberOfSectors() const { return nSectors; } /// Phi
    unsigned int ReturnNumberOfWedges() const  { return nWedges; } /// Eta

    /// Return the value of the magnetic field
    double ReturnMagneticField() const { return mMagneticField; }

    /// Fit the Track
    void FitTrack( TTTrack< T > &track ) const;

}; /// Close class

/*! \brief   Implementation of methods
 *  \details Here, in the header file, the methods which do not depend
 *           on the specific type <T> that can fit the template.
 *           Other methods, with type-specific features, are implemented
 *           in the source file.
 */

/// Find the patterns
template< >
void TTTrackAlgorithm_AMBE< Ref_PixelDigi_ >::PatternFinding( std::vector< TTTrack< Ref_PixelDigi_ > > &output,
                                                              edm::Handle< std::vector< TTStub< Ref_PixelDigi_ > > > &input ) const;

/// Fit the track
template< typename T >
void TTTrackAlgorithm_AMBE< T >::FitTrack( TTTrack< T > &track ) const
{
  std::cerr << "HOUGH!!!" << std::endl;
}





/*! \class   ES_TTTrackAlgorithm_AMBE
 *  \brief   Class to declare the algorithm to the framework
 *
 *  \author Nicola Pozzobon
 *  \date   2013, Jul 18
 *
 */

template< typename T >
class ES_TTTrackAlgorithm_AMBE : public edm::ESProducer
{
  private:
    /// Data members
    boost::shared_ptr< TTTrackAlgorithm< T > > _theAlgo;

    /// Number of Sectors
    unsigned int mSectors;
    unsigned int mWedges;
    std::string mBKName;
    int mThresh;

  public:
    /// Constructor
    ES_TTTrackAlgorithm_AMBE( const edm::ParameterSet & p )
      : mSectors( p.getParameter< int >("NumSectors") ),
        mWedges( p.getParameter< int >("NumWedges") ),
        mBKName( p.getParameter< std::string >("inputBankFile") ),
        mThresh( p.getParameter< int >("threshold") )
    {
      setWhatProduced( this );
    }

    /// Destructor
    virtual ~ES_TTTrackAlgorithm_AMBE() {}

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
        new TTTrackAlgorithm_AMBE< T >( &(*StackedTrackerGeomHandle),
                                        mMagneticFieldRounded,
                                        mSectors,
                                        mWedges,
                                        mBKName,
                                        mThresh );

      _theAlgo = boost::shared_ptr< TTTrackAlgorithm< T > >( TTTrackAlgo );
      return _theAlgo;
    }

};

#endif

