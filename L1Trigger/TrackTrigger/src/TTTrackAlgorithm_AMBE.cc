/*! \brief   Implementation of methods of TTTrackAlgorithm_AMBE
 *  \details Here, in the source file, the methods which do depend
 *           on the specific type <T> that can fit the template.
 *
 *  \author Seb Viret
 *  \author Guillaume Baulieu
 *  \author Nicola Pozzobon
 *  \date   2013, Oct 28
 *
 */

#include "L1Trigger/TrackTrigger/interface/TTTrackAlgorithm_AMBE.h"

/// Pattern finding
template< >
void TTTrackAlgorithm_AMBE< Ref_PixelDigi_ >::PatternFinding( std::vector< TTTrack< Ref_PixelDigi_ > > &output,
                                                              edm::Handle< std::vector< TTStub< Ref_PixelDigi_ > > > &input ) const
{

  /// STEP 0
  /// Prepare output
  output.clear();

  int layer  = 0;
  int ladder = 0;
  int module = 0;
  int disk   = 0;
  int lad_cor= 0;

  /// STEP 1
  /// Loop over input stubs

  //  std::cout << "Start the loop on stubs to transform them into sstrips" << std::endl;

  std::vector<Hit*> m_hits;

  for(unsigned int i=0;i<m_hits.size();i++)
  {
    delete m_hits[i];
  }
  
  m_hits.clear();

  typename std::vector< TTStub< Ref_PixelDigi_ > >::const_iterator inputIter;
  unsigned int j = 0; /// Counter needed to build the edm::Ptr to the TTStub
  for ( inputIter = input->begin();
        inputIter != input->end();
        ++inputIter )
  {

    /// Make the pointer to be put in the Track
    edm::Ptr< TTStub< Ref_PixelDigi_ > > tempStubPtr( input, j++ );

    /// Calculate average coordinates col/row for inner/outer Cluster
    /// These are already corrected for being at the center of each pixel
    MeasurementPoint mp0 = tempStubPtr->getClusterPtr(0)->findAverageLocalCoordinates();
    GlobalPoint posStub  = theStackedTracker->findGlobalPosition( &(*tempStubPtr) );

    StackedTrackerDetId detIdStub( tempStubPtr->getDetId() );

    //   bool isPS = TTTrackAlgorithm< Ref_PixelDigi_ >::theStackedTracker->isPSModule( detIdStub );

    const GeomDetUnit* det0 = TTTrackAlgorithm< Ref_PixelDigi_ >::theStackedTracker->idToDetUnit( detIdStub, 0 );
    const GeomDetUnit* det1 = TTTrackAlgorithm< Ref_PixelDigi_ >::theStackedTracker->idToDetUnit( detIdStub, 1 );

    /// Find pixel pitch and topology related information
    const PixelGeomDetUnit* pix0 = dynamic_cast< const PixelGeomDetUnit* >( det0 );
    const PixelGeomDetUnit* pix1 = dynamic_cast< const PixelGeomDetUnit* >( det1 );
    const PixelTopology* top0    = dynamic_cast< const PixelTopology* >( &(pix0->specificTopology()) );
    const PixelTopology* top1    = dynamic_cast< const PixelTopology* >( &(pix1->specificTopology()) );
    

    /// Stop if the clusters are not in the same z-segment
    int cols0   = top0->ncolumns();
    int cols1   = top1->ncolumns();
    int ratio   = cols0/cols1; /// This assumes the ratio is integer!
    int segment = floor( mp0.y() / ratio );

    // Here we rearrange the number in order to be compatible with the AM emulator

    if ( detIdStub.isBarrel() ) 
    {
      layer  = detIdStub.iLayer()+4;
      ladder = detIdStub.iPhi()-1;
      module = detIdStub.iZ()-1;
    }
    else if ( detIdStub.isEndcap() )
    {
      layer  = 10+detIdStub.iZ()+abs(detIdStub.iSide()-2)*7;

      if (layer>10 && layer<=17) disk=(layer-10)%8;
      if (layer>17 && layer<=24) disk=(layer-17)%8;
      if (disk>=5) lad_cor = (disk-4)%4;

      ladder = detIdStub.iRing()-1+lad_cor;
      module = detIdStub.iPhi()-1;
    }    
    
    module = CMSPatternLayer::getModuleCode(layer,module);
    if(module<0)// the stub is on the third Z position on the other side of the tracker -> out of range
      continue;

    ladder = CMSPatternLayer::getLadderCode(layer, ladder);
    
    int strip  =  mp0.x();
    int tp     = -1;
    float eta  = 0;
    float phi0 = 0;
    float spt  = 0;
    float x    = posStub.x();
    float y    = posStub.y();
    float z    = posStub.z();
    float x0   = 0.;
    float y0   = 0.;
    float z0   = 0.;
    float ip   = sqrt(x0*x0+y0*y0);

    Hit* h = new Hit(layer,ladder, module, segment, strip, j, tp, spt, ip, eta, phi0, x, y, z, x0, y0, z0);
    m_hits.push_back(h);
  } /// End of loop over input stubs


  /// STEP 2
  /// PAssing the superstrips into the AM chip

  //  std::cout << "AM chip processing" << std::endl;

  vector<Sector*> patternsSectors = m_pf->find(m_hits); // AM PR is done here....


  /// STEP 3
  /// Collect the info and store the track seed stuff

  //  std::cout << "AM chip processing" << std::endl;

  vector<Hit*> hits; 
  for(unsigned int i=0;i<patternsSectors.size();i++)
  {      
    vector<GradedPattern*> pl = patternsSectors[i]->getPatternTree()->getLDPatterns();

    if (pl.size()==0) continue; // No patterns

    int secID = patternsSectors[i]->getOfficialID();

    //    std::cout<<"Found "<<pl.size()<<" patterns in sector " << secID<<std::endl;

    //delete the GradedPattern objects
    for(unsigned j=0;j<pl.size();j++) 
    {
      hits.clear();
      hits = pl[j]->getHits();

      /// Create the Seed in the form of a Track and store it in the output
      std::vector< edm::Ptr< TTStub< Ref_PixelDigi_ > > > tempVec;
     
      for(unsigned k=0;k<hits.size();k++) 
      {
        edm::Ptr< TTStub< Ref_PixelDigi_ > > tempStubPtr( input, hits[k]->getID());
        tempVec.push_back( tempStubPtr );
      }

      TTTrack< Ref_PixelDigi_ > tempTrack( tempVec );
      tempTrack.setSector( secID );
      output.push_back( tempTrack );

      delete pl[j];
    }

    //delete the Sectors
    delete patternsSectors[i];
  }
}
