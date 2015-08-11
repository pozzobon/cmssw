#ifndef Geometry_GEMGeometry_ME0Chamber_h
#define Geometry_GEMGeometry_ME0Chamber_h

/** \class ME0Chamber
 *
 *  Model of a ME0 chamber.
 *   
 *  A chamber is a GeomDet.
 *
 *  \author S. Dildick
 *  \edited by D. Nash
 */

#include "Geometry/CommonDetUnit/interface/GeomDet.h"
#include "DataFormats/MuonDetId/interface/ME0DetId.h"

class ME0EtaPartition;

class ME0Chamber : public GeomDet {
public:
  /// Constructor
  ME0Chamber(ME0DetId id, const ReferenceCountingPointer<BoundPlane>& plane);

  /// Destructor
  virtual ~ME0Chamber();

  /// Return the ME0DetId of this chamber
  ME0DetId id() const;

  // Which subdetector
  virtual SubDetector subDetector() const {return GeomDetEnumerators::ME0;}

  /// equal if the id is the same
  bool operator==(const ME0Chamber& ch) const;

  /// Add EtaPartition to the chamber which takes ownership
  void add(ME0EtaPartition* roll);

  /// Return the rolls in the chamber
  virtual std::vector<const GeomDet*> components() const;

  /// Return the sub-component (roll) with a given id in this chamber
  virtual const GeomDet* component(DetId id) const;

  /// Return the eta partition corresponding to the given id 
  const ME0EtaPartition* etaPartition(ME0DetId id) const;

  const ME0EtaPartition* etaPartition(int isl) const;
  
  /// Return the eta partitions
  const std::vector<const ME0EtaPartition*>& etaPartitions() const;

  /// Retunr numbers of eta partitions
  int nEtaPartitions() const;

private:

  ME0DetId detId_;

  // vector of eta partitions for a chamber
  std::vector<const ME0EtaPartition*> etaPartitions_;

};
#endif
