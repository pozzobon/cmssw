import FWCore.ParameterSet.Config as cms

### Needed to access DTConfigManagerRcd and by DTTrig
#from L1TriggerConfig.DTTPGConfigProducers.L1DTTPGConfig_cff import *

DTTPG = cms.EDProducer( "DTTPG",
  #OnlyRPhi = cms.untracked.bool(True),
  #SingleChamberTest = cms.untracked.bool(True),
  FirstBX = cms.untracked.int32(0),
  LastBX = cms.untracked.int32(3564)
)
