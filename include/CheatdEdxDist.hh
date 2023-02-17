#ifndef CheatdEdxDist_h
#define CheatdEdxDist_h 1
#include <iomanip>
#include <EVENT/LCRelation.h>
#include "marlin/Processor.h"
#include <EVENT/LCCollection.h>
#include <EVENT/MCParticle.h>
#include <EVENT/ReconstructedParticle.h>
#include <IMPL/LCCollectionVec.h>
#include <IMPL/TrackImpl.h>
#include <IMPL/TrackStateImpl.h>
#include <EVENT/Track.h>
#include <EVENT/Vertex.h>
#include <IMPL/ReconstructedParticleImpl.h>
#include <IMPL/LCRelationImpl.h>
#include <UTIL/LCRelationNavigator.h>
#include <UTIL/PIDHandler.h>
#include <cmath>
#include <string>
#include <vector>

#include <TF1.h>
#include <TRandom1.h>
#include <TFile.h>
#include <TTree.h>
#include "MathOperator.hh"
// #include <TLorentzVector.h>
using namespace lcio;
using namespace marlin;

class CheatdEdxDist : public Processor
{

public:
  virtual Processor *newProcessor() { return new CheatdEdxDist; }

  CheatdEdxDist();
  virtual ~CheatdEdxDist();

  /** Called at the begin of the job before anything is read.
   */
  virtual void init();

  /** Called for every run.
   */
  virtual void processRunHeader(LCRunHeader *run);

  /** Called for every event - the working horse.
   */
  virtual void processEvent(LCEvent *evt);

  virtual void check(LCEvent *evt);

  virtual void end();

private:
  MCParticle *getMCParticle(EVENT::ReconstructedParticle *particle, LCCollection *myRelCollection);

  std::vector<std::string> _methodstorun;
  std::vector<std::string> _dEdxNames;
  std::string _methodstorun_version="_v3";
  std::string _versionPID_nominal="";


  std::string _MCColName;
  std::string _MCRelColName;
  std::string _colName;

  std::vector<float> _dedx_dist_kaon_mean={0.0260616,0.0028051,-0.0276483,-0.0354468,-0.0328854,-0.0299818,-0.0394224,-0.074456,-0.112405,-0.0592645};
  std::vector<float> _dedx_dist_kaon_rms={1.19158,1.18987,1.18354,1.17574,1.16956,1.16489,1.16028,1.16037,1.08117,0.856091};
  std::vector<float> _dedx_dist_pion_mean={3.90534,3.88268,3.85414,3.85134,3.87144,3.9142,3.96996,4.04427,3.65627,2.68418};
  std::vector<float> _dedx_dist_pion_rms={1.02647,1.02467,1.01896,1.01646,1.01421,1.01079,1.01255,1.02192,1.10763,1.23959};
  std::vector<float> _dedx_dist_proton_mean={-2.22942,-2.24991,-2.27241,-2.2816,-2.29309,-2.31275,-2.35251,-2.44654,-2.27726,-1.618};
  std::vector<float> _dedx_dist_proton_rms={1.35219,1.35673,1.35233,1.34314,1.33618,1.32997,1.32287,1.32804,1.22818,0.990744};
  std::vector<float> _dedx_dist_muon_mean={4.39099,4.36903,4.33255,4.32964,4.35646,4.39512,4.45686,4.53304,4.08812,2.98312};
  std::vector<float> _dedx_dist_muon_rms={1.06788,1.06773,1.06846,1.06721,1.06307,1.06725,1.07258,1.09079,1.20165,1.37217};
  std::vector<float> _dedx_dist_electron_mean={7.63539,7.63139,7.5973,7.59836,7.64307,7.71311,7.81468,7.95942,7.22736,5.14705};
  std::vector<float> _dedx_dist_electron_rms={1.5517,1.49618,1.5434,1.55364,1.51213,1.5219,1.56808,1.65505,1.85577,2.34292};
  std::vector<float> _dedx_dist_others_mean={-3.19445,-3.11617,-3.16022,-3.20967,-3.20955,-3.2378,-3.28648,-3.38983,-3.1024,-2.15415};
  std::vector<float> _dedx_dist_others_rms={1.3387,1.3472,1.37304,1.36929,1.35874,1.34594,1.33698,1.30841,1.224,1.05959};


  std::vector<std::vector<float>> mean;
  std::vector<std::vector<float>> rms;

  float improvement=0.;// from 0.1 improvement means a 10% improvement

};

#endif
