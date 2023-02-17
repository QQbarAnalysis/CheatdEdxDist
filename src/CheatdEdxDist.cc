#include "CheatdEdxDist.hh"

// ----- include for verbosity dependent logging ---------
// #include "marlin/VerbosityLevels.h"
// #include "marlin/StringParameters.h"
// #define SLM streamlog_out(MESSAGE)

// using namespace std;
// using namespace lcio ;
// using namespace marlin ;
using EVENT::LCCollection;
using EVENT::MCParticle;
using EVENT::ReconstructedParticle;
using EVENT::Track;
using EVENT::Vertex;
using HiddenVAnalysis::MathOperator;
using IMPL::LCRelationImpl;
using IMPL::ReconstructedParticleImpl;
using IMPL::TrackImpl;
using IMPL::TrackStateImpl;
using std::string;
using std::vector;
using UTIL::LCRelationNavigator;

CheatdEdxDist aCheatdEdxDist;

CheatdEdxDist::CheatdEdxDist() : Processor("CheatdEdxDist")
{

	// modify processor description
	_description = "";

	registerProcessorParameter("_dedx_dist_kaon_mean",
							   "Mean of Kaon-dEdx dist with kaon particles (from 0-1 in costheta)",
							   _dedx_dist_kaon_mean,
							   _dedx_dist_kaon_mean);
	registerProcessorParameter("_dedx_dist_kaon_rms",
							   "RMS of Kaon-dEdx dist with kaon particles (from 0-1 in costheta)",
							   _dedx_dist_kaon_rms,
							   _dedx_dist_kaon_rms);

	registerProcessorParameter("_dedx_dist_proton_mean",
							   "Mean of Kaon-dEdx dist with proton particles (from 0-1 in costheta)",
							   _dedx_dist_proton_mean,
							   _dedx_dist_proton_mean);
	registerProcessorParameter("_dedx_dist_proton_rms",
							   "RMS of Kaon-dEdx dist with proton particles (from 0-1 in costheta)",
							   _dedx_dist_proton_rms,
							   _dedx_dist_proton_rms);

	registerProcessorParameter("_dedx_dist_pion_mean",
							   "Mean of Kaon-dEdx dist with pion particles (from 0-1 in costheta)",
							   _dedx_dist_pion_mean,
							   _dedx_dist_pion_mean);
	registerProcessorParameter("_dedx_dist_pion_rms",
							   "RMS of Kaon-dEdx dist with pion particles (from 0-1 in costheta)",
							   _dedx_dist_pion_rms,
							   _dedx_dist_pion_rms);

	registerProcessorParameter("_dedx_dist_muon_mean",
							   "Mean of Kaon-dEdx dist with muon particles (from 0-1 in costheta)",
							   _dedx_dist_muon_mean,
							   _dedx_dist_muon_mean);
	registerProcessorParameter("_dedx_dist_muon_rms",
							   "RMS of Kaon-dEdx dist with muon particles (from 0-1 in costheta)",
							   _dedx_dist_muon_rms,
							   _dedx_dist_muon_rms);

	registerProcessorParameter("_dedx_dist_electron_mean",
							   "Mean of Kaon-dEdx dist with electron particles (from 0-1 in costheta)",
							   _dedx_dist_electron_mean,
							   _dedx_dist_electron_mean);
	registerProcessorParameter("_dedx_dist_electron_rms",
							   "RMS of Kaon-dEdx dist with electron particles (from 0-1 in costheta)",
							   _dedx_dist_electron_rms,
							   _dedx_dist_electron_rms);

	registerProcessorParameter("_dedx_dist_others_mean",
							   "Mean of Kaon-dEdx dist with others particles (from 0-1 in costheta)",
							   _dedx_dist_others_mean,
							   _dedx_dist_others_mean);
	registerProcessorParameter("_dedx_dist_others_rms",
							   "RMS of Kaon-dEdx dist with others particles (from 0-1 in costheta)",
							   _dedx_dist_others_rms,
							   _dedx_dist_others_rms);

	// input collections			std::string _methodstorun_version="_v3";

	registerProcessorParameter("improvement",
							   "improvement of the dEdx calculation, 0.1 means 10% improvement",
							   improvement,
							   improvement);
	registerProcessorParameter("_methodstorun_version",
							   "tag for the new PID algorithm dEdxPID_tag",
							   _methodstorun_version,
							   _methodstorun_version);

	registerProcessorParameter("_versionPID_nominal",
							   "tag of the default new PID algorithm dEdxPID_tag",
							   _versionPID_nominal,
							   _versionPID_nominal);
	registerInputCollection(LCIO::RECONSTRUCTEDPARTICLE,
							"PFOCollection",
							"PFO collection name",
							_colName,
							string("PandoraPFOs"));
	registerInputCollection(LCIO::MCPARTICLE,
							"MCCollectionName",
							"Name of the MC collection",
							_MCColName,
							std::string("MCParticlesSkimmed"));
	registerInputCollection(LCIO::MCPARTICLE,
							"MCRelCollectionName",
							"Name of the MC-PFO relation collection",
							_MCRelColName,
							std::string("RecoMCTruthLink"));
}

CheatdEdxDist::~CheatdEdxDist() {}

void CheatdEdxDist::init()
{

	_methodstorun.push_back("dEdxPID");
	_dEdxNames.push_back("kaon_dEdxdistance");

	mean.push_back(_dedx_dist_kaon_mean);
	mean.push_back(_dedx_dist_pion_mean);
	mean.push_back(_dedx_dist_proton_mean);
	mean.push_back(_dedx_dist_muon_mean);
	mean.push_back(_dedx_dist_electron_mean);
	mean.push_back(_dedx_dist_others_mean);

	rms.push_back(_dedx_dist_kaon_rms);
	rms.push_back(_dedx_dist_pion_rms);
	rms.push_back(_dedx_dist_proton_rms);
	rms.push_back(_dedx_dist_muon_rms);
	rms.push_back(_dedx_dist_electron_rms);
	rms.push_back(_dedx_dist_others_rms);

	printParameters();
}

// *****************************************************************************************
// PFO to MC relations
MCParticle *CheatdEdxDist::getMCParticle(EVENT::ReconstructedParticle *particle, LCCollection *myRelCollection)
{

	LCRelationNavigator navigator(myRelCollection);
	streamlog_out(DEBUG) << "getPDG: using Reconstructed collection \n";

	vector<LCObject *> obj = navigator.getRelatedToObjects(particle);
	vector<float> weights = navigator.getRelatedToWeights(particle);
	MCParticle *winner = NULL;
	MCParticle *winner2 = NULL;

	float maxweight = 0.50;
	for (unsigned int i = 0; i < obj.size(); i++)
	{
		if (particle->getTracks().size() != 0)
		{
			weights[i] = fmod(weights[i], 10000) / 1000.;
			streamlog_out(DEBUG) << "getPDG : PFO is a TRACK \n";
		}
		else
		{
			streamlog_out(DEBUG) << "getPDG : PFO is a not a TRACK \n";
		}
		// considerations on the track weigth from Mikael:
		// The track-part of the weight of a PFO should just be a copy of the weight in the
		// corresponding track-object (in permil) so (pfowgt%10000)/1000.0 should
		// be the same as the track weight. Be aware of the navigators in two
		// direction (RecoMCTruthLink and MCTruthRecoLink): They link the same
		// objects, but the weights are different: One tells you which fraction of
		// the PFO observation comes from which MCParticle, the other tells you
		// which fraction of all signal produced of the MCParticle goes into the
		// PFO in question....

		winner2 = dynamic_cast<MCParticle *>(obj[i]);

		streamlog_out(DEBUG) << "getMCParticle " << i << " max weight " << maxweight << " weight " << weights[i] << " MCPDG " << winner2->getPDG() << " mom " << MathOperator::getModule(winner2->getMomentum()) << " px: " << winner2->getMomentum()[0] << " py: " << winner2->getMomentum()[1] << " En\
d point (in m), x, y, z " << winner2->getEndpoint()[0]
							 << " " << winner2->getEndpoint()[1] << " " << winner2->getEndpoint()[2] << "\n";
		if (weights[i] > maxweight)
		{
			winner = dynamic_cast<MCParticle *>(obj[i]);
			maxweight = weights[i];
		}
	}
	if (!winner)
		return NULL;
	streamlog_out(DEBUG) << "end getPDG " << maxweight << " " << winner->getPDG() << "\n";
	return winner;
}

void CheatdEdxDist::processRunHeader(LCRunHeader *run)
{
}

void CheatdEdxDist::processEvent(LCEvent *evt)
{

	try
	{
		LCCollection *mccol = evt->getCollection(_MCColName);
		LCCollection *pfocol = evt->getCollection(_colName);
		LCCollection *myRelCollection = evt->getCollection(_MCRelColName);
		PIDHandler pidh(evt->getCollection(_colName));
		pidh.addAlgorithm(_methodstorun.at(0) + _methodstorun_version, _dEdxNames);
		int algoID = pidh.getAlgorithmID(_methodstorun.at(0) + _methodstorun_version);

		for (int i = 0; i < pfocol->getNumberOfElements(); i++)
		{

			ReconstructedParticle *pfo_rp = dynamic_cast<ReconstructedParticle *>(pfocol->getElementAt(i));
			MCParticle *mcpart = getMCParticle(pfo_rp, myRelCollection);

			if (mcpart == NULL)
			{
				std::vector<float> valores;
				valores.push_back(-99999);
				pidh.setParticleID(pfo_rp, 0, -99999, -99999, algoID, valores);
				continue;
			}
			int pdg_mc = mcpart->getPDG();
			float dedx_nominal = 0;
			try
			{
				int algoID_nominal = pidh.getAlgorithmID("dEdxPID" + _versionPID_nominal);
				const ParticleID &pid_nominal = pidh.getParticleID(pfo_rp, algoID_nominal);
				int index_dedxdist_nominal = pidh.getParameterIndex(algoID_nominal, "kaon_dEdxdistance");
				std::vector<float> params_nominal = pid_nominal.getParameters();
				for(int ip=0; ip<params_nominal.size();ip++) std::cout<<ip<<" "<<index_dedxdist_nominal<<std::endl;
				if(params_nominal.size()>0) 	dedx_nominal = params_nominal.at(index_dedxdist_nominal);
			}
			catch (DataNotAvailableException &e)
			{
				streamlog_out(DEBUG) << "dEdx distance of nominal dEdxPID" + _versionPID_nominal + " doesn't EXIST \n";
				streamlog_out(DEBUG) << e.what();
			}

			if (dedx_nominal == 0)
				continue;


			float costheta1 = -2.0;
			vector<float> d1 = MathOperator::getDirection(mcpart->getMomentum());

			costheta1 = std::cos(MathOperator::getAngles(d1)[1]);
			int icostheta = int(fabs(costheta1) * 10);

			// cases
			int id[5] = {321, 211, 2212, 13, 11};

			if (pfo_rp->getTracks().size() == 1)
			{

				TRandom *r1 = new TRandom1();
				float dedxdist = -9999;
				bool hadronfound = false;
				for (int ir = 0; ir < 5; ir++)
				{
					if (fabs(pdg_mc) == id[ir])
					{
						dedxdist = r1->Gaus(mean.at(ir).at(icostheta), rms.at(ir).at(icostheta)*(1.0-improvement));
						hadronfound = true;
						continue;
					}
				}
				if (hadronfound == false)
					dedxdist = r1->Gaus(mean.at(5).at(icostheta), rms.at(5).at(icostheta)*(1.0-improvement));

				std::vector<float> valores;
				valores.push_back(dedxdist);
				pidh.setParticleID(pfo_rp, 0, pdg_mc, 1, algoID, valores);
			}
			else
			{
				std::vector<float> valores;
				valores.push_back(-9999);
				pidh.setParticleID(pfo_rp, 0, pdg_mc, 0, algoID, valores);
			}
		}
	}
	catch (DataNotAvailableException &e)
	{
		streamlog_out(DEBUG) << "Whoops!....\n";
		streamlog_out(DEBUG) << e.what();
	}
}

void CheatdEdxDist::check(LCEvent *evt)
{
	// nothing to check here - could be used to fill checkplots in reconstruction processor
}

void CheatdEdxDist::end()
{
}
