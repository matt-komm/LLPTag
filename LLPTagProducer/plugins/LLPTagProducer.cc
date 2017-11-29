// system include files
#include <memory>
#include <vector>
#include <unordered_map>
// user include files

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Utilities/interface/Exception.h"

#include "LLPTag/DataFormats/interface/DisplacedGenVertex.h"

#include "SimGeneral/HepPDTRecord/interface/ParticleDataTable.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "DataFormats/Math/interface/deltaR.h"


class LLPTagProducer:
    public edm::stream::EDProducer<>
    
{
    private:    
    
        
 

    public:
        explicit LLPTagProducer(const edm::ParameterSet&);
        ~LLPTagProducer();

        void produce(edm::Event& iEvent, const edm::EventSetup& iSetup);
        static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

};


//
// constructors and destructor

//
LLPTagProducer::LLPTagProducer(const edm::ParameterSet& iConfig)
{
}


LLPTagProducer::~LLPTagProducer()
{
}


// ------------ method called to produce the data  ------------
void
LLPTagProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
}



// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
LLPTagProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}



//define this as a plug-in
DEFINE_FWK_MODULE(LLPTagProducer);
