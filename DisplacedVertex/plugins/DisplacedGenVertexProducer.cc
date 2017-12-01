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

#include "XTag/DataFormats/interface/DisplacedGenVertex.h"

#include "SimGeneral/HepPDTRecord/interface/ParticleDataTable.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "DataFormats/Math/interface/deltaR.h"

using xtag::DisplacedGenVertex;
using xtag::DisplacedGenVertexCollection;

class DisplacedGenVertexProducer:
    public edm::stream::EDProducer<>
    
{
    private:    
        static double distance(const reco::Candidate::Point& p1, const reco::Candidate::Point& p2)
        {
            return std::sqrt((p1-p2).mag2());
        }
        
        static int getHadronFlavor(const reco::GenParticle& genParticle)
        {
            int absPdgId = std::abs(genParticle.pdgId());
            int nq3 = (absPdgId/     10)%10; //quark content
            int nq2 = (absPdgId/    100)%10; //quark content
            int nq1 = (absPdgId/   1000)%10; //quark content
            int nL  = (absPdgId/  10000)%10;
            int n   = (absPdgId/1000000)%10;
            return std::max({nq1,nq2,nq3})+n*10000+(n>0 and nL==9)*100;
        }
        
        static bool displacedDecay(const reco::GenParticle& genParticle)
        {
            for (unsigned int idaughter = 0; idaughter<genParticle.numberOfDaughters(); ++idaughter)
            {
                if (distance(genParticle.daughter(idaughter)->vertex(),genParticle.vertex())>10e-10)
                {
                    return true;
                }
            }
            return false;
        }
        
    
        edm::InputTag _genParticleInputTag;
        edm::EDGetTokenT<edm::View<reco::GenParticle>> _genParticleToken;
        
        edm::InputTag _genJetInputTag;
        edm::EDGetTokenT<edm::View<reco::GenJet>> _genJetToken;
        
        
        virtual void produce(edm::Event& iEvent, const edm::EventSetup& iSetup) override;
            
            
        template<typename HANDLE, typename TYPE>
        void writeValueMap(edm::Event &out, const HANDLE& handle, const std::vector<TYPE> values, const std::string &name) const 
        {
             typedef edm::ValueMap<TYPE> Map;
             std::unique_ptr<Map> map(new Map());
             typename Map::Filler filler(*map);
             filler.insert(handle, values.begin(), values.end());
             filler.fill();
             out.put(std::move(map), name);
        }
 

    public:
        explicit DisplacedGenVertexProducer(const edm::ParameterSet&);
        ~DisplacedGenVertexProducer();

        static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

};


//
// constructors and destructor

//
DisplacedGenVertexProducer::DisplacedGenVertexProducer(const edm::ParameterSet& iConfig):
    _genParticleInputTag(iConfig.getParameter<edm::InputTag>("srcGenParticles")),
    _genParticleToken(consumes<edm::View<reco::GenParticle>>(_genParticleInputTag)),
    _genJetInputTag(iConfig.getParameter<edm::InputTag>("srcGenJets")),
    _genJetToken(consumes<edm::View<reco::GenJet>>(_genJetInputTag))
{
    //produces<edm::ValueMap<unsigned int>>("vertexGroupIndex");
    produces<std::vector<DisplacedGenVertex>>();
    /*
    produces<edm::ValueMap<double>>("distance");
    produces<edm::ValueMap<double>>("ctau");
    */
}


DisplacedGenVertexProducer::~DisplacedGenVertexProducer()
{
}


// ------------ method called to produce the data  ------------
void
DisplacedGenVertexProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    edm::Handle<edm::View<reco::GenParticle>> genParticleCollection;
    iEvent.getByToken(_genParticleToken, genParticleCollection);
    
    edm::Handle<edm::View<reco::GenJet>> genJetCollection;
    iEvent.getByToken(_genJetToken, genJetCollection);
    
    std::unordered_map<size_t,unsigned int> genParticleToVertexGroupMap;
    std::shared_ptr<reco::Candidate::Point> hardInteractionVertex(nullptr);
    
    std::unique_ptr<DisplacedGenVertexCollection> displacedGenVertices(new DisplacedGenVertexCollection());
    edm::RefProd<DisplacedGenVertexCollection> displacedGenVertices_refProd = iEvent.getRefBeforePut<DisplacedGenVertexCollection>();
    
    for (unsigned int igenParticle = 0; igenParticle < genParticleCollection->size(); ++igenParticle)
    {
        const reco::GenParticle& genParticle = genParticleCollection->at(igenParticle);
        if (genParticle.status()==22)
        {
            if (!hardInteractionVertex)
            {
                hardInteractionVertex.reset(new reco::Candidate::Point(genParticle.vertex()));
            }
            else if (distance(*hardInteractionVertex,genParticle.vertex())>10e-10)
            {
                throw cms::Exception("DisplacedGenVertexProducer: multiple hard interaction vertices found!");
            }
        }
        
        //group particles by vertex position
        bool inserted = false;
        for (unsigned int ivertex = 0; ivertex<displacedGenVertices->size(); ++ivertex)
        {
            DisplacedGenVertex& displacedGenVertex = displacedGenVertices->at(ivertex);
            if (distance(displacedGenVertex.vertex,genParticle.vertex())<10e-10)
            {
                displacedGenVertex.genParticles.push_back(genParticleCollection->ptrAt(igenParticle));
                genParticleToVertexGroupMap[(size_t)&genParticle]=ivertex;
                inserted=true;
                break;
            }
        }
        if (not inserted)
        {
            DisplacedGenVertex displacedGenVertex;
            displacedGenVertex.vertex = genParticle.vertex();
            displacedGenVertex.genParticles.push_back(genParticleCollection->ptrAt(igenParticle));
            displacedGenVertices->push_back(displacedGenVertex);
            genParticleToVertexGroupMap[(size_t)&genParticle]=displacedGenVertices->size()-1;
        }
    }
    
    for (unsigned int ivertex = 0; ivertex<displacedGenVertices->size(); ++ivertex)
    {
        if (hardInteractionVertex and distance(displacedGenVertices->at(ivertex).vertex,*hardInteractionVertex)<10e-10)
        {
            displacedGenVertices->at(ivertex).isHardInteraction=true;
            break;
        }
    }
    
    //use longlived particles to link vertex groups
    for (unsigned int igenParticle = 0; igenParticle < genParticleCollection->size(); ++igenParticle)
    {
        const reco::GenParticle& genParticle = genParticleCollection->at(igenParticle);
        //if (not ((isGluino(genParticle) or isGluinoHadron(genParticle) or isHadron(genParticle,5) or isHadron(genParticle,4)) and genParticle.isLastCopy() and displacedDecay(genParticle)))
        if (genParticle.mass()<=0)
        {
            continue;
        }
        if (not (genParticle.isLastCopy() and displacedDecay(genParticle)))
        {
            continue;
        }

        unsigned int originVertexGroupIndex = genParticleToVertexGroupMap.at((size_t)&genParticle);
        //book keep summed p4 per vertex group
        std::unordered_map<unsigned int,reco::Candidate::LorentzVector> momentumDistribution;

        for (unsigned int idaughter = 0; idaughter<genParticle.numberOfDaughters(); ++idaughter)
        {
            const reco::Candidate* daughter = genParticle.daughter(idaughter);
            unsigned int daughterVertexGroupIndex = genParticleToVertexGroupMap.at((size_t)daughter);
            if (originVertexGroupIndex!=daughterVertexGroupIndex)
            {
                momentumDistribution[daughterVertexGroupIndex]+=daughter->p4();
            }
        }
        if (genParticle.numberOfDaughters()==0) throw cms::Exception("DisplacedGenVertexProducer: Particle has no daughters in vertex groups");
        //find vertex which shares most of the invariant mass with the long lived particle
        double maxMassRatio = -1;
        int decayVertexIndex = -1;
        for (auto idMomentumPair: momentumDistribution)
        {
            double massRatio = idMomentumPair.second.mass()/genParticle.mass()+1e-8;
            if (massRatio>maxMassRatio)
            {
                maxMassRatio=massRatio;
                decayVertexIndex = idMomentumPair.first;
            }
        }
        if (decayVertexIndex<0) throw cms::Exception("DisplacedGenVertexProducer: A long lived particle should always connect two vertices");
        //make link
        displacedGenVertices->at(originVertexGroupIndex).sharedMassFraction = maxMassRatio;
        displacedGenVertices->at(originVertexGroupIndex).daughterVertices.push_back(edm::Ref<DisplacedGenVertexCollection>(displacedGenVertices_refProd,decayVertexIndex));
        edm::Ref<DisplacedGenVertexCollection> motherRef(displacedGenVertices_refProd,originVertexGroupIndex);
        displacedGenVertices->at(decayVertexIndex).motherVertex = std::move(edm::Ptr<DisplacedGenVertex>(motherRef.id(),motherRef.key(),motherRef.productGetter()));
        //store long lived particle in daughter vertex
        displacedGenVertices->at(decayVertexIndex).motherLongLivedParticle = std::move(edm::Ptr<reco::GenParticle>(genParticleCollection,igenParticle));
    }
    
    
    //calculate overlap of gen particles in gen jets with displaced vertices and associate them
    //IMPORTANT: the following will only work if the GenJets were constructed from the SAME GenParticle collection as the DisplacedGenVertices!!!
    //Note: cannot use ghost tagging since jets may not be pointing along long lived particle direction

    for (unsigned int ijet = 0; ijet < genJetCollection->size(); ++ijet)
    {
        std::vector<unsigned int> matchedVerticesIndices;
        for (const reco::GenParticle* genParticle: genJetCollection->at(ijet).getGenConstituents())
        {
            auto foundGenParticleIt = genParticleToVertexGroupMap.find((size_t)genParticle); 
            if (foundGenParticleIt!=genParticleToVertexGroupMap.end())
            {
                if (std::find(matchedVerticesIndices.begin(),matchedVerticesIndices.end(),foundGenParticleIt->second)==matchedVerticesIndices.end())
                {
                    matchedVerticesIndices.push_back(foundGenParticleIt->second);
                }
            }
        }
        if (matchedVerticesIndices.size()==0)
        {
            continue;
        }
        else if (matchedVerticesIndices.size()==1)
        {
            displacedGenVertices->at(matchedVerticesIndices.front()).genJets.push_back(genJetCollection->ptrAt(ijet));
        }
        else if (matchedVerticesIndices.size()>=2)
        {
            auto displacedVertexIndexIt = std::max_element(matchedVerticesIndices.begin(),matchedVerticesIndices.end(),
                [&displacedGenVertices](const unsigned int& iv1, const unsigned int& iv2)
                {

                    const DisplacedGenVertex& v1 = displacedGenVertices->at(iv1);
                    const DisplacedGenVertex& v2 = displacedGenVertices->at(iv2);
                    
                    if (v1.motherLongLivedParticle.isNull() and v2.motherLongLivedParticle.isNull())
                    {
                        return true; //arbitrary order
                    }
                    else if (v1.motherLongLivedParticle.isNonnull() and v2.motherLongLivedParticle.isNull())
                    {
                        return true;
                    }
                    else if (v1.motherLongLivedParticle.isNull() and v2.motherLongLivedParticle.isNonnull())
                    {
                        return false;
                    }
                    
                    return getHadronFlavor(*v1.motherLongLivedParticle)>getHadronFlavor(*v2.motherLongLivedParticle);
                }
            );
            
            displacedGenVertices->at(*displacedVertexIndexIt).genJets.push_back(genJetCollection->ptrAt(ijet)); //no cross cleaning -> multiple jets can originate from one vertex
            
        }
    }
    
    //TODO: collaps vertices if llp is stable particle e.g. proton/pion/kaon
    
    iEvent.put(std::move(displacedGenVertices));
}



// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
DisplacedGenVertexProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}



//define this as a plug-in
DEFINE_FWK_MODULE(DisplacedGenVertexProducer);
