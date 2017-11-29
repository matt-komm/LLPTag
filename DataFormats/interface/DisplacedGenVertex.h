#ifndef CMGTools_RA1_DisplacedGenVertex_h
#define CMGTools_RA1_DisplacedGenVertex_h

#include "DataFormats/Common/interface/RefProd.h"

#include "DataFormats/Common/interface/Ref.h"
#include "DataFormats/Common/interface/RefVector.h"

#include "DataFormats/Common/interface/Ptr.h"
#include "DataFormats/Common/interface/PtrVector.h"

#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "DataFormats/JetReco/interface/GenJet.h"


#include <vector>


// a simple class
struct DisplacedGenVertex
{
    bool isHardInteraction;
    double sharedMassFraction; //fraction of invariant LLP mass carried by its first daughters in vertex
    reco::Candidate::Point vertex;
    edm::Ptr<DisplacedGenVertex> motherVertex;
    edm::RefVector<std::vector<DisplacedGenVertex>> daughterVertices;
    
    edm::PtrVector<reco::GenParticle> genParticles;
    edm::Ptr<reco::GenParticle> motherLongLivedParticle;
    
    edm::PtrVector<reco::GenJet> genJets;
    
    DisplacedGenVertex():
        isHardInteraction(false),
        sharedMassFraction(0),
        vertex(0,0,0)
    {
    }
    
    double d3d() const;
    double dx() const;
    double dy() const;
    double dz() const;
    double dxy() const;
};

typedef std::vector<DisplacedGenVertex> DisplacedGenVertexCollection;

#endif
