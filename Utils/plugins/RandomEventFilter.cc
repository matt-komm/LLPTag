#include <vector>
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/stream/EDFilter.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

class RandomEventFilter: 
    public edm::stream::EDFilter<>
{
    public:
        explicit RandomEventFilter(const edm::ParameterSet &);
        ~RandomEventFilter();

    private:
        virtual bool filter(edm::Event& event, const edm::EventSetup& setup);
        const double percentage_;
        const bool invert_;
};


RandomEventFilter::RandomEventFilter(const edm::ParameterSet& config):
    percentage_(config.getParameter<double>("percentage")),
    invert_(config.getParameter<bool>("invert"))
{
}

RandomEventFilter::~RandomEventFilter()
{
}

bool RandomEventFilter::filter(edm::Event& event, const edm::EventSetup& setup) 
{
    size_t hash =  std::hash<unsigned long long>{}(event.id().event()) ^ ( std::hash<unsigned long long>{}(event.id().luminosityBlock()) << 1);
    if (hash%100000<(percentage_*100000.))
    {
        return !invert_;
    }
    return invert_;
}

DEFINE_FWK_MODULE(RandomEventFilter);
