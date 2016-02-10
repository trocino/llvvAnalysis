// -*- C++ -*-
//
// Package:    GenHTFilter
// Class:      GenHTFilter
// 
/**\class GenHTFilter GenHTFilter.cc FinalStateAnalysis/GenHTFilter/src/GenHTFilter.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Nick Smith
//         Created:  Fri Jun 12 02:15:46 CDT 2015
// $Id$
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Math/interface/deltaR.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"

//
// class declaration
//

class GenHTFilter : public edm::EDFilter {
   public:
      explicit GenHTFilter(const edm::ParameterSet&);
      ~GenHTFilter();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void beginJob() ;
      virtual bool filter(edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
      
      virtual bool beginRun(edm::Run&, edm::EventSetup const&);
      virtual bool endRun(edm::Run&, edm::EventSetup const&);
      virtual bool beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);
      virtual bool endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&);

      // ----------member data ---------------------------

   edm::EDGetTokenT<LHEEventProduct> lheEventProductToken_;
   double htMin_, htMax_;
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
GenHTFilter::GenHTFilter(const edm::ParameterSet& iConfig) :
      lheEventProductToken_(consumes<LHEEventProduct>(iConfig.getUntrackedParameter<edm::InputTag>("lheEventProduct", edm::InputTag("externalLHEProducer")))),
      htMin_(iConfig.getParameter<double>("htMin")),
      htMax_(iConfig.getParameter<double>("htMax"))
{

}


GenHTFilter::~GenHTFilter()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called on each new Event  ------------
bool
GenHTFilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

   Handle<LHEEventProduct> lheProduct;
   iEvent.getByToken(lheEventProductToken_, lheProduct);

   const lhef::HEPEUP& lheeventinfo = lheProduct->hepeup();

   float sumpt=0;
   for (int i = 0; i < lheeventinfo.NUP ; ++i) {
         if (lheeventinfo.ISTUP[i] !=1 ||((abs(lheeventinfo.IDUP[i])>5&&lheeventinfo.IDUP[i]!=21) ))  continue;
         double px=lheeventinfo.PUP.at(i)[0];
         double py=lheeventinfo.PUP.at(i)[1];
         double pt=sqrt(px*px+py*py);
         sumpt+=pt;
   }

   if ( sumpt < htMin_ || sumpt > htMax_ ) return false;
   return true;
}

// ------------ method called once each job just before starting event loop  ------------
void 
GenHTFilter::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
GenHTFilter::endJob() {
}

// ------------ method called when starting to processes a run  ------------
bool 
GenHTFilter::beginRun(edm::Run&, edm::EventSetup const&)
{ 
  return true;
}

// ------------ method called when ending the processing of a run  ------------
bool 
GenHTFilter::endRun(edm::Run&, edm::EventSetup const&)
{
  return true;
}

// ------------ method called when starting to processes a luminosity block  ------------
bool 
GenHTFilter::beginLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
  return true;
}

// ------------ method called when ending the processing of a luminosity block  ------------
bool 
GenHTFilter::endLuminosityBlock(edm::LuminosityBlock&, edm::EventSetup const&)
{
  return true;
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
GenHTFilter::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}
//define this as a plug-in
DEFINE_FWK_MODULE(GenHTFilter);
