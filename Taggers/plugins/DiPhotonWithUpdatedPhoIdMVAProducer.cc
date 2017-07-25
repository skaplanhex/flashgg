#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "flashgg/MicroAOD/interface/PhotonIdUtils.h"
#include "flashgg/DataFormats/interface/DiPhotonCandidate.h"
#include "RecoEgamma/EgammaTools/interface/EffectiveAreas.h"
// #include "RecoEgamma/EgammaTools/plugins/EGExtraInfoModifierFromDB.cc"
#include "CommonTools/CandAlgos/interface/ModifyObjectValueBase.h"
#include "flashgg/Taggers/src/IsolationCorrection.C"

#include "TFile.h"
#include "TGraph.h"

using namespace std;
using namespace edm;
using namespace reco;

namespace flashgg {

    class DiPhotonWithUpdatedPhoIdMVAProducer : public edm::EDProducer
    {
    public:
        DiPhotonWithUpdatedPhoIdMVAProducer( const edm::ParameterSet & );
        void produce( edm::Event &, const edm::EventSetup & ) override;

        void storePhotonRegressions(flashgg::DiPhotonCandidate & diph, const std::string & label);
        void updatePhotonRegressions(flashgg::DiPhotonCandidate & diph);

    private:
        float correctPhoton( flashgg::Photon & ph );
        void correctPhoton_non5x5( flashgg::Photon & ph ) ;
        void storeRegression(flashgg::Photon & cand, const std::string & label);

        edm::EDGetTokenT<edm::View<flashgg::DiPhotonCandidate> > token_;
        edm::EDGetTokenT<double> rhoToken_;
        PhotonIdUtils phoTools_;
        edm::FileInPath phoIdMVAweightfileEB_, phoIdMVAweightfileEE_, correctionFile_, non5x5correctionFile_;
        shared_ptr<ModifyObjectValueBase> regress_;
        bool reRunRegression_, reRunRegressionOnData_;
        bool correctInputs_;
        bool debug_;
        //        std::vector<TGraph*> corrections_;
        std::vector<std::unique_ptr<TGraph> > corrections_;

        bool doNon5x5transformation_;
        std::vector<std::unique_ptr<TGraph> > non5x5corrections_;

        bool useNewPhoId_;

        EffectiveAreas _effectiveAreas;
        vector<double> _phoIsoPtScalingCoeff;
        double _phoIsoCutoff;

        bool keepInitialEnergy_;
        bool _doIsoCorrection;
        unique_ptr<IsolationCorrection> _isoCorrector;
    };

    DiPhotonWithUpdatedPhoIdMVAProducer::DiPhotonWithUpdatedPhoIdMVAProducer( const edm::ParameterSet &ps ) :
        token_(consumes<edm::View<flashgg::DiPhotonCandidate> >(ps.getParameter<edm::InputTag>("src"))),
        rhoToken_( consumes<double>( ps.getParameter<edm::InputTag>( "rhoFixedGridCollection" ) ) ),
        regress_(0),
        debug_( ps.getParameter<bool>( "Debug" ) ),
        _effectiveAreas((ps.getParameter<edm::FileInPath>("effAreasConfigFile")).fullPath()),
        _phoIsoPtScalingCoeff(ps.getParameter<std::vector<double >>("phoIsoPtScalingCoeff")),
        _phoIsoCutoff(ps.getParameter<double>("phoIsoCutoff")),
        keepInitialEnergy_(ps.getParameter<bool>("keepInitialEnergy")),
        _doIsoCorrection(ps.getParameter<bool>("doIsoCorrection"))
    {
        if (_doIsoCorrection) {
            _isoCorrector = make_unique<IsolationCorrection>(ps.getParameter<edm::FileInPath>("isoCorrectionFile").fullPath().c_str());
        }

        useNewPhoId_ = ps.getParameter<bool>( "useNewPhoId" );
        phoIdMVAweightfileEB_ = ps.getParameter<edm::FileInPath>( "photonIdMVAweightfile_EB" );
        phoIdMVAweightfileEE_ = ps.getParameter<edm::FileInPath>( "photonIdMVAweightfile_EE" );
        if(useNewPhoId_){
            phoIdMVAweightfileEB_ = ps.getParameter<edm::FileInPath>( "photonIdMVAweightfile_EB_new" );
            phoIdMVAweightfileEE_ = ps.getParameter<edm::FileInPath>( "photonIdMVAweightfile_EE_new" );
        }
        phoTools_.setupMVA( phoIdMVAweightfileEB_.fullPath(), phoIdMVAweightfileEE_.fullPath(), useNewPhoId_ );

        correctInputs_ = ps.existsAs<edm::FileInPath>("correctionFile") ? true: false;
        if (correctInputs_) {
            correctionFile_ = ps.getParameter<edm::FileInPath>( "correctionFile" );
            TFile* f = TFile::Open(correctionFile_.fullPath().c_str());
            corrections_.emplace_back((TGraph*)((TGraph*) f->Get("transffull5x5R9EB"))->Clone() );
            corrections_.emplace_back((TGraph*)((TGraph*) f->Get("transfEtaWidthEB"))->Clone() );
            corrections_.emplace_back((TGraph*)((TGraph*) f->Get("transfS4EB"))->Clone() );
            corrections_.emplace_back((TGraph*)((TGraph*) f->Get("transffull5x5sieieEB"))->Clone() );
            corrections_.emplace_back((TGraph*)((TGraph*) f->Get("transffull5x5R9EE"))->Clone() );
            corrections_.emplace_back((TGraph*)((TGraph*) f->Get("transfEtaWidthEE"))->Clone() );
            corrections_.emplace_back((TGraph*)((TGraph*) f->Get("transfS4EE"))->Clone() );
            corrections_.emplace_back((TGraph*)((TGraph*) f->Get("transffull5x5sieieEE"))->Clone() );
            f->Close();
        }

        doNon5x5transformation_ =ps.getParameter<bool>( "doNon5x5transformation" );
        if (doNon5x5transformation_) {
            non5x5correctionFile_ = ps.getParameter<edm::FileInPath>( "non5x5correctionFile" );
            TFile* non5x5_f = TFile::Open(non5x5correctionFile_.fullPath().c_str());
            non5x5corrections_.emplace_back((TGraph*)((TGraph*) non5x5_f->Get("transfr9EB"))->Clone() );
            non5x5corrections_.emplace_back((TGraph*)((TGraph*) non5x5_f->Get("transfsieieEB"))->Clone() );
            non5x5corrections_.emplace_back((TGraph*)((TGraph*) non5x5_f->Get("transfsipipEB"))->Clone() );
            non5x5corrections_.emplace_back((TGraph*)((TGraph*) non5x5_f->Get("transfsieipEB"))->Clone() );
            non5x5corrections_.emplace_back((TGraph*)((TGraph*) non5x5_f->Get("transfr9EE"))->Clone() );
            non5x5corrections_.emplace_back((TGraph*)((TGraph*) non5x5_f->Get("transfsieieEE"))->Clone() );
            non5x5corrections_.emplace_back((TGraph*)((TGraph*) non5x5_f->Get("transfsipipEE"))->Clone() );
            non5x5corrections_.emplace_back((TGraph*)((TGraph*) non5x5_f->Get("transfsieipEE"))->Clone() );
            non5x5_f->Close();
        }

        edm::ConsumesCollector sumes(consumesCollector());
        reRunRegression_ = ps.getParameter<bool>("reRunRegression");
        if( reRunRegression_ ) {
            reRunRegressionOnData_ = ps.getParameter<bool>("reRunRegressionOnData");
            // regress_ = new EGExtraInfoModifierFromDB(ps.getParameter<edm::ParameterSet>("regressionConfig"));
            regress_.reset(ModifyObjectValueFactory::get()->create( "EGExtraInfoModifierFromDB", ps.getParameter<edm::ParameterSet>("regressionConfig") )); 
            regress_->setConsumes(sumes);
        }

        produces<std::vector<flashgg::DiPhotonCandidate> >();
    }

    void DiPhotonWithUpdatedPhoIdMVAProducer::updatePhotonRegressions(flashgg::DiPhotonCandidate & diph)
    {
            regress_->modifyObject(diph.getLeadingPhoton());
            regress_->modifyObject(diph.getSubLeadingPhoton());
    }

    void DiPhotonWithUpdatedPhoIdMVAProducer::storePhotonRegressions(flashgg::DiPhotonCandidate & diph, const std::string & label)
    {
            flashgg::Photon p,p2;
            storeRegression(p, label);
            storeRegression(p2, label);
    }

    void DiPhotonWithUpdatedPhoIdMVAProducer::storeRegression(flashgg::Photon & ph, const std::string & label)
    {
        ph.addUserFloat(label + "_regr_E", ph.energyCorrections().regression2Energy);
        ph.addUserFloat(label + "_regr_E_err", ph.energyCorrections().regression2EnergyError);
    }

    float DiPhotonWithUpdatedPhoIdMVAProducer::correctPhoton( flashgg::Photon & ph ) 
    {
        if (this->debug_) {
            std::cout << ph.full5x5_r9() << std::endl;
            std::cout << ph.old_r9() << std::endl;
        }
        size_t corr_index =  ph.isEB() ? 0 : 4;
        reco::Photon::ShowerShape newShowerShapes = ph.full5x5_showerShapeVariables();
        ph.addUserFloat("uncorr_r9",ph.full5x5_r9());
        ph.addUserFloat("uncorr_etaWidth",ph.superCluster()->etaWidth());
        ph.addUserFloat("uncorr_s4",ph.s4());
        ph.addUserFloat("uncorr_sigmaIetaIeta",ph.full5x5_sigmaIetaIeta());

        newShowerShapes.e3x3 = corrections_[corr_index+0]->Eval(ph.full5x5_r9())*ph.superCluster()->rawEnergy();
        newShowerShapes.sigmaIetaIeta = corrections_[corr_index+3]->Eval(ph.full5x5_sigmaIetaIeta());

        float correctedEtaWidth = corrections_[corr_index+1]->Eval(ph.superCluster()->etaWidth());
        ph.getSuperCluster()->setEtaWidth(correctedEtaWidth);
        ph.setS4(corrections_[corr_index+2]->Eval(ph.s4()));

        ph.full5x5_setShowerShapeVariables(newShowerShapes);
        return correctedEtaWidth;
    }

    void DiPhotonWithUpdatedPhoIdMVAProducer::correctPhoton_non5x5( flashgg::Photon & ph ) 
    {
        // if (this->debug_) {
        //     std::cout << ph.full5x5_r9() << std::endl;
        //     std::cout << ph.old_r9() << std::endl;
        // }
        size_t corr_index =  ph.isEB() ? 0 : 4;
        reco::Photon::ShowerShape non5x5ShowerShapes = ph.showerShapeVariables();

        ph.addUserFloat("uncorr_non5x5_r9",ph.old_r9());
        ph.addUserFloat("uncorr_non5x5_sigmaIetaIeta",non5x5ShowerShapes.sigmaIetaIeta);
        ph.addUserFloat("uncorr_non5x5_sigmaIphiIphi",non5x5ShowerShapes.sigmaIphiIphi);
        ph.addUserFloat("uncorr_non5x5_sigmaIetaIphi",non5x5ShowerShapes.sigmaIetaIphi);

        non5x5ShowerShapes.e3x3          = non5x5corrections_[corr_index+0]->Eval(ph.old_r9())*ph.superCluster()->rawEnergy();
        non5x5ShowerShapes.sigmaIetaIeta = non5x5corrections_[corr_index+1]->Eval(non5x5ShowerShapes.sigmaIetaIeta);
        non5x5ShowerShapes.sigmaIphiIphi = non5x5corrections_[corr_index+2]->Eval(non5x5ShowerShapes.sigmaIphiIphi);
        non5x5ShowerShapes.sigmaIetaIphi = non5x5corrections_[corr_index+3]->Eval(non5x5ShowerShapes.sigmaIetaIphi);

        ph.setShowerShapeVariables(non5x5ShowerShapes);
        return;
    }

    void DiPhotonWithUpdatedPhoIdMVAProducer::produce( edm::Event &evt, const edm::EventSetup & es)
    {
        edm::Handle<edm::View<flashgg::DiPhotonCandidate> > objects;
        evt.getByToken( token_, objects );

        edm::Handle<double> rhoHandle;
        evt.getByToken( rhoToken_, rhoHandle );
        const double rhoFixedGrd = *( rhoHandle.product() );

        if( reRunRegression_ ) {
            if( reRunRegressionOnData_ || ! evt.isRealData() ) {
                regress_->setEvent(evt);
                regress_->setEventContent(es);
            }
        }
        
        unique_ptr<std::vector<flashgg::DiPhotonCandidate> > out_obj( new std::vector<flashgg::DiPhotonCandidate>() );

        for (const auto & obj : *objects) {}
        evt.put( std::move(out_obj) );
    }
}

typedef flashgg::DiPhotonWithUpdatedPhoIdMVAProducer FlashggDiPhotonWithUpdatedPhoIdMVAProducer;
DEFINE_FWK_MODULE( FlashggDiPhotonWithUpdatedPhoIdMVAProducer );

// Local Variables:
// mode:c++
// indent-tabs-mode:nil
// tab-width:4
// c-basic-offset:4
// End:
// vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
