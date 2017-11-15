#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/EDMException.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "flashgg/DataFormats/interface/DiPhotonCandidate.h"
#include "flashgg/MicroAOD/interface/PhotonIdUtils.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "flashgg/MicroAOD/interface/VertexSelectorBase.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "flashgg/DataFormats/interface/VertexCandidateMap.h"
#include "DataFormats/EgammaCandidates/interface/Conversion.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"

#include <map>

using namespace edm;
using namespace std;

bool compareSumPts( const std::pair<double,int>& a, const std::pair<double,int>& b){
    return a.first > b.first;
}

namespace flashgg {

    class DiPhotonProducer : public EDProducer
    {

    public:
        DiPhotonProducer( const ParameterSet & );
    private:
        void produce( Event &, const EventSetup & ) override;
        EDGetTokenT<View<reco::Vertex> > vertexToken_;
        EDGetTokenT<View<pat::Photon> > photonToken_;
        EDGetTokenT< VertexCandidateMap > vertexCandidateMapToken_;
        unique_ptr<VertexSelectorBase> vertexSelector_;
        EDGetTokenT<View<reco::Conversion> > conversionToken_;
        EDGetTokenT<reco::BeamSpot>  beamSpotToken_;
        EDGetTokenT<View<reco::Conversion> > conversionTokenSingleLeg_;
        bool useSingleLeg_;
        unsigned int maxJetCollections_;
        EDGetTokenT<View<reco::GenParticle> >      genPartToken_;
    };

    DiPhotonProducer::DiPhotonProducer( const ParameterSet &iConfig ) :
        vertexToken_( consumes<View<reco::Vertex> >( iConfig.getParameter<InputTag> ( "VertexTag" ) ) ),
        photonToken_( consumes<View<pat::Photon> >( iConfig.getParameter<InputTag> ( "PhotonTag" ) ) ),
        vertexCandidateMapToken_( consumes<VertexCandidateMap>( iConfig.getParameter<InputTag>( "VertexCandidateMapTag" ) ) ),
        conversionToken_( consumes<View<reco::Conversion> >( iConfig.getParameter<InputTag>( "ConversionTag" ) ) ),
        beamSpotToken_( consumes<reco::BeamSpot >( iConfig.getParameter<InputTag>( "beamSpotTag" ) ) ),
        conversionTokenSingleLeg_( consumes<View<reco::Conversion> >( iConfig.getParameter<InputTag>( "ConversionTagSingleLeg" ) ) ),
        maxJetCollections_( iConfig.getParameter<unsigned int>( "MaxJetCollections" ) ),
        genPartToken_( consumes<View<reco::GenParticle> >( iConfig.getParameter<InputTag> ( "GenParticleTag" ) ) )
    {
        const std::string &VertexSelectorName = iConfig.getParameter<std::string>( "VertexSelectorName" );
        // VertexSelectorBase *vsb = new VertexSelectorBase(iConfig); // already contains the needed VertexSelectorName parameter
        // vertexSelector_.reset( vsb );
        vertexSelector_.reset( FlashggVertexSelectorFactory::get()->create( VertexSelectorName, iConfig ) );
        useSingleLeg_ = iConfig.getParameter<bool>( "useSingleLeg" );
        produces<vector<flashgg::DiPhotonCandidate> >();
        produces< vector<pat::Photon> >(); // SK Note: changed this code to return a vector with diphotonColl[0] since that's what's used, and it more easily plugs into the Alabama/Rutgers code.
    }

    void DiPhotonProducer::produce( Event &evt, const EventSetup & )
    {
        Handle<View<reco::Vertex> > primaryVertices;
        evt.getByToken( vertexToken_, primaryVertices );
        // const PtrVector<reco::Vertex>& pvPointers = primaryVertices->ptrVector();

        Handle<View<pat::Photon> > photons;
        evt.getByToken( photonToken_, photons );
        //  const PtrVector<pat::Photon>& photonPointers = photons->ptrVector();

        Handle<VertexCandidateMap> vertexCandidateMap;
        evt.getByToken( vertexCandidateMapToken_, vertexCandidateMap );

        Handle<View<reco::Conversion> > conversions;
        evt.getByToken( conversionToken_, conversions );
        // const PtrVector<reco::Conversion>& conversionPointers = conversions->ptrVector();

        Handle<reco::BeamSpot> recoBeamSpotHandle;
        evt.getByToken( beamSpotToken_, recoBeamSpotHandle );
        math::XYZPoint vertexPoint;
        //    float beamsig;
        if( recoBeamSpotHandle.isValid() ) {
            vertexPoint = recoBeamSpotHandle->position();
            //      beamsig = recoBeamSpotHandle->sigmaZ();
        }

        math::XYZPoint higgsVtx;
        if( ! evt.isRealData() ) {
            Handle<View<reco::GenParticle> > genParticles;
            evt.getByToken( genPartToken_, genParticles );
            for( unsigned int genLoop = 0 ; genLoop < genParticles->size(); genLoop++ ) {
                int pdgid = genParticles->ptrAt( genLoop )->pdgId();
                if( pdgid == 25 || pdgid == 22 ) {
                    higgsVtx = genParticles->ptrAt( genLoop )->vertex();
                    break;
                }
            }
        }

        Handle<View<reco::Conversion> > conversionsSingleLeg;
        evt.getByToken( conversionTokenSingleLeg_, conversionsSingleLeg );
        //const PtrVector<reco::Conversion>& conversionPointersSingleLeg = conversionsSingleLeg->ptrVector();

        unique_ptr<vector<DiPhotonCandidate> > diPhotonColl( new vector<DiPhotonCandidate> );
        unique_ptr<vector<pat::Photon> > selectedPhotons( new vector<pat::Photon> );
//    cout << "evt.id().event()= " << evt.id().event() << "\tevt.isRealData()= " << evt.isRealData() << "\tphotons->size()= " << photons->size() << "\tprimaryVertices->size()= " << primaryVertices->size() << endl;

        vector< std::pair<double,int> > sumPtVec;
        int pairNum = 0;
        for( unsigned int i = 0 ; i < photons->size() ; i++ ) {

            Ptr<pat::Photon> pp1 = photons->ptrAt( i );
            for( unsigned int j = i + 1 ; j < photons->size() ; j++ ) {
                Ptr<pat::Photon> pp2 = photons->ptrAt( j );

                Ptr<reco::Vertex> pvx = vertexSelector_->select( pp1, pp2, primaryVertices->ptrs(), *vertexCandidateMap, conversions->ptrs(), conversionsSingleLeg->ptrs(),
                                        vertexPoint, useSingleLeg_ );
                // Finding and storing the vertex index to check if it corresponds to the primary vertex.
                // This could be moved within the vertexSelector, but would need rewriting some interface
                int ivtx = 0;
                for( unsigned int k = 0; k < primaryVertices->size() ; k++ )
                    if( pvx == primaryVertices->ptrAt( k ) ) {
                        ivtx = k;
                        break;
                    }

                DiPhotonCandidate dipho( pp1, pp2, pvx );
                dipho.setVertexIndex( ivtx );
                dipho.setGenPV( higgsVtx );

                // Obviously the last selection has to be for this diphoton or this is wrong
                vertexSelector_->writeInfoFromLastSelectionTo( dipho );

                // store the diphoton into the collection
                dipho.makePhotonsPersistent();
                sumPtVec.push_back( std::make_pair(dipho.leadingPhoton()->pt() + dipho.subLeadingPhoton()->pt(),pairNum) );
                diPhotonColl->push_back( dipho );
                // selectedPhotons->push_back(*pp1);
                // selectedPhotons->push_back(*pp2);
                selectedPhotons->push_back(*dipho.leadingPhoton());
                selectedPhotons->push_back(*dipho.subLeadingPhoton());
                pairNum++;
            }
        }
        // Sort the final collection (descending) and put it in the event
        std::sort( diPhotonColl->begin(), diPhotonColl->end(), greater<DiPhotonCandidate>() );
        // now need to sort the sumPtVec in descending order of sumPt.  The "second" will say which photons from selectedPhotons to put now in the collection to be added to the event
        std::sort( sumPtVec.begin(), sumPtVec.end(), compareSumPts);
        unique_ptr<vector<pat::Photon> > selectedPhotons2( new vector<pat::Photon> );
        // std::cout << "====flashgg DiPhotonProducer====" << std::endl;
        for (unsigned int i=0; i < sumPtVec.size(); i++){
            int idx = sumPtVec.at(i).second;
            int idxInColl = idx * 2;
            selectedPhotons2->push_back( selectedPhotons->at(idxInColl) );
            selectedPhotons2->push_back( selectedPhotons->at(idxInColl+1) );
            // std::cout << "pair# = " << idx << std::endl;
            // std::cout << "pt1 = " << selectedPhotons->at(idxInColl).pt() << std::endl;
            // std::cout << "pt2 = " << selectedPhotons->at(idxInColl+1).pt() << std::endl;
            // std::cout << "sumPt from sumPtVec Photons = " << selectedPhotons->at(idxInColl).pt() + selectedPhotons->at(idxInColl+1).pt() << std::endl;
            // std::cout << "sumPt from sumPtVec = " << sumPtVec.at(i).first << std::endl;
            // std::cout << "mass from sumPtVec = " << (selectedPhotons->at(idxInColl).p4() + selectedPhotons->at(idxInColl+1).p4()).M() << std::endl;
            // std::cout << "sumPt from diPhotonColl = " << diPhotonColl->at(i).sumPt() << std::endl;
            // std::cout << "mass from diPhotonColl = " << (diPhotonColl->at(i).leadingPhoton()->p4() + diPhotonColl->at(i).subLeadingPhoton()->p4()).M() << std::endl;
            // std::cout << " " << std::endl;

        }


        // std::cout << "====flashgg DiPhotonProducer====" << std::endl;
        // for(unsigned int i=0; i < diPhotonColl->size(); i++){
        //     DiPhotonCandidate selectedDiPhotonCand = diPhotonColl->at(i);
        //     selectedDiPhotonCand.makePhotonsPersistent();
        //     pat::Photon leadingPho = selectedDiPhotonCand.getLeadingPhoton();
        //     pat::Photon subLeadingPho = selectedDiPhotonCand.getSubLeadingPhoton();
        //     std::cout << "diPhotonColl[" << i << "] (scEta1,scEta2,sumPt): " << leadingPho.superCluster()->eta() << ", " << subLeadingPho.superCluster()->eta() << ", " << selectedDiPhotonCand.sumPt() << std::endl;
        // }


        // SK Note: amending this to add vector with the two photons in diPhotonColl[0]
        // unique_ptr<vector<pat::Photon> > selectedPhotons( new vector<pat::Photon> );
        // // std::vector<pat::Photon> selectedPhotons;
        // if (diPhotonColl->size() > 0){
        //     DiPhotonCandidate selectedDiPhotonCand = diPhotonColl->at(0);
        //     selectedDiPhotonCand.makePhotonsPersistent();
        //     pat::Photon leadingPho = selectedDiPhotonCand.getLeadingPhoton();
        //     pat::Photon subLeadingPho = selectedDiPhotonCand.getSubLeadingPhoton();
        //     selectedPhotons->push_back(leadingPho);
        //     selectedPhotons->push_back(subLeadingPho);
        // }
        evt.put( std::move( diPhotonColl ) );
        evt.put( std::move( selectedPhotons2 ) );

    }
}

typedef flashgg::DiPhotonProducer FlashggDiPhotonProducer;
DEFINE_FWK_MODULE( FlashggDiPhotonProducer );
// Local Variables:
// mode:c++
// indent-tabs-mode:nil
// tab-width:4
// c-basic-offset:4
// End:
// vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4

