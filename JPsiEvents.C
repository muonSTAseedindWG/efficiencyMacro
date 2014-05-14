
TChain *chain = new TChain("eventsTree");

std::vector<float> *T_Gen_Muon_PtElephat;
std::vector<float> *T_Gen_Muon_Pt;
std::vector<float> *T_Gen_Muon_Eta;
std::vector<float> *T_Gen_Muon_Phi;
std::vector<float> *T_Gen_Muon_Energy;
std::vector<float> *T_Gen_Muon_Px;
std::vector<float> *T_Gen_Muon_Py;
std::vector<float> *T_Gen_Muon_Pz;
std::vector<float> *T_Gen_Muon_tpPt;
std::vector<float> *T_Muon_Eta;
std::vector<float> *T_Muon_Phi;
std::vector<float> *T_Muon_Pt;
std::vector<bool> *T_Muon_IsGlobalMuon;
std::vector<bool> *T_Muon_IsTrackerMuon;
std::vector<int> *T_Gen_Muon_PDGid;
std::vector<int> *T_Gen_Muon_status;
std::vector<int> *T_Gen_Muon_MotherID;
std::vector<int> *T_Gen_Muon_FoundSTA;
std::vector<int> *T_Gen_Muon_FoundSTAseed;
std::vector<int> *T_Gen_Muon_STAseedcrudeMaching;
std::vector<int> *T_Seed_Muon_nHits;
std::vector<float> *T_Gen_Muon_STAseedEta;
std::vector<float> *T_Gen_Muon_STAseedPhi;
std::vector<float> *T_Seed_Muon_Eta;
std::vector<float> *T_Seed_Muon_Phi;
std::vector<float> *T_Hits_Muon_STAseedEta;
std::vector<float> *T_Hits_Muon_STAseedPhi;
std::vector<int> *T_Seed_Muon_refFirstHit;
int T_Event_EventNumber;



TFile *myOutFile;

TTree *mytree;
float mu_eta;
float mu_pt;
int   foundSTAseed;
int   foundSTAseedother;
int   foundSTAseedcrude;
int   foundSTA;
int   STAseedmaxNbOfSharedHits;
int   nbOfSeedInEvent;
int   STAseedNbOfSeedWithSharedHits;
int   STAseedtotNumberOfSharedHits;
int   STAseedtotNumberHits;
bool  foundRECO;
bool  foundRECOglobal;
bool  foundRECOtracker;
float pair_mass;
float pair_dR;
float pair_dEta;
float pair_dPhi;
float pair_pt;
float seedLargestDphi;
float seedLowestDphi;
float seedLargestDeta;
float seedLowestDeta;
int eventNumber;


float deltaR(float phi1, float phi2, float eta1, float eta2)
{
    float dphi=deltaPhi(phi1,phi2);
    float deta=fabs(eta1-eta2);
    float dr = sqrt(dphi*dphi+ deta*deta);
    return dr;
}

float deltaPhi(float phi1, float phi2)
{
    float dphi;
    if(phi1<0) phi1+=2*TMath::Pi();
    if(phi2<0) phi2+=2*TMath::Pi();
    dphi=fabs(phi1-phi2);
    if(dphi>2*TMath::Pi()) dphi-=2*TMath::Pi();
    if(dphi>TMath::Pi()) dphi=2*TMath::Pi()-dphi;
    return dphi;
}

int findTheRecoSeed(int iteOfMatchedSeed){ // find the index of the reco seed corresponding to the matched seed
    float matchSeedEta = T_Gen_Muon_STAseedEta->at(iteOfMatchedSeed);
    float matchSeedPhi = T_Gen_Muon_STAseedPhi->at(iteOfMatchedSeed);
    int nbOfSeed = T_Seed_Muon_nHits->size();
    int rightSeedIndex = -1;
    for (int i = 0 ; i < nbOfSeed ; i++){
        int refFirstHit = T_Seed_Muon_refFirstHit->at(i);
        float recoSeedEta = T_Hits_Muon_STAseedEta->at(refFirstHit);
        float recoSeedPhi = T_Hits_Muon_STAseedPhi->at(refFirstHit);
        float deltaRseed = deltaR(recoSeedPhi, matchSeedPhi, recoSeedEta,matchSeedEta);
        if (deltaRseed<0.00001) rightSeedIndex = i;
    }
    if (rightSeedIndex==-1) cout << "warning : no reco seed found ! ! ! ! " << endl;
    return rightSeedIndex;
}

int findTheNumberOfCommunHits(int iteSeedOne, int iteSeedTwo){
    if ( T_Seed_Muon_nHits->size()<2) {
        cout << "warning : only 1 seed found ! " << endl;
        return 0;
    }
    int nbOfCommunHit = 0;
    int refFirstHit1 = T_Seed_Muon_refFirstHit->at(iteSeedOne);
    int nbOfHit1 = T_Seed_Muon_nHits->at(iteSeedOne);
    
    int refFirstHit2 = T_Seed_Muon_refFirstHit->at(iteSeedTwo);
    int nbOfHit2 = T_Seed_Muon_nHits->at(iteSeedTwo);
    for (int i = refFirstHit1; i < (refFirstHit1+nbOfHit1) ; i++){
        float recoSeedEta1 = T_Hits_Muon_STAseedEta->at(i);
        float recoSeedPhi1 = T_Hits_Muon_STAseedPhi->at(i);
        for (int j = refFirstHit2; j < (refFirstHit2+nbOfHit2) ; j++){
            float recoSeedEta2 = T_Hits_Muon_STAseedEta->at(j);
            float recoSeedPhi2 = T_Hits_Muon_STAseedPhi->at(j);
            float deltaRHits = deltaR(recoSeedPhi1, recoSeedPhi2, recoSeedEta1,recoSeedEta2);
            if (deltaRHits < 0.00001) {
                nbOfCommunHit++;
            }
        }
    }
    return nbOfCommunHit;
}


int findMaxNumberOfHit(int iteSeed, int *outNbOfSeed){
    int theNumberOfComHits;
    int theMaxNumberOfComHits=0;
    int nbOfSeedWithCommunHits = 0 ;
    for (int i =0 ; i < T_Seed_Muon_nHits->size(); i++){
        if (i==iteSeed) continue;
        theNumberOfComHits = findTheNumberOfCommunHits(iteSeed, i);
        if (theNumberOfComHits > theMaxNumberOfComHits) {
            theMaxNumberOfComHits = theNumberOfComHits;
            nbOfSeedWithCommunHits++;
        }
    }
    *outNbOfSeed = nbOfSeedWithCommunHits;
    return theMaxNumberOfComHits;
}

int findNbOfSharedHits(int iteSeed){
    int nbOfCommunHit = 0;
    int refFirstHit = T_Seed_Muon_refFirstHit->at(iteSeed);
    int nbOfHit = T_Seed_Muon_nHits->at(iteSeed);
    int totNbOfHits = T_Hits_Muon_STAseedEta->size();
    int totNbOfHits2 = T_Hits_Muon_STAseedPhi->size();
    for (int i = refFirstHit ; i < (refFirstHit+nbOfHit) ; i++){
        float recoHitEta1 = T_Hits_Muon_STAseedEta->at(i);
        float recoHitPhi1 = T_Hits_Muon_STAseedPhi->at(i);
        bool foundAHit = false;
        for (int j = 0 ; j < totNbOfHits ; j++){
            if ((j>=refFirstHit)&&(j<(refFirstHit+nbOfHit))) continue;
            float recoHitEta2 = T_Hits_Muon_STAseedEta->at(j);
            float recoHitPhi2 = T_Hits_Muon_STAseedPhi->at(j);
            float deltaRHits = deltaR(recoHitPhi1, recoHitPhi2, recoHitEta1, recoHitEta2);
            if (deltaRHits<0.0001) foundAHit=true;
        }
        if (foundAHit) nbOfCommunHit++;
    }
    return nbOfCommunHit;
}

float computeTheLargestdPhi(int iteSeed){
    int refFirstHit = T_Seed_Muon_refFirstHit->at(iteSeed);
    int nbOfHit = T_Seed_Muon_nHits->at(iteSeed);
    float largestDphi = 0;
    if (nbOfHit<2) return -1;
    float phiRef = T_Hits_Muon_STAseedPhi->at(refFirstHit);
    for (int i = refFirstHit+1 ; i < (refFirstHit+nbOfHit) ; i++){
        float thePhi = T_Hits_Muon_STAseedPhi->at(i);
        float theDeltaPhi = fabs(phiRef-thePhi);
        if (theDeltaPhi>largestDphi) largestDphi = theDeltaPhi;
    }
    return largestDphi;
}

float computeTheLowestdPhi(int iteSeed){///function cut on paste so largest means lowest inside
    int refFirstHit = T_Seed_Muon_refFirstHit->at(iteSeed);
    int nbOfHit = T_Seed_Muon_nHits->at(iteSeed);
    float largestDphi = 100;
    if (nbOfHit<2) return -1;
    float phiRef = T_Hits_Muon_STAseedPhi->at(refFirstHit);
    for (int i = refFirstHit+1 ; i < (refFirstHit+nbOfHit) ; i++){
        float thePhi = T_Hits_Muon_STAseedPhi->at(i);
        float theDeltaPhi = fabs(phiRef-thePhi);
        if (theDeltaPhi<largestDphi) largestDphi = theDeltaPhi;
    }
    return largestDphi;
}

float computeTheLargestdEta(int iteSeed){
    int refFirstHit = T_Seed_Muon_refFirstHit->at(iteSeed);
    int nbOfHit = T_Seed_Muon_nHits->at(iteSeed);
    float largestDeta = 0;
    if (nbOfHit<2) return -1;
    float etaRef = T_Hits_Muon_STAseedEta->at(refFirstHit);
    for (int i = refFirstHit+1 ; i < (refFirstHit+nbOfHit) ; i++){
        float theEta = T_Hits_Muon_STAseedEta->at(i);
        float theDeltaEta = fabs(etaRef-theEta);
        if (theDeltaEta>largestDeta) largestDeta = theDeltaEta;
    }
    return largestDeta;
}

float computeTheLowestdEta(int iteSeed){///function cut on paste so largest means lowest inside
    int refFirstHit = T_Seed_Muon_refFirstHit->at(iteSeed);
    int nbOfHit = T_Seed_Muon_nHits->at(iteSeed);
    float largestDeta = 100;
    if (nbOfHit<2) return -1;
    float etaRef = T_Hits_Muon_STAseedEta->at(refFirstHit);
    for (int i = refFirstHit+1 ; i < (refFirstHit+nbOfHit) ; i++){
        float theEta = T_Hits_Muon_STAseedEta->at(i);
        float theDeltaEta = fabs(etaRef-theEta);
        if (theDeltaEta<largestDeta) largestDeta = theDeltaEta;
    }
    return largestDeta;
}

JPsiEvents(TString minitreeFile){
    
    myOutFile = new TFile("minitree_"+nameSample+".root","RECREATE");
    
    mytree = new TTree("eventsTree","");
    mytree->Branch("mu_eta", &mu_eta, "mu_eta/F");
    mytree->Branch("mu_pt", &mu_pt, "mu_pt/F");
    mytree->Branch("pair_mass", &pair_mass, "pair_mass/F");
    mytree->Branch("pair_pt", &pair_pt, "pair_pt/F");
    mytree->Branch("pair_dR", &pair_dR, "pair_dR/F");
    mytree->Branch("pair_dEta", &pair_dEta, "pair_dEta/F");
    mytree->Branch("pair_dPhi", &pair_dPhi, "pair_dPhi/F");
    mytree->Branch("foundSTAseed", &foundSTAseed, "foundSTAseed/I");
    mytree->Branch("STAseedmaxNbOfSharedHits", &STAseedmaxNbOfSharedHits, "STAseedmaxNbOfSharedHits/I");
    mytree->Branch("nbOfSeedInEvent", &nbOfSeedInEvent, "nbOfSeedInEvent/I");
    mytree->Branch("STAseedNbOfSeedWithSharedHits", &STAseedNbOfSeedWithSharedHits, "STAseedNbOfSeedWithSharedHits/I");
    mytree->Branch("STAseedtotNumberOfSharedHits", &STAseedtotNumberOfSharedHits, "STAseedtotNumberOfSharedHits/I");
    mytree->Branch("STAseedtotNumberHits", &STAseedtotNumberHits, "STAseedtotNumberHits/I");
    mytree->Branch("foundSTAseedcrude", &foundSTAseedcrude, "foundSTAseedcrude/I");
    mytree->Branch("foundSTA", &foundSTA, "foundSTA/I");
    mytree->Branch("foundRECO", &foundRECO, "foundRECO/B");
    mytree->Branch("foundRECOglobal", &foundRECOglobal, "foundRECOglobal/B");
    mytree->Branch("foundRECOtracker", &foundRECOtracker, "foundRECOtracker/B");
    mytree->Branch("eventNumber", &eventNumber, "eventNumber/I");
    mytree->Branch("seedLargestDphi", &seedLargestDphi, "seedLargestDphi/F");
    mytree->Branch("seedLowestDphi", &seedLowestDphi, "seedLowestDphi/F");
    mytree->Branch("seedLargestDeta", &seedLargestDeta, "seedLargestDeta/F");
    mytree->Branch("seedLowestDeta", &seedLowestDeta, "seedLowestDeta/F");
    
    
    chain->Add(minitreeFile);
    chain->SetBranchAddress("T_Gen_Muon_Eta",&T_Gen_Muon_Eta);
    chain->SetBranchAddress("T_Gen_Muon_Phi",&T_Gen_Muon_Phi);
    chain->SetBranchAddress("T_Gen_Muon_Pt",&T_Gen_Muon_Pt);
    chain->SetBranchAddress("T_Gen_Muon_Energy",&T_Gen_Muon_Energy);
    chain->SetBranchAddress("T_Gen_Muon_Px",&T_Gen_Muon_Px);
    chain->SetBranchAddress("T_Gen_Muon_Py",&T_Gen_Muon_Py);
    chain->SetBranchAddress("T_Gen_Muon_Pz",&T_Gen_Muon_Pz);
    chain->SetBranchAddress("T_Gen_Muon_tpPt",&T_Gen_Muon_tpPt);
    chain->SetBranchAddress("T_Gen_Muon_PDGid",&T_Gen_Muon_PDGid);
    chain->SetBranchAddress("T_Gen_Muon_status",&T_Gen_Muon_status);
    chain->SetBranchAddress("T_Gen_Muon_MotherID",&T_Gen_Muon_MotherID);
    chain->SetBranchAddress("T_Gen_Muon_FoundSTA",&T_Gen_Muon_FoundSTA);
    chain->SetBranchAddress("T_Gen_Muon_FoundSTAseed",&T_Gen_Muon_FoundSTAseed);
    chain->SetBranchAddress("T_Gen_Muon_STAseedcrudeMaching",&T_Gen_Muon_STAseedcrudeMaching);
    chain->SetBranchAddress("T_Gen_Muon_STAseedEta",&T_Gen_Muon_STAseedEta);
    chain->SetBranchAddress("T_Gen_Muon_STAseedPhi",&T_Gen_Muon_STAseedPhi);
    chain->SetBranchAddress("T_Seed_Muon_nHits",&T_Seed_Muon_nHits);
    chain->SetBranchAddress("T_Seed_Muon_Eta",&T_Seed_Muon_Eta);
    chain->SetBranchAddress("T_Seed_Muon_Phi",&T_Seed_Muon_Phi);
    chain->SetBranchAddress("T_Seed_Muon_refFirstHit",&T_Seed_Muon_refFirstHit);
    chain->SetBranchAddress("T_Hits_Muon_STAseedEta",&T_Hits_Muon_STAseedEta);
    chain->SetBranchAddress("T_Hits_Muon_STAseedPhi",&T_Hits_Muon_STAseedPhi);
    
    chain->SetBranchAddress("T_Muon_Eta",&T_Muon_Eta);
    chain->SetBranchAddress("T_Muon_Phi",&T_Muon_Phi);
    chain->SetBranchAddress("T_Muon_Pt",&T_Muon_Pt);
    chain->SetBranchAddress("T_Muon_IsGlobalMuon",&T_Muon_IsGlobalMuon);
    chain->SetBranchAddress("T_Muon_IsTrackerMuon",&T_Muon_IsTrackerMuon);
    chain->SetBranchAddress("T_Event_EventNumber",&T_Event_EventNumber);
    
    
    
    
    
    int NbEntries = chain->GetEntries();
    cout << "nb tot events=" << NbEntries << endl;
    for (int i=0 ;i<NbEntries ; i++){
        chain->GetEntry(i);
        TLorentzVector sumMuons;

        if ((i%1000)==0) cout << "i=" << i << endl;

        
        //start to loop on double muon events
        for (int j = 0  ; j < T_Gen_Muon_Eta->size() ; j++){
            if (!((fabs(T_Gen_Muon_PDGid->at(j))==13)&&(T_Gen_Muon_status->at(j)==1))) continue; //keep only muon comming from JPSI
            TLorentzVector* muon1 = new TLorentzVector(T_Gen_Muon_Px->at(j), T_Gen_Muon_Py->at(j), T_Gen_Muon_Pz->at(j), T_Gen_Muon_Energy->at(j));
            for (int k= 0 ; k < T_Gen_Muon_Eta->size() ; k++){
                if (j==k) continue;
                if (!((fabs(T_Gen_Muon_PDGid->at(k))==13)&&(T_Gen_Muon_status->at(k)==1))) continue;
                TLorentzVector* muon2 = new TLorentzVector(T_Gen_Muon_Px->at(k), T_Gen_Muon_Py->at(k), T_Gen_Muon_Pz->at(k), T_Gen_Muon_Energy->at(k));
                sumMuons = *muon1 + *muon2;
                mu_eta = T_Gen_Muon_Eta->at(k);
                mu_pt  = T_Gen_Muon_Pt->at(k);
                
                pair_mass = sumMuons.M();
                pair_pt = sumMuons.Pt();
                pair_dR = deltaR(T_Gen_Muon_Phi->at(j), T_Gen_Muon_Phi->at(k), T_Gen_Muon_Eta->at(j),  T_Gen_Muon_Eta->at(k));
                pair_dEta = deltaR(0, T_Gen_Muon_Phi->at(k), 0,  T_Gen_Muon_Eta->at(k));
                pair_dPhi = deltaR(T_Gen_Muon_Phi->at(j), T_Gen_Muon_Phi->at(k), 0,  0);
                if ((T_Gen_Muon_tpPt->at(j)==-1)||(T_Gen_Muon_tpPt->at(k)==-1)) {
                    cout << "did not found any tp corresponding the gen muon, will be unsefull to check one day why it happens... " << endl;
                    continue;
                }
                nbOfSeedInEvent = T_Seed_Muon_nHits->size();
                foundSTAseed = T_Gen_Muon_FoundSTAseed->at(k);
                int recoSeedFound;
                int theNumberOfHitsWithCommunHits;
                int maxNumberOfCommunHits;
                int theNumberOfSharedHits;
                if (foundSTAseed) {
                    recoSeedFound = findTheRecoSeed(k);
                    maxNumberOfCommunHits = findMaxNumberOfHit(recoSeedFound, &theNumberOfHitsWithCommunHits);
                    theNumberOfSharedHits = findNbOfSharedHits(recoSeedFound);
                    STAseedmaxNbOfSharedHits = maxNumberOfCommunHits;
                    STAseedNbOfSeedWithSharedHits = theNumberOfHitsWithCommunHits;
                    STAseedtotNumberOfSharedHits = theNumberOfSharedHits;
                    STAseedtotNumberHits = T_Seed_Muon_nHits->at(recoSeedFound);
                    seedLargestDphi = computeTheLargestdPhi(recoSeedFound);
                    seedLowestDphi  = computeTheLowestdPhi(recoSeedFound);
                    seedLargestDeta = computeTheLargestdEta(recoSeedFound);
                    seedLowestDeta  = computeTheLowestdEta(recoSeedFound);

                }
                else{
                    STAseedmaxNbOfSharedHits = -1;
                    STAseedNbOfSeedWithSharedHits = -1;
                    STAseedtotNumberOfSharedHits = -1;
                }
                foundSTA  = T_Gen_Muon_FoundSTA->at(k);
                foundSTAseedcrude = T_Gen_Muon_STAseedcrudeMaching->at(k);
                int nbRecoMuons = T_Muon_Pt->size();
                int iteMinDr = -1;
                float minDr = 1000;
                for (int iteReco = 0 ; iteReco < nbRecoMuons ; iteReco++){
                    float theDeltaR = deltaR(  T_Gen_Muon_Phi->at(k), T_Muon_Phi->at(iteReco), T_Gen_Muon_Eta->at(k), T_Muon_Eta->at(iteReco));
                    if (theDeltaR<minDr){
                        minDr = theDeltaR;
                        iteMinDr = iteReco;
                    }
                }
                if (minDr < 0.1){//match the reco muon with the gen info !
                    foundRECO = 1;
                    foundRECOglobal = T_Muon_IsGlobalMuon->at(iteMinDr);
                    foundRECOtracker = T_Muon_IsTrackerMuon->at(iteMinDr);
                }
                else {
                    foundRECO = 0;
                    foundRECOglobal = 0;
                    foundRECOtracker = 0;
                }
                eventNumber = T_Event_EventNumber;
                mytree->Fill();
            }
        }


    }
    myOutFile->Write();
    myOutFile->Close();
    delete myOutFile;
}
