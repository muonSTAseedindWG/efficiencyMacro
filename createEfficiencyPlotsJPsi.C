
TChain *chain = new TChain("eventsTree");
TFile *myFile;

TH1F *giveEfficiency(TString nomPlot, TString variable, TString theCut, TString denomCut, float *theBins, int nbBins){
    cout << "will do " << nomPlot << endl;
    TH1F *denom = new TH1F("denom", "", nbBins, theBins);
    TH1F *num = new TH1F("num", "", nbBins, theBins);
    TString baseCut = "mu_pt>2&&abs(mu_eta)<2.4";
    chain->Draw(variable+">>denom",baseCut+"&&"+denomCut);
    chain->Draw(variable+">>num",baseCut+"&&"+denomCut+"&&"+theCut+"==1");
    
    TH1F *efficiency = (TH1F*) denom->Clone(nomPlot);
    efficiency->Sumw2();
    efficiency->Divide(num, denom, 1,1);

    delete denom;
    delete num;
    return efficiency;
}


TH1F *drawTheHisto(TString nomPlot, TString variable, TString denomCut, float *theBins, int nbBins){
    TH1F *num = new TH1F(nomPlot, "", nbBins, theBins);
    TString baseCut = "mu_pt>2&&abs(mu_eta)<2.4";
    chain->Draw(variable+">>"+nomPlot,baseCut+"&&"+denomCut);
    return num;
}

TH1F *drawTheHisto(TString nomPlot, TString variable, TString denomCut, float theMin, float theMax, int nbBins){
    TH1F *num = new TH1F(nomPlot, "", nbBins, theMin, theMax);
    TString baseCut = "mu_pt>2&&abs(mu_eta)<2.4";
    chain->Draw(variable+">>"+nomPlot,baseCut+"&&"+denomCut);
    return num;
}

void drawSetHisto(TString nomPlot, TString variable, TString denomCut, float theMin, float theMax, int nbBins){
    
    myFile->cd();
    
    TH1F *all = drawTheHisto( nomPlot, variable, denomCut, theMin, theMax, nbBins);
    all->Write();
    delete all;
 
    TH1F *passing = drawTheHisto( nomPlot+"_pass", variable, denomCut+"&&foundSTA==1", theMin, theMax, nbBins);
    passing->Write();
    delete passing;
    
    TH1F *failling = drawTheHisto( nomPlot+"_fail", variable, denomCut+"&&foundSTA==0", theMin, theMax, nbBins);
    failling->Write();
    delete failling;

}


createEfficiencyPlotsJPsi(TString dataset){
    chain->Add("minitree_"+dataset+".root");
    
    myFile = new TFile("histo_"+dataset+".root","RECREATE");
    TString nameCut[6] = {"STAseed",     "STAseedcrude",     "STA",     "RECO",      "Glb",           "Trk"};
    TString theCut[6]  = {"foundSTAseed","foundSTAseedcrude","foundSTA","foundRECO","foundRECOglobal","foundRECOtracker"};
    
    float etaBins[11] = {0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.8, 2.1, 2.4};
    int nbEtaBins = 10;

    float deltaRbins[6] = {0.05,0.1,0.2,0.3,0.5,1};
    int nbRbins = 5;
    
    //float PtBinsbins[9] = {20, 30, 40, 50, 60, 70, 80, 90, 100};
    float PtBinsbins[9] = {0, 10, 20, 30, 40, 50, 60, 70, 100};
    int nbPtBins = 8;
    
    //float PtMuBinsbins[7] = {0, 10, 20, 30, 50, 80, 100};
    float PtMuBinsbins[7] = {0, 3,  5, 7, 10, 15, 50};
    int nbPtMuBins = 6;
    
    float sharedHitsBins[7] = {0,1,2,3,4,5,6};
    int nbOfsharedHitsBins = 6;
    
    float seedInEventBins[9] = {0,1,2,3,4,5,7,10,15};
    int nbOfseedInEventBins = 8;
    
    float phiBins[7] = {0.05, 0.1, 0.15, 0.20, 0.25, 0.5, 1};
    int nbOfphiBins = 6;

    
    for(int i = 0 ; i < 6 ; i++){
        TH1F* histoEffEta = giveEfficiency("effVsEta_"+nameCut[i], "abs(mu_eta)", theCut[i],"1", etaBins, nbEtaBins);
        histoEffEta->Write();
        TH1F* histoEffEtaHigh = giveEfficiency("effVsEta_"+nameCut[i]+"_high", "abs(mu_eta)", theCut[i],"mu_pt>5", etaBins, nbEtaBins);
        histoEffEtaHigh->Write();
        TH1F* histoEffdR = giveEfficiency("effVsdR_"+nameCut[i], "pair_dR", theCut[i],"1", deltaRbins, nbRbins);
        histoEffdR->Write();
        TH1F* histoEffPt = giveEfficiency("effVsPt_"+nameCut[i], "pair_pt", theCut[i], "1", PtBinsbins, nbPtBins);
        histoEffPt->Write();
        TH1F* histoEffPt = giveEfficiency("effVsMuPt_"+nameCut[i], "mu_pt", theCut[i], "1", PtMuBinsbins, nbPtMuBins);
        histoEffPt->Write();
        TH1F* histoEffDPhi = giveEfficiency("effVsMudPhi_"+nameCut[i], "pair_dPhi", theCut[i], "1", phiBins, nbOfphiBins);
        histoEffDPhi->Write();
    }
    
    TH1F *histoEffNbHitShared = giveEfficiency("effVsNbSharedHits","STAseedtotNumberOfSharedHits", "foundSTA","foundSTAseed",sharedHitsBins, nbOfsharedHitsBins);
    histoEffNbHitShared->Write();
    
    TH1F *histoEffNbMaxHitShared = giveEfficiency("effVsNbMaxSharedHits","STAseedmaxNbOfSharedHits", "foundSTA","foundSTAseed",sharedHitsBins, nbOfsharedHitsBins);
    histoEffNbMaxHitShared->Write();
    
    TH1F *histoEffNbHit = giveEfficiency("histoEffNbHit","STAseedtotNumberHits", "foundSTA","foundSTAseed",sharedHitsBins, nbOfsharedHitsBins);
    histoEffNbHit->Write();
    
    TH1F *histoEffNbSeedShared = giveEfficiency("histoEffNbSeedShared","STAseedNbOfSeedWithSharedHits", "foundSTA","foundSTAseed",sharedHitsBins, nbOfsharedHitsBins);
    histoEffNbSeedShared->Write();
    
    
    TH1F *histoEffNbofRemainHit = giveEfficiency("histoEffNbHitRemain","(STAseedtotNumberHits-STAseedtotNumberOfSharedHits)", "foundSTA","foundSTAseed",sharedHitsBins, nbOfsharedHitsBins);
    histoEffNbofRemainHit->Write();
    
    TH1F *histoEffNbofRemainMaxHit = giveEfficiency("histoEffNbHitRemainMax","(STAseedtotNumberHits-STAseedmaxNbOfSharedHits)", "foundSTA","foundSTAseed",sharedHitsBins, nbOfsharedHitsBins);
    histoEffNbofRemainMaxHit->Write();
    
    TH1F *histoEffNbofSeed = giveEfficiency("histoEffNbofSeed","nbOfSeedInEvent", "foundSTA","foundSTAseed",seedInEventBins, nbOfseedInEventBins);
    histoEffNbofSeed->Write();
    
    
    drawSetHisto("deltaR","pair_dR","foundSTAseed",0,2, 50);
    
    drawSetHisto("nbOfSharedHitsPass","STAseedtotNumberOfSharedHits","foundSTAseed",0,7,7);
   
    drawSetHisto("nbMaxOfSharedHitsPass","STAseedmaxNbOfSharedHits","foundSTAseed",0,7,7);
    
    drawSetHisto("nbHitRemain","(STAseedtotNumberHits-STAseedtotNumberOfSharedHits)","foundSTAseed",0,7,7);
    
    drawSetHisto("nbSeedInTheEvent","nbOfSeedInEvent","foundSTAseed",0,20,20);
    
    drawSetHisto("nbOfSeedWithSharedHits","STAseedNbOfSeedWithSharedHits","foundSTAseed",0,7,7);
    
    drawSetHisto("nbOfHitsPerSeed","STAseedtotNumberHits","foundSTAseed",0,7,7);
    drawSetHisto("nbOfHitsPerSeedLowPt","STAseedtotNumberHits","foundSTAseed&&(mu_pt<5)",0,7,7);
    drawSetHisto("nbOfHitsPerSeedMediumPt","STAseedtotNumberHits","foundSTAseed&&(mu_pt<10)&&(mu_pt>5)",0,7,7);
    drawSetHisto("nbOfHitsPerSeedHighPt","STAseedtotNumberHits","foundSTAseed&&(mu_pt>10)",0,7,7);
 
    drawSetHisto("SeedLargestDphi","seedLargestDphi","foundSTAseed",0,0.3,60);
    drawSetHisto("SeedLargestDphiLowPt","seedLargestDphi","foundSTAseed&&(mu_pt<5)",0,0.3,60);
    drawSetHisto("SeedLargestDphiMediumPt","seedLargestDphi","foundSTAseed&&(mu_pt<10)&&(mu_pt>5)",0,0.3,60);
    drawSetHisto("SeedLargestDphiHighPt","seedLargestDphi","foundSTAseed&&(mu_pt>10)",0,0.3,60);
    
    drawSetHisto("SeedLargestDphiDT","seedLargestDphi","foundSTAseed&&abs(mu_eta)<0.8",0,0.3,60);
    drawSetHisto("SeedLargestDphiLowPtDT","seedLargestDphi","foundSTAseed&&(mu_pt<5)&&abs(mu_eta)<0.8",0,0.3,60);
    drawSetHisto("SeedLargestDphiMediumPtDT","seedLargestDphi","foundSTAseed&&(mu_pt<10)&&(mu_pt>5)&&abs(mu_eta)<0.8",0,0.3,60);
    drawSetHisto("SeedLargestDphiHighPtDT","seedLargestDphi","foundSTAseed&&(mu_pt>10)&&abs(mu_eta)<0.8",0,0.3,60);
  
    drawSetHisto("SeedLargestDphiCSC","seedLargestDphi","foundSTAseed&&abs(mu_eta)>0.8",0,0.3,60);
    drawSetHisto("SeedLargestDphiLowPtCSC","seedLargestDphi","foundSTAseed&&(mu_pt<5)&&abs(mu_eta)>0.8",0,0.3,60);
    drawSetHisto("SeedLargestDphiMediumPtCSC","seedLargestDphi","foundSTAseed&&(mu_pt<10)&&(mu_pt>5)&&abs(mu_eta)>0.8",0,0.3,60);
    drawSetHisto("SeedLargestDphiHighPtCSC","seedLargestDphi","foundSTAseed&&(mu_pt>10)&&abs(mu_eta)>0.8",0,0.3,60);
    
    drawSetHisto("SeedLowestDphi","seedLowestDphi","foundSTAseed",0,0.3,60);
    drawSetHisto("SeedLowestDphiLowPt","seedLowestDphi","foundSTAseed&&(mu_pt<5)",0,0.3,60);
    drawSetHisto("SeedLowestDphiMediumPt","seedLowestDphi","foundSTAseed&&(mu_pt<10)&&(mu_pt>5)",0,0.3,60);
    drawSetHisto("SeedLowestDphiHighPt","seedLowestDphi","foundSTAseed&&(mu_pt>10)",0,0.3,60);
    
    drawSetHisto("SeedLowestDphiDT","seedLowestDphi","foundSTAseed&&abs(mu_eta)<0.8",0,0.3,60);
    drawSetHisto("SeedLowestDphiLowPtDT","seedLowestDphi","foundSTAseed&&(mu_pt<5)&&abs(mu_eta)<0.8",0,0.3,60);
    drawSetHisto("SeedLowestDphiMediumPtDT","seedLowestDphi","foundSTAseed&&(mu_pt<10)&&(mu_pt>5)&&abs(mu_eta)<0.8",0,0.3,60);
    drawSetHisto("SeedLowestDphiHighPtDT","seedLowestDphi","foundSTAseed&&(mu_pt>10)&&abs(mu_eta)<0.8",0,0.3,60);
    
    drawSetHisto("SeedLowestDphiCSC","seedLowestDphi","foundSTAseed&&abs(mu_eta)>0.8",0,0.3,60);
    drawSetHisto("SeedLowestDphiLowPtCSC","seedLowestDphi","foundSTAseed&&(mu_pt<5)&&abs(mu_eta)>0.8",0,0.3,60);
    drawSetHisto("SeedLowestDphiMediumPtCSC","seedLowestDphi","foundSTAseed&&(mu_pt<10)&&(mu_pt>5)&&abs(mu_eta)>0.8",0,0.3,60);
    
    
    drawSetHisto("PtDistri","mu_pt","1",0,100,100);
    drawSetHisto("EtaDistri","mu_eta","1",-2.4,2.4,100);
    drawSetHisto("PtJpsiDistri","pair_pt","1",0,100,100);
    
    
    drawSetHisto("DphiCloseBy0","pair_dPhi","pair_pt>20&&pair_pt<30",0,0.5,50);
    drawSetHisto("DphiCloseBy1","pair_dPhi","pair_pt>30&&pair_pt<40",0,0.5,50);
    drawSetHisto("DphiCloseBy2","pair_dPhi","pair_pt>40&&pair_pt<50",0,0.5,50);
    drawSetHisto("DphiCloseBy3","pair_dPhi","pair_pt>50&&pair_pt<60",0,0.5,50);
    drawSetHisto("DphiCloseBy4","pair_dPhi","pair_pt>60&&pair_pt<70",0,0.5,50);
    drawSetHisto("DphiCloseBy5","pair_dPhi","pair_pt>70&&pair_pt<80",0,0.5,50);
    drawSetHisto("DphiCloseBy6","pair_dPhi","pair_pt>80&&pair_pt<90",0,0.5,50);


    
    myFile->Close();
    delete myFile;
    
}
