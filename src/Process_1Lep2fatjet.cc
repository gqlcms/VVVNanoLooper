#include "Process_1Lep2fatjet.h"
#include "TMath.h"

TLorentzVector getNeutrinoP4(double& MetPt, double& MetPhi, TLorentzVector& lep, int lepType){
        float MW_=80.385;

        double leppt = lep.Pt();
        double lepphi = lep.Phi();
        double lepeta = lep.Eta();
        double lepenergy = lep.Energy();

        double metpt = MetPt;
        double metphi = MetPhi;

        double  px = metpt*cos(metphi);
        double  py = metpt*sin(metphi);
        double  pz = 0;
        double  pxl= leppt*cos(lepphi);
        double  pyl= leppt*sin(lepphi);
        double  pzl= leppt*sinh(lepeta);
        double  El = lepenergy;
        double  a = pow(MW_,2) + pow(px+pxl,2) + pow(py+pyl,2) - px*px - py*py - El*El + pzl*pzl;
        double  b = 2.*pzl;
        double  A = b*b -4.*El*El;
        double  B = 2.*a*b;
        double  C = a*a-4.*(px*px+py*py)*El*El;

        ///////////////////////////pz for fnal
        double M_mu =  0;

        //if(lepType==1)M_mu=0.105658367;//mu
        //if(lepType==0)M_mu=0.00051099891;//electron

        int type=2; // use the small abs real root

        a = MW_*MW_ - M_mu*M_mu + 2.0*pxl*px + 2.0*pyl*py;
        A = 4.0*(El*El - pzl*pzl);
        B = -4.0*a*pzl;
        C = 4.0*El*El*(px*px + py*py) - a*a;

        double tmproot = B*B - 4.0*A*C;

        if (tmproot<0) {
            //std::cout << "Complex root detected, taking real part..." << std::endl;
            pz = - B/(2*A); // take real part of complex roots
        }
        else {
            double tmpsol1 = (-B + sqrt(tmproot))/(2.0*A);
            double tmpsol2 = (-B - sqrt(tmproot))/(2.0*A);
            //std::cout << " Neutrino Solutions: " << tmpsol1 << ", " << tmpsol2 << std::endl;

            if (type == 0 ) {
                // two real roots, pick the one closest to pz of muon
                if (TMath::Abs(tmpsol2-pzl) < TMath::Abs(tmpsol1-pzl)) { pz = tmpsol2; }
                else { pz = tmpsol1; }
                // if pz is > 300 pick the most central root
                if ( abs(pz) > 300. ) {
                    if (TMath::Abs(tmpsol1)<TMath::Abs(tmpsol2) ) { pz = tmpsol1; }
                    else { pz = tmpsol2; }
                }
            }
            if (type == 1 ) {
                // two real roots, pick the one closest to pz of muon
                if (TMath::Abs(tmpsol2-pzl) < TMath::Abs(tmpsol1-pzl)) { pz = tmpsol2; }
                else {pz = tmpsol1; }
            }
            if (type == 2 ) {
                // pick the most central root.
                if (TMath::Abs(tmpsol1)<TMath::Abs(tmpsol2) ) { pz = tmpsol1; }
                else { pz = tmpsol2; }
            }
            /*if (type == 3 ) {
             // pick the largest value of the cosine
             TVector3 p3w, p3mu;
             p3w.SetXYZ(pxl+px, pyl+py, pzl+ tmpsol1);
             p3mu.SetXYZ(pxl, pyl, pzl );

             double sinthcm1 = 2.*(p3mu.Perp(p3w))/MW_;
             p3w.SetXYZ(pxl+px, pyl+py, pzl+ tmpsol2);
             double sinthcm2 = 2.*(p3mu.Perp(p3w))/MW_;

             double costhcm1 = sqrt(1. - sinthcm1*sinthcm1);
             double costhcm2 = sqrt(1. - sinthcm2*sinthcm2);

             if ( costhcm1 > costhcm2 ) { pz = tmpsol1; otherSol_ = tmpsol2; }
             else { pz = tmpsol2;otherSol_ = tmpsol1; }

             }*///end of type3

        }//endl of if real root

        //dont correct pt neutrino
        TLorentzVector outP4;
        outP4.SetPxPyPzE(px,py,pz,sqrt(px*px+py*py+pz*pz));
        return outP4;

    }

void Process_1Lep2fatjet_test(){
    for (unsigned int imu = 0; imu < nt.Muon_p4().size(); ++imu){
        ana.tx.pushbackToBranch<float>("t_1Lep2fatjet_Muon_pt_test",nt.Muon_p4()[imu].pt());
        ana.tx.pushbackToBranch<float>("t_1Lep2fatjet_Muon_eta_test",nt.Muon_p4()[imu].eta());
        ana.tx.pushbackToBranch<float>("t_1Lep2fatjet_Muon_iso_test",nt.Muon_pfRelIso03_all()[imu]);
        ana.tx.pushbackToBranch<float>("t_1Lep2fatjet_Muon_iso2_test",nt.Muon_pfRelIso04_all()[imu]);
        ana.tx.pushbackToBranch<float>("t_1Lep2fatjet_Muon_iso3_test",nt.Muon_tkRelIso	()[imu]);
        
    }
    for (unsigned int imu = 0; imu < nt.Muon_p4().size(); ++imu)
    {
        int temp_id = nt.Muon_highPtId()[imu]; // have to switch to int
        // cout<<"nt.Muon_highPtId()[imu]:"<<temp_id<<" id origin:  "<<nt.Muon_highPtId()[imu]<<endl;
        if (not (temp_id==2 )) continue;
        if (not (nt.Muon_p4()[imu].pt()>50.0 )) continue;
        // std::cout<<"find one highPt muon"<<std::endl;
        ana.tx.pushbackToBranch<float>("t_1Lep2fatjet_highptMuon_pt_test",nt.Muon_p4()[imu].pt());
        ana.tx.pushbackToBranch<float>("t_1Lep2fatjet_highptMuon_eta_test",nt.Muon_p4()[imu].eta());
    }
    for (unsigned int usenumber1 = 0; usenumber1 < nt.FatJet_p4().size(); ++usenumber1)
    {
        ana.tx.pushbackToBranch<float>("t_1Lep2fatjet_fatjet_pt_test",nt.FatJet_p4()[usenumber1].pt());
        ana.tx.pushbackToBranch<float>("t_1Lep2fatjet_fatjet_eta_test",nt.FatJet_p4()[usenumber1].eta());
    }
}


void Process_1Lep2fatjet_Lepton(){
    // cout << "num of muon"  <<nt.Muon_p4().size()<<endl; 
    for (unsigned int imu = 0; imu < nt.Muon_p4().size(); ++imu)
    {
        int temp_id = nt.Muon_highPtId()[imu]; // have to switch to int
        // cout<<"nt.Muon_highPtId()[imu]:"<<temp_id<<" id origin:  "<<nt.Muon_highPtId()[imu]<<endl;
        if (not (temp_id==2 )) continue;
        if (not (nt.Muon_p4()[imu].pt()>20.0 )) continue;
        if (not (fabs(nt.Muon_p4()[imu].eta())<2.4 )) continue;
        if (not (nt.Muon_tkRelIso()[imu]<0.1 )) continue;
        // std::cout<<"find one highPt muon"<<std::endl;
        ana.tx.pushbackToBranch<int>("1Lep2fatjet_LooseMuon_highPtId_idxs",imu);
        ana.tx.pushbackToBranch<int>("t_1Lep2fatjet_LooseMuon_highPtId_idxs",imu);
    }
    for (unsigned int imu = 0; imu < nt.Muon_p4().size(); ++imu)
    {
        int temp_id = nt.Muon_highPtId()[imu]; // have to switch to int
        // cout<<"nt.Muon_highPtId()[imu]:"<<temp_id<<" id origin:  "<<nt.Muon_highPtId()[imu]<<endl;
        if (not (temp_id==2 )) continue;
        if (not (nt.Muon_p4()[imu].pt()>55.0 )) continue;
        if (not (fabs(nt.Muon_p4()[imu].eta())<2.4 )) continue;
        if (not (nt.Muon_tkRelIso()[imu]<0.1 )) continue;
        // std::cout<<"find one highPt muon"<<std::endl;
        ana.tx.pushbackToBranch<int>("1Lep2fatjet_CoodMuon_highPtId_idxs",imu);
        ana.tx.pushbackToBranch<int>("t_1Lep2fatjet_GoodMuon_highPtId_idxs",imu);
    }
    for (unsigned int iel = 0; iel < nt.Electron_p4().size(); ++iel)
    {
        if (not (nt.Electron_cutBased_HEEP()[iel])) continue;
        if (not (nt.Electron_p4()[iel].pt()>35.0 )) continue;
        if (not (fabs(nt.Electron_p4()[iel].eta())<2.4 )) continue;
        ana.tx.pushbackToBranch<int>("1Lep2fatjet_LooseElectron_cutBased_HEEP_idxs",iel);
        ana.tx.pushbackToBranch<int>("t_1Lep2fatjet_LooseElectron_cutBased_HEEP_idxs",iel);
    }
    for (unsigned int iel = 0; iel < nt.Electron_p4().size(); ++iel)
    {
        if (not (nt.Electron_cutBased_HEEP()[iel])) continue;
        if (not (nt.Electron_p4()[iel].pt()>55.0 )) continue;
        if (not (fabs(nt.Electron_p4()[iel].eta())<2.4 )) continue;
        ana.tx.pushbackToBranch<int>("1Lep2fatjet_GoodElectron_cutBased_HEEP_idxs",iel);
        ana.tx.pushbackToBranch<int>("t_1Lep2fatjet_GoodElectron_cutBased_HEEP_idxs",iel);
    }

}

float Process_1Lep2fatjet_Fatjet_SD_sub(Int_t FatJet_subJetIdx1,Int_t FatJet_subJetIdx2){
    TLorentzVector subjet1,subjet2,sum_p4;
    subjet1.SetPtEtaPhiE(nt.SubJet_p4()[FatJet_subJetIdx1].pt(),nt.SubJet_p4()[FatJet_subJetIdx1].eta(),nt.SubJet_p4()[FatJet_subJetIdx1].phi(),nt.SubJet_p4()[FatJet_subJetIdx1].energy());
    subjet2.SetPtEtaPhiE(nt.SubJet_p4()[FatJet_subJetIdx2].pt(),nt.SubJet_p4()[FatJet_subJetIdx2].eta(),nt.SubJet_p4()[FatJet_subJetIdx2].phi(),nt.SubJet_p4()[FatJet_subJetIdx2].energy());
    sum_p4 = subjet1+subjet2;
    return sum_p4.M();
}

float Process_1Lep2fatjet_Fatjet_SD_sub_no_JEC(Int_t FatJet_subJetIdx1,Int_t FatJet_subJetIdx2){
    TLorentzVector subjet1,subjet2,sum_p4;
    if(FatJet_subJetIdx1>=0&&FatJet_subJetIdx2>=0){
    float pt1 = nt.SubJet_p4()[FatJet_subJetIdx1].pt()*(1-nt.SubJet_rawFactor()[FatJet_subJetIdx1]);
    float pt2 = nt.SubJet_p4()[FatJet_subJetIdx2].pt()*(1-nt.SubJet_rawFactor()[FatJet_subJetIdx2]);
    float mass1 = nt.SubJet_p4()[FatJet_subJetIdx1].mass()*(1-nt.SubJet_rawFactor()[FatJet_subJetIdx1]);
    float mass2 = nt.SubJet_p4()[FatJet_subJetIdx2].mass()*(1-nt.SubJet_rawFactor()[FatJet_subJetIdx2]);
    subjet1.SetPtEtaPhiM(pt1,nt.SubJet_p4()[FatJet_subJetIdx1].eta(),nt.SubJet_p4()[FatJet_subJetIdx1].phi(),mass1);
    subjet2.SetPtEtaPhiM(pt2,nt.SubJet_p4()[FatJet_subJetIdx2].eta(),nt.SubJet_p4()[FatJet_subJetIdx2].phi(),mass2);
    sum_p4 = subjet1+subjet2;
    return sum_p4.M();
    }
    else{
        return -99998.0;
    }
}

void Process_1Lep2fatjet_Fatjet(){
    float fjSFvlc(1.), fjSFvlu(1.), fjSFvld(1.), fjSFlc(1.), fjSFlu(1.), fjSFld(1.), fjSFmc(1.), fjSFmu(1.), fjSFmd(1.), fjSFtc(1.), fjSFtu(1.), fjSFtd(1.);
    for (unsigned int usenumber1 = 0; usenumber1 < nt.FatJet_p4().size(); ++usenumber1)
    {
      float fjWPvloose = 0.274; //https://twiki.cern.ch/twiki/bin/view/CMS/DeepAK8Tagging2018WPsSFs
      float fjWPloose  = 0.506;
      float fjWPmedium = 0.731;
      float fjWPtight  = 0.828;
      if(nt.year() == 2017){
        fjWPvloose = 0.258;
        fjWPloose  = 0.506;
        fjWPmedium = 0.739;
        fjWPtight  = 0.838;
      }
      if(nt.year() == 2018){
        fjWPvloose = 0.245;
        fjWPloose  = 0.479;
        fjWPmedium = 0.704;
        fjWPtight  = 0.806;
      }

        // TODO: What is POG recommendation? do we use nt.FatJet_jetId()?
        // Figure this out
        // For now, accept anything above 250 GeV (TODO: is 250 GeV also ok?)
        if (not (nt.FatJet_p4()[usenumber1].pt() > 200.))
            continue;
        if (not (abs(nt.FatJet_p4()[usenumber1].eta()) < 2.4))
            continue;

        // Because every muon and electron shows up in PF FatJet collections
        // Need to check against leptons
        bool is_overlapping_with_a_lepton = false;

        // Overlap check against leptons (electrons)
        for (unsigned int ilep = 0; ilep < ana.tx.getBranchLazy<vector<int>>("1Lep2fatjet_GoodElectron_cutBased_HEEP_idxs").size(); ++ilep)
        {
            int ilep_idx = ana.tx.getBranchLazy<vector<int>>("1Lep2fatjet_GoodElectron_cutBased_HEEP_idxs")[ilep];
            if (RooUtil::Calc::DeltaR(nt.FatJet_p4()[usenumber1], nt.Electron_p4()[ilep_idx]) < 1.0)//SHOULD THIS BE 0.8? - try it
                {
                    is_overlapping_with_a_lepton = true;
                    cout<<"ele clean jet pt"<<nt.FatJet_p4()[usenumber1].pt()<<"ele pt "<<nt.Electron_p4()[ilep_idx].pt()<<endl;
                    break;
                }
        }

        for (unsigned int ilep = 0; ilep < ana.tx.getBranchLazy<vector<int>>("1Lep2fatjet_GoodMuon_highPtId_idxs").size(); ++ilep)
        {
            int ilep_idx = ana.tx.getBranchLazy<vector<int>>("1Lep2fatjet_GoodMuon_highPtId_idxs")[ilep];
            if (RooUtil::Calc::DeltaR(nt.FatJet_p4()[usenumber1], nt.Muon_p4()[ilep_idx]) < 1.0)//SHOULD THIS BE 0.8? - try it
                {
                    is_overlapping_with_a_lepton = true;
                    // cout<<"mu clean jet pt"<<nt.FatJet_p4()[usenumber1].pt()<<"  pt "<<nt.Muon_p4()[ilep_idx].pt()<<endl;
                    break;
                }
        }

        if (is_overlapping_with_a_lepton)
            continue;

 
        // For now, accept anything that reaches this point

    ana.tx.pushbackToBranch<int>("t_1Lep2fatjet_fatjet_idxs_new", usenumber1);
    ana.tx.pushbackToBranch<LorentzVector>("t_1Lep2fatjet_fatjet_p4_new", nt.FatJet_p4()[usenumber1]);
    }

    // leading pt jet
    int usenumber3 = -1; double pt_larger=0;
    for (unsigned int inum = 0; inum < ana.tx.getBranchLazy<vector<int>>("t_1Lep2fatjet_fatjet_idxs_new").size(); ++inum){
        if(nt.FatJet_p4()[ana.tx.getBranchLazy<vector<int>>("t_1Lep2fatjet_fatjet_idxs_new")[inum]].pt() > pt_larger && fabs(nt.FatJet_p4()[ana.tx.getBranchLazy<vector<int>>("t_1Lep2fatjet_fatjet_idxs_new")[inum]].eta())<2.4 && inum<4) {
            pt_larger = nt.FatJet_p4()[ana.tx.getBranchLazy<vector<int>>("t_1Lep2fatjet_fatjet_idxs_new")[inum]].pt(); 
            usenumber3 = ana.tx.getBranchLazy<vector<int>>("t_1Lep2fatjet_fatjet_idxs_new")[inum]; 
            continue;
        }
    }
    if (usenumber3>-1) {
    ana.tx.setBranch<float>("jetAK8puppi_ptJEC", nt.FatJet_p4()[usenumber3].pt());
    ana.tx.setBranch<float>("jetAK8puppi_eta", nt.FatJet_p4()[usenumber3].eta());
    ana.tx.setBranch<float>("jetAK8puppi_phi", nt.FatJet_p4()[usenumber3].phi());
    ana.tx.setBranch<int>("jetAK8puppi_tightid", nt.FatJet_jetId()[usenumber3]&2);
    
    ana.tx.setBranch<float>("jetAK8puppi_sd_JECcorr", nt.FatJet_msoftdrop()[usenumber3]);
    ana.tx.setBranch<float>("jetAK8puppi_mass", nt.FatJet_mass()[usenumber3]);
    // cout<<"sd_sub"<<Process_1Lep2fatjet_Fatjet_SD_sub(nt.FatJet_subJetIdx1()[usenumber3],nt.FatJet_subJetIdx2()[usenumber3])<<endl;
    // ana.tx.setBranch<float>("jetAK8puppi_sd_sub", Process_1Lep2fatjet_Fatjet_SD_sub(nt.FatJet_subJetIdx1()[usenumber3],nt.FatJet_subJetIdx2()[usenumber3]));
    ana.tx.setBranch<float>("jetAK8puppi_sd", Process_1Lep2fatjet_Fatjet_SD_sub_no_JEC(nt.FatJet_subJetIdx1()[usenumber3],nt.FatJet_subJetIdx2()[usenumber3]));
    ana.tx.setBranch<float>("Jet_rawFactor", nt.Jet_rawFactor()[usenumber3]);
    ana.tx.setBranch<float>("jetAK8puppi_dnnDecorrW", nt.FatJet_deepTagMD_WvsQCD()[usenumber3]);
    ana.tx.setBranch<float>("jetAK8puppi_dnnW", nt.FatJet_deepTag_WvsQCD()[usenumber3]);
    ana.tx.setBranch<float>("jetAK8puppi_dnnDecorrZ", nt.FatJet_deepTagMD_ZvsQCD()[usenumber3]);
    ana.tx.setBranch<float>("jetAK8puppi_dnnZ", nt.FatJet_deepTag_ZvsQCD()[usenumber3]);
    ana.tx.setBranch<float>("jetAK8puppi_dnnDecorrTop", nt.FatJet_deepTagMD_TvsQCD()[usenumber3]);
    ana.tx.setBranch<float>("jetAK8puppi_dnnTop", nt.FatJet_deepTag_TvsQCD()[usenumber3]);
    // ana.tx.pushbackToBranch<float>("Common_fatjet_deepMD_bb", nt.FatJet_deepTagMD_bbvsLight()[usenumber3]);
    }
    
    // second pt jet
    int usenumber2 = -1; 
    pt_larger=0;
    for (unsigned int inum = 0; inum < ana.tx.getBranchLazy<vector<int>>("t_1Lep2fatjet_fatjet_idxs_new").size(); ++inum){
        if(nt.FatJet_p4()[ana.tx.getBranchLazy<vector<int>>("t_1Lep2fatjet_fatjet_idxs_new")[inum]].pt() > pt_larger && fabs(nt.FatJet_p4()[ana.tx.getBranchLazy<vector<int>>("t_1Lep2fatjet_fatjet_idxs_new")[inum]].eta())<2.4 && inum<4 && ana.tx.getBranchLazy<vector<int>>("t_1Lep2fatjet_fatjet_idxs_new")[inum] != usenumber3) {
            pt_larger = nt.FatJet_p4()[ana.tx.getBranchLazy<vector<int>>("t_1Lep2fatjet_fatjet_idxs_new")[inum]].pt(); 
            usenumber2 = ana.tx.getBranchLazy<vector<int>>("t_1Lep2fatjet_fatjet_idxs_new")[inum]; 
            continue;
        }
    }
    if (usenumber2>-1) {
    ana.tx.setBranch<float>("jetAK8puppi_ptJEC_2", nt.FatJet_p4()[usenumber2].pt());
    ana.tx.setBranch<float>("jetAK8puppi_eta_2", nt.FatJet_p4()[usenumber2].eta());
    ana.tx.setBranch<float>("jetAK8puppi_phi_2", nt.FatJet_p4()[usenumber2].phi());
    ana.tx.setBranch<int>("jetAK8puppi_tightid_2", nt.FatJet_jetId()[usenumber2]&2);

    ana.tx.setBranch<float>("jetAK8puppi_sd_JECcorr_2", nt.FatJet_msoftdrop()[usenumber2]);
    ana.tx.setBranch<float>("jetAK8puppi_sd_2", Process_1Lep2fatjet_Fatjet_SD_sub_no_JEC(nt.FatJet_subJetIdx1()[usenumber2],nt.FatJet_subJetIdx2()[usenumber2]));
    ana.tx.setBranch<float>("Jet_rawFactor_2", nt.Jet_rawFactor()[usenumber2]);
    ana.tx.setBranch<float>("jetAK8puppi_dnnDecorrW_2", nt.FatJet_deepTagMD_WvsQCD()[usenumber2]);
    ana.tx.setBranch<float>("jetAK8puppi_dnnW_2", nt.FatJet_deepTag_WvsQCD()[usenumber2]);
    ana.tx.setBranch<float>("jetAK8puppi_dnnDecorrZ_2", nt.FatJet_deepTagMD_ZvsQCD()[usenumber2]);
    ana.tx.setBranch<float>("jetAK8puppi_dnnZ_2", nt.FatJet_deepTag_ZvsQCD()[usenumber2]);
    ana.tx.setBranch<float>("jetAK8puppi_dnnDecorrTop_2", nt.FatJet_deepTagMD_TvsQCD()[usenumber2]);
    ana.tx.setBranch<float>("jetAK8puppi_dnnTop_2", nt.FatJet_deepTag_TvsQCD()[usenumber2]);
    }

    // third pt jet
    int usenumber1 = -1; 
    pt_larger=0;
    for (unsigned int inum = 0; inum < ana.tx.getBranchLazy<vector<int>>("t_1Lep2fatjet_fatjet_idxs_new").size(); ++inum){
        if(nt.FatJet_p4()[ana.tx.getBranchLazy<vector<int>>("t_1Lep2fatjet_fatjet_idxs_new")[inum]].pt() > pt_larger && fabs(nt.FatJet_p4()[ana.tx.getBranchLazy<vector<int>>("t_1Lep2fatjet_fatjet_idxs_new")[inum]].eta())<2.4 && inum<4 && ana.tx.getBranchLazy<vector<int>>("t_1Lep2fatjet_fatjet_idxs_new")[inum] != usenumber3 && ana.tx.getBranchLazy<vector<int>>("t_1Lep2fatjet_fatjet_idxs_new")[inum] != usenumber2) {
            pt_larger = nt.FatJet_p4()[ana.tx.getBranchLazy<vector<int>>("t_1Lep2fatjet_fatjet_idxs_new")[inum]].pt(); 
            usenumber1 = ana.tx.getBranchLazy<vector<int>>("t_1Lep2fatjet_fatjet_idxs_new")[inum]; 
            continue;
        }
    }
    if (usenumber1>-1) {
    ana.tx.setBranch<float>("jetAK8puppi_ptJEC_3", nt.FatJet_p4()[usenumber1].pt());
    ana.tx.setBranch<float>("jetAK8puppi_eta_3", nt.FatJet_p4()[usenumber1].eta());
    ana.tx.setBranch<float>("jetAK8puppi_phi_3", nt.FatJet_p4()[usenumber1].phi());
    ana.tx.setBranch<int>("jetAK8puppi_tightid_3", nt.FatJet_jetId()[usenumber1]&2);

    ana.tx.setBranch<float>("jetAK8puppi_sd_JECcorr_3", nt.FatJet_msoftdrop()[usenumber1]);
    ana.tx.setBranch<float>("jetAK8puppi_sd_3", Process_1Lep2fatjet_Fatjet_SD_sub_no_JEC(nt.FatJet_subJetIdx1()[usenumber1],nt.FatJet_subJetIdx2()[usenumber1]));
    ana.tx.setBranch<float>("Jet_rawFactor_3", nt.Jet_rawFactor()[usenumber1]);
    ana.tx.setBranch<float>("jetAK8puppi_dnnDecorrW_3", nt.FatJet_deepTagMD_WvsQCD()[usenumber1]);
    ana.tx.setBranch<float>("jetAK8puppi_dnnW_3", nt.FatJet_deepTag_WvsQCD()[usenumber1]);
    ana.tx.setBranch<float>("jetAK8puppi_dnnDecorrZ_3", nt.FatJet_deepTagMD_ZvsQCD()[usenumber1]);
    ana.tx.setBranch<float>("jetAK8puppi_dnnZ_3", nt.FatJet_deepTag_ZvsQCD()[usenumber1]);
    ana.tx.setBranch<float>("jetAK8puppi_dnnDecorrTop_3", nt.FatJet_deepTagMD_TvsQCD()[usenumber1]);
    ana.tx.setBranch<float>("jetAK8puppi_dnnTop_3", nt.FatJet_deepTag_TvsQCD()[usenumber1]);
    }
}

void Process_1Lep2fatjet()
{
    //==============================================
    // Process_1Lep2fatjet:
    // This function gets called during the event looping.
    // This is where one sets the variables used for the category 1Lep2fatjet.
    //==============================================

    if (ana.test_1Lep2fatjet){
    Process_1Lep2fatjet_test();
    }

    Process_1Lep2fatjet_Lepton();
    Process_1Lep2fatjet_Fatjet();

    
    
    
    
    // MET et is the missing pt
    ana.tx.setBranch<float>("MET_et",nt.MET_pt() );
    ana.tx.setBranch<float>("MET_phi",nt.MET_phi() );
    
    
    // Met
    double MET_et,MET_phi;
    MET_et = nt.MET_pt();
    MET_phi = nt.MET_phi();
    ana.tx.setBranch<float>("MET_et", MET_et );
    ana.tx.setBranch<float>("MET_phi", MET_phi );


    // lepton
    
    float ptlep1, etalep1, philep1, energylep1;
    int lep;
    
    if ( ana.tx.getBranchLazy<vector<int>>("1Lep2fatjet_LooseElectron_cutBased_HEEP_idxs").size() == 1 ){
        ptlep1 = nt.Electron_p4()[ ana.tx.getBranchLazy<vector<int>>("1Lep2fatjet_LooseElectron_cutBased_HEEP_idxs")[0]].pt();
        etalep1 = nt.Electron_p4()[ana.tx.getBranchLazy<vector<int>>("1Lep2fatjet_LooseElectron_cutBased_HEEP_idxs")[0]].eta();
        philep1 = nt.Electron_p4()[ana.tx.getBranchLazy<vector<int>>("1Lep2fatjet_LooseElectron_cutBased_HEEP_idxs")[0]].phi();
        energylep1 = nt.Electron_p4()[ana.tx.getBranchLazy<vector<int>>("1Lep2fatjet_LooseElectron_cutBased_HEEP_idxs")[0]].energy();
        lep = nt.Electron_pdgId()[ana.tx.getBranchLazy<vector<int>>("1Lep2fatjet_LooseElectron_cutBased_HEEP_idxs")[0]];
    }
    if ( ana.tx.getBranchLazy<vector<int>>("1Lep2fatjet_LooseMuon_highPtId_idxs").size() == 1 ){
        ptlep1 = nt.Muon_p4()[ana.tx.getBranchLazy<vector<int>>("1Lep2fatjet_LooseMuon_highPtId_idxs")[0]].pt();
        etalep1 = nt.Muon_p4()[ana.tx.getBranchLazy<vector<int>>("1Lep2fatjet_LooseMuon_highPtId_idxs")[0]].eta();
        philep1 = nt.Muon_p4()[ana.tx.getBranchLazy<vector<int>>("1Lep2fatjet_LooseMuon_highPtId_idxs")[0]].phi();
        energylep1 = nt.Muon_p4()[ana.tx.getBranchLazy<vector<int>>("1Lep2fatjet_LooseMuon_highPtId_idxs")[0]].energy();
        lep = nt.Muon_pdgId()[ana.tx.getBranchLazy<vector<int>>("1Lep2fatjet_LooseMuon_highPtId_idxs")[0]];
    }
    ana.tx.setBranch<float>("ptlep1", ptlep1 );
    ana.tx.setBranch<float>("etalep1", etalep1 );
    ana.tx.setBranch<float>("philep1", philep1 );
    ana.tx.setBranch<float>("energylep1", energylep1 );
    ana.tx.setBranch<int>("lep", lep );


    // leptonic w
    TLorentzVector  glepton,neutrino,neutrinoP4,WLeptonic;
    glepton.SetPtEtaPhiE(ptlep1, etalep1, philep1, energylep1);
    int leptontype = 1;
    neutrino = getNeutrinoP4(MET_et, MET_phi, glepton, leptontype);
    neutrinoP4.SetPtEtaPhiE(neutrino.Pt(),neutrino.Eta(),neutrino.Phi(),neutrino.Energy());
    WLeptonic = glepton+neutrinoP4;
    float ptVlepJEC, yVlepJEC, phiVlepJEC,massVlepJEC,mtVlepJEC;
    ptVlepJEC    = WLeptonic.Pt();
    yVlepJEC     = WLeptonic.Eta();
    phiVlepJEC   = WLeptonic.Phi();
    massVlepJEC  = WLeptonic.M();
    mtVlepJEC    = WLeptonic.Mt();
    ana.tx.setBranch<float>("ptVlepJEC", ptVlepJEC );
    ana.tx.setBranch<float>("yVlepJEC", yVlepJEC );
    ana.tx.setBranch<float>("phiVlepJEC", phiVlepJEC );
    ana.tx.setBranch<float>("massVlepJEC", massVlepJEC );
    ana.tx.setBranch<float>("mtVlepJEC", mtVlepJEC );

    // neutrino
    float ptlep2, etalep2, philep2, energylep2;
    ptlep2 = neutrinoP4.Pt();
    etalep2 = neutrinoP4.Eta();
    philep2 = neutrinoP4.Phi();
    energylep2 = neutrinoP4.Energy();
    ana.tx.setBranch<float>("ptlep2", ptlep2 );
    ana.tx.setBranch<float>("etalep2", etalep2 );
    ana.tx.setBranch<float>("philep2", philep2 );
    ana.tx.setBranch<float>("energylep2", energylep2 );



    
    for (unsigned int ijet = 0; ijet < nt.Jet_p4().size(); ++ijet)
    {
      
        //  do we use nt.Jet_jetId()?
        //   Jet ID flags bit1 is loose (always false in 2017 since it does not exist), bit2 is tight, bit3 is tightLepVeto

        if ( not (nt.Jet_p4()[ijet].pt() > 20.) ){
            int _Jet_jetId = nt.Jet_jetId()[ijet];
            if (_Jet_jetId !=2)
                continue;
        }
        if (not (abs(nt.Jet_p4()[ijet].eta()) < 5.0))
            continue;

        // Because every muon and electron shows up in PF FatJet collections
        // Need to check against leptons
        bool is_overlapping_with_a_lepton = false;

        // Overlap check against leptons (electrons)
        for (unsigned int ilep = 0; ilep < ana.tx.getBranchLazy<vector<int>>("1Lep2fatjet_LooseElectron_cutBased_HEEP_idxs").size(); ++ilep)
        {
            int ilep_idx = ana.tx.getBranchLazy<vector<int>>("1Lep2fatjet_LooseElectron_cutBased_HEEP_idxs")[ilep];
            if (RooUtil::Calc::DeltaR(nt.Jet_p4()[ijet], nt.Electron_p4()[ilep_idx]) < 1.0)//SHOULD THIS BE 0.8? - try it
                {
                    is_overlapping_with_a_lepton = true;
                    // cout<<"ele clean jet pt"<<nt.FatJet_p4()[usenumber1].pt()<<"ele pt "<<nt.Electron_p4()[ilep_idx].pt()<<endl;
                    break;
                }
        }

        for (unsigned int ilep = 0; ilep < ana.tx.getBranchLazy<vector<int>>("1Lep2fatjet_LooseMuon_highPtId_idxs").size(); ++ilep)
        {
            int ilep_idx = ana.tx.getBranchLazy<vector<int>>("1Lep2fatjet_LooseMuon_highPtId_idxs")[ilep];
            if (RooUtil::Calc::DeltaR(nt.Jet_p4()[ijet], nt.Muon_p4()[ilep_idx]) < 1.0)//SHOULD THIS BE 0.8? - try it
                {
                    is_overlapping_with_a_lepton = true;
                    // cout<<"mu clean jet pt"<<nt.FatJet_p4()[usenumber1].pt()<<"  pt "<<nt.Muon_p4()[ilep_idx].pt()<<endl;
                    break;
                }
        }

        if (is_overlapping_with_a_lepton)
            continue;

 
        // For now, accept anything that reaches this point

    ana.tx.pushbackToBranch<int>("Common_jet_idxs_new", ijet);
    ana.tx.pushbackToBranch<LorentzVector>("Common_jet_p4_new", nt.Jet_p4()[ijet]);
    ana.tx.pushbackToBranch<int>("ak4jet_hf", nt.Jet_hadronFlavour()[ijet]);
    ana.tx.pushbackToBranch<int>("ak4jet_pf", nt.Jet_partonFlavour()[ijet]);
    ana.tx.pushbackToBranch<float>("ak4jet_pt", nt.Jet_p4()[ijet].pt());
    ana.tx.pushbackToBranch<float>("ak4jet_eta", nt.Jet_p4()[ijet].eta());
    ana.tx.pushbackToBranch<float>("ak4jet_phi", nt.Jet_p4()[ijet].phi());
    ana.tx.pushbackToBranch<float>("ak4jet_e", nt.Jet_p4()[ijet].energy());
    ana.tx.pushbackToBranch<float>("Jet_btagDeepB", nt.Jet_btagDeepB()[ijet]);
    ana.tx.pushbackToBranch<float>("Jet_btagDeepC", nt.Jet_btagDeepC()[ijet]);

    }


    // Gen level info
    
    if (not nt.isData())
    {
        // Z boson info, last copy
        for(size_t ik=0; ik<nt.nGenPart();ik++)
        {
            if (abs(nt.GenPart_pdgId()[ik]) == 23 )
            {
                if (not (nt.GenPart_statusFlags()[ik]&(1<<13))) continue; // isLastCopy
                int indexgenzl;
                indexgenzl = ik;
                float ptgenzl,etagenzl,phigenzl,massgenzl;
                ptgenzl = nt.GenPart_pt()[ik];
                etagenzl = nt.GenPart_eta()[ik];
                phigenzl = nt.GenPart_phi()[ik];
                massgenzl = nt.GenPart_p4()[ik].energy();
                ana.tx.pushbackToBranch<float>("ptgenzl",ptgenzl);
                ana.tx.pushbackToBranch<float>("etagenzl",etagenzl);
                ana.tx.pushbackToBranch<float>("phigenzl",phigenzl);
                ana.tx.pushbackToBranch<float>("massgenzl",massgenzl);
                vector<float> daughter_index;
                for (size_t id=0; id<nt.nGenPart();id++){
                    // if (nt.GenPart_genPartIdxMother()[id].size() >0){
                    //     for (size_t im=0; im<nt.GenPart_genPartIdxMother()[id].size();im++){
                    //         if (nt.GenPart_genPartIdxMother()[id][im] == indexgenzl){
                    //             daughter_index.push_back(id);
                    //         }

                    //     }
                    // }
                    if (nt.GenPart_genPartIdxMother()[id] == indexgenzl){
                        daughter_index.push_back(id);
                    }
                }

                float genz_q1_pt,genz_q1_eta,genz_q1_phi,genz_q1_e,genz_q1_pdg,genz_q2_pt,genz_q2_eta,genz_q2_phi,genz_q2_e,genz_q2_pdg;
                genz_q1_pt = nt.GenPart_pt()[daughter_index[0]];
                genz_q1_eta = nt.GenPart_eta()[daughter_index[0]];
                genz_q1_phi = nt.GenPart_phi()[daughter_index[0]];
                genz_q1_e = nt.GenPart_p4()[daughter_index[0]].energy();
                genz_q1_pdg = nt.GenPart_pdgId()[daughter_index[0]];
                genz_q2_pt = nt.GenPart_pt()[daughter_index[1]];
                genz_q2_eta = nt.GenPart_eta()[daughter_index[1]];
                genz_q2_phi = nt.GenPart_phi()[daughter_index[1]];
                genz_q2_e = nt.GenPart_p4()[daughter_index[1]].energy();
                genz_q2_pdg = nt.GenPart_pdgId()[daughter_index[1]];

                ana.tx.pushbackToBranch<float>("genz_q1_pt",genz_q1_pt);
                ana.tx.pushbackToBranch<float>("genz_q1_eta",genz_q1_eta);
                ana.tx.pushbackToBranch<float>("genz_q1_phi",genz_q1_phi);
                ana.tx.pushbackToBranch<float>("genz_q1_e",genz_q1_e);
                ana.tx.pushbackToBranch<float>("genz_q1_pdg",genz_q1_pdg);
                ana.tx.pushbackToBranch<float>("genz_q2_pt",genz_q2_pt);
                ana.tx.pushbackToBranch<float>("genz_q2_eta",genz_q2_eta);
                ana.tx.pushbackToBranch<float>("genz_q2_phi",genz_q2_phi);
                ana.tx.pushbackToBranch<float>("genz_q2_e",genz_q2_e);
                ana.tx.pushbackToBranch<float>("genz_q2_pdg",genz_q2_pdg);
               
            }
        }

        // W boson info, first copy
        for(size_t ik=0; ik<nt.nGenPart();ik++)
        {
            if (abs(nt.GenPart_pdgId()[ik]) == 24 )
            {
                if (not (nt.GenPart_statusFlags()[ik]&(1<<13))) continue; // isLastCopy
                int indexgenwl;
                indexgenwl = ik;
                float ptgenwl,etagenwl,phigenwl,massgenwl;
                ptgenwl = nt.GenPart_pt()[ik];
                etagenwl = nt.GenPart_eta()[ik];
                phigenwl = nt.GenPart_phi()[ik];
                massgenwl = nt.GenPart_p4()[ik].energy();
                ana.tx.pushbackToBranch<float>("ptgenwl",ptgenwl);
                ana.tx.pushbackToBranch<float>("etagenwl",etagenwl);
                ana.tx.pushbackToBranch<float>("phigenwl",phigenwl);
                ana.tx.pushbackToBranch<float>("massgenwl",massgenwl);
                vector<float> daughter_index;
                for (size_t id=0; id<nt.nGenPart();id++){
                    // if (nt.GenPart_genPartIdxMother()[id].size() >0){
                    //     for (size_t im=0; im<nt.GenPart_genPartIdxMother()[id].size();im++){
                    //         if (nt.GenPart_genPartIdxMother()[id][im] == indexgenwl){
                    //             daughter_index.push_back(id);
                    //         }

                    //     }
                    // }
                    if (nt.GenPart_genPartIdxMother()[id] == indexgenwl){
                        daughter_index.push_back(id);
                    }
                }

                float genw_q1_pt,genw_q1_eta,genw_q1_phi,genw_q1_e,genw_q1_pdg,genw_q2_pt,genw_q2_eta,genw_q2_phi,genw_q2_e,genw_q2_pdg;
                // cout<<"gen info"<< nt.GenPart_pt()[daughter_index[0]]<<endl;
                // genw_q1_pt = nt.GenPart_pt()[daughter_index[0]];
                genw_q1_pt = nt.GenPart_p4()[daughter_index[0]].pt();
                genw_q1_eta = nt.GenPart_eta()[daughter_index[0]];
                genw_q1_phi = nt.GenPart_phi()[daughter_index[0]];
                genw_q1_e = nt.GenPart_p4()[daughter_index[0]].energy();
                genw_q1_pdg = nt.GenPart_pdgId()[daughter_index[0]];
                genw_q2_pt = nt.GenPart_pt()[daughter_index[1]];
                genw_q2_eta = nt.GenPart_eta()[daughter_index[1]];
                genw_q2_phi = nt.GenPart_phi()[daughter_index[1]];
                genw_q2_e = nt.GenPart_p4()[daughter_index[1]].energy();
                genw_q2_pdg = nt.GenPart_pdgId()[daughter_index[1]];

                ana.tx.pushbackToBranch<float>("genw_q1_pt",genw_q1_pt);
                ana.tx.pushbackToBranch<float>("genw_q1_eta",genw_q1_eta);
                ana.tx.pushbackToBranch<float>("genw_q1_phi",genw_q1_phi);
                ana.tx.pushbackToBranch<float>("genw_q1_e",genw_q1_e);
                ana.tx.pushbackToBranch<float>("genw_q1_pdg",genw_q1_pdg);
                ana.tx.pushbackToBranch<float>("genw_q2_pt",genw_q2_pt);
                ana.tx.pushbackToBranch<float>("genw_q2_eta",genw_q2_eta);
                ana.tx.pushbackToBranch<float>("genw_q2_phi",genw_q2_phi);
                ana.tx.pushbackToBranch<float>("genw_q2_e",genw_q2_e);
                ana.tx.pushbackToBranch<float>("genw_q2_pdg",genw_q2_pdg);
            }
        }


        // top info


    }
    


    
    ana.tx.fill();
    // Set variables used in this category.
    // If histograms are booked with these variables the histograms will be filled automatically.
    // Please follow the convention of <category>_<varname> structure.

    // Example of reading from Nano
    // std::vector<LorentzVector> electron_p4s = nt.Electron_p4(); // nt is a global variable that accesses NanoAOD
    // std::vector<float> electron_mvaTTH = nt.Electron_mvaTTH(); // electron ttH MVA scores from NanoAOD
    // Semi-complete list of NanoAOD for 102X can be found here: https://cms-nanoaod-integration.web.cern.ch/integration/master-102X/mc102X_doc.html
    // Also consult here: https://github.com/cmstas/NanoTools/blob/d641a6d6c1aa9ecc8094a1af73d5e1bd7d6502ab/NanoCORE/Nano.h#L4875 (if new variables are added they may show up in master)
}
