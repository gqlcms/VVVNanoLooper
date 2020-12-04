#include "Begin_1Lep2fatjet.h"

void Begin_1Lep2fatjet_test(){
    ana.tx.createBranch<vector<int>>("t_1Lep2fatjet_Muon_idxs_test");
    ana.tx.createBranch<vector<float>>("t_1Lep2fatjet_Muon_pt_test");
    ana.tx.createBranch<vector<float>>("t_1Lep2fatjet_Muon_eta_test");
    ana.tx.createBranch<vector<float>>("t_1Lep2fatjet_Muon_iso_test");
    ana.tx.createBranch<vector<float>>("t_1Lep2fatjet_Muon_iso2_test");
    ana.tx.createBranch<vector<float>>("t_1Lep2fatjet_Muon_iso3_test");
    ana.tx.createBranch<vector<float>>("t_1Lep2fatjet_highptMuon_pt_test");
    ana.tx.createBranch<vector<float>>("t_1Lep2fatjet_highptMuon_eta_test");

    ana.tx.createBranch<vector<float>>("t_1Lep2fatjet_fatjet_pt_test");
    ana.tx.createBranch<vector<float>>("t_1Lep2fatjet_fatjet_eta_test");



}

void Begin_1Lep2fatjet()
{
    //==============================================
    // Begin_1Lep2fatjet:
    // This function gets called prior to the event looping.
    // This is where one declares variables, histograms, event selections for the category 1Lep2fatjet.
    //==============================================


    // Create variables used in this category.
    // Please follow the convention of <category>_<varname> structure.

    if (ana.test_1Lep2fatjet){
        Begin_1Lep2fatjet_test();
    }

    ana.tx.createBranch<LorentzVector>("1Lep2fatjet_lep_p4");
    ana.tx.createBranch<int>("1Lep2fatjet_lep_pdgid");
    ana.tx.createBranch<vector<int>>("1Lep2fatjet_LooseMuon_highPtId_idxs");
    ana.tx.createBranch<vector<int>>("1Lep2fatjet_LooseElectron_cutBased_HEEP_idxs");

    
    ana.tx.createBranch<vector<int>>("t_1Lep2fatjet_LooseMuon_highPtId_idxs");
    ana.tx.createBranch<vector<int>>("t_1Lep2fatjet_LooseElectron_cutBased_HEEP_idxs");

    ana.tx.createBranch<vector<int>>("1Lep2fatjet_GoodMuon_highPtId_idxs");
    ana.tx.createBranch<vector<int>>("1Lep2fatjet_GoodElectron_cutBased_HEEP_idxs");

    
    ana.tx.createBranch<vector<int>>("t_1Lep2fatjet_GoodMuon_highPtId_idxs");
    ana.tx.createBranch<vector<int>>("t_1Lep2fatjet_GoodElectron_cutBased_HEEP_idxs");

    ana.tx.createBranch<vector<int>>("t_1Lep2fatjet_fatjet_idxs_new");
    ana.tx.createBranch<vector<LorentzVector>>("t_1Lep2fatjet_fatjet_p4_new");
    ana.tx.createBranch<float>("Jet_rawFactor");

    ana.tx.createBranch<float>("jetAK8puppi_ptJEC");
    ana.tx.createBranch<float>("jetAK8puppi_eta");
    ana.tx.createBranch<float>("jetAK8puppi_phi");
    ana.tx.createBranch<int>("jetAK8puppi_tightid");

    ana.tx.createBranch<float>("jetAK8puppi_ptJEC_2");
    ana.tx.createBranch<float>("jetAK8puppi_eta_2");
    ana.tx.createBranch<float>("jetAK8puppi_phi_2");
    ana.tx.createBranch<int>("jetAK8puppi_tightid_2");

    ana.tx.createBranch<float>("jetAK8puppi_ptJEC_3");
    ana.tx.createBranch<float>("jetAK8puppi_eta_3");
    ana.tx.createBranch<float>("jetAK8puppi_phi_3");
    ana.tx.createBranch<int>("jetAK8puppi_tightid_3");

    ana.tx.createBranch<float>("jetAK8puppi_sd");
    ana.tx.createBranch<float>("jetAK8puppi_sd_JECcorr");
    // ana.tx.createBranch<float>("jetAK8puppi_sd_sub");
    ana.tx.createBranch<float>("jetAK8puppi_mass");
    ana.tx.createBranch<float>("jetAK8puppi_dnnDecorrW");
    ana.tx.createBranch<float>("jetAK8puppi_dnnW");
    ana.tx.createBranch<float>("jetAK8puppi_dnnDecorrZ");
    ana.tx.createBranch<float>("jetAK8puppi_dnnZ");
    ana.tx.createBranch<float>("jetAK8puppi_dnnDecorrTop");
    ana.tx.createBranch<float>("jetAK8puppi_dnnTop");
    
    
    ana.tx.createBranch<float>("jetAK8puppi_sd_JECcorr_2");
    ana.tx.createBranch<float>("jetAK8puppi_sd_2");
    ana.tx.createBranch<float>("Jet_rawFactor_2");
    ana.tx.createBranch<float>("jetAK8puppi_dnnDecorrW_2");
    ana.tx.createBranch<float>("jetAK8puppi_dnnW_2");
    ana.tx.createBranch<float>("jetAK8puppi_dnnDecorrZ_2");
    ana.tx.createBranch<float>("jetAK8puppi_dnnZ_2");
    ana.tx.createBranch<float>("jetAK8puppi_dnnDecorrTop_2");
    ana.tx.createBranch<float>("jetAK8puppi_dnnTop_2");

    ana.tx.createBranch<float>("jetAK8puppi_sd_JECcorr_3");
    ana.tx.createBranch<float>("jetAK8puppi_sd_3");
    ana.tx.createBranch<float>("Jet_rawFactor_3");
    ana.tx.createBranch<float>("jetAK8puppi_dnnDecorrW_3");
    ana.tx.createBranch<float>("jetAK8puppi_dnnW_3");
    ana.tx.createBranch<float>("jetAK8puppi_dnnDecorrZ_3");
    ana.tx.createBranch<float>("jetAK8puppi_dnnZ_3");
    ana.tx.createBranch<float>("jetAK8puppi_dnnDecorrTop_3");
    ana.tx.createBranch<float>("jetAK8puppi_dnnTop_3");

    ana.tx.createBranch<vector<int>>("Common_jet_idxs_new");
    ana.tx.createBranch<vector<LorentzVector>>("Common_jet_p4_new");
    ana.tx.createBranch<vector<int>>("ak4jet_hf");
    ana.tx.createBranch<vector<int>>("ak4jet_pf");
    ana.tx.createBranch<vector<float>>("ak4jet_pt");
    ana.tx.createBranch<vector<float>>("ak4jet_eta");
    ana.tx.createBranch<vector<float>>("ak4jet_phi");
    ana.tx.createBranch<vector<float>>("ak4jet_e");
    ana.tx.createBranch<vector<float>>("Jet_btagDeepB");
    ana.tx.createBranch<vector<float>>("Jet_btagDeepC");


    ana.tx.createBranch<float>("ptlep1");
    ana.tx.createBranch<float>("etalep1");
    ana.tx.createBranch<float>("philep1");
    ana.tx.createBranch<float>("energylep1");
    ana.tx.createBranch<int>("lep");

    ana.tx.createBranch<float>("ptlep2");
    ana.tx.createBranch<float>("etalep2" );
    ana.tx.createBranch<float>("philep2" );
    ana.tx.createBranch<float>("energylep2" );

    ana.tx.createBranch<float>("ptVlepJEC" );
    ana.tx.createBranch<float>("yVlepJEC" );
    ana.tx.createBranch<float>("phiVlepJEC" );
    ana.tx.createBranch<float>("massVlepJEC" );
    ana.tx.createBranch<float>("mtVlepJEC" );


    // gen level info
    ana.tx.createBranch<vector<float>>("ptgenwl");
    ana.tx.createBranch<vector<float>>("etagenwl");
    ana.tx.createBranch<vector<float>>("phigenwl");
    ana.tx.createBranch<vector<float>>("massgenwl");
    ana.tx.createBranch<vector<float>>("ptgenzl");
    ana.tx.createBranch<vector<float>>("etagenzl");
    ana.tx.createBranch<vector<float>>("phigenzl");
    ana.tx.createBranch<vector<float>>("massgenzl");

    ana.tx.createBranch<vector<float>>("genz_q1_pt");
    ana.tx.createBranch<vector<float>>("genz_q1_eta");
    ana.tx.createBranch<vector<float>>("genz_q1_phi");
    ana.tx.createBranch<vector<float>>("genz_q1_e");
    ana.tx.createBranch<vector<float>>("genz_q1_pdg");
    ana.tx.createBranch<vector<float>>("genz_q2_pt");
    ana.tx.createBranch<vector<float>>("genz_q2_eta");
    ana.tx.createBranch<vector<float>>("genz_q2_phi");
    ana.tx.createBranch<vector<float>>("genz_q2_e");
    ana.tx.createBranch<vector<float>>("genz_q2_pdg");

    ana.tx.createBranch<vector<float>>("genw_q1_pt");
    ana.tx.createBranch<vector<float>>("genw_q1_eta");
    ana.tx.createBranch<vector<float>>("genw_q1_phi");
    ana.tx.createBranch<vector<float>>("genw_q1_e");
    ana.tx.createBranch<vector<float>>("genw_q1_pdg");
    ana.tx.createBranch<vector<float>>("genw_q2_pt");
    ana.tx.createBranch<vector<float>>("genw_q2_eta");
    ana.tx.createBranch<vector<float>>("genw_q2_phi");
    ana.tx.createBranch<vector<float>>("genw_q2_e");
    ana.tx.createBranch<vector<float>>("genw_q2_pdg");


    
    

    ana.tx.createBranch<float>("MET_et");
    ana.tx.createBranch<float>("MET_phi");

    
    

    ana.cutflow.getCut("CommonCut");
    ana.cutflow.addCutToLastActiveCut("Cut_1Lep2fatjet_Has1Lepton_new", [&](){
        if (not ((ana.tx.getBranchLazy<vector<int>>("1Lep2fatjet_LooseMuon_highPtId_idxs").size()+ana.tx.getBranchLazy<vector<int>>("1Lep2fatjet_LooseElectron_cutBased_HEEP_idxs").size()) == 1))
           return false;
        return true;

    }, UNITY);
    ana.cutflow.addCutToLastActiveCut("Cut_1Lep2fatjet_Has2FatJet_new", [&]() { return (ana.tx.getBranchLazy<vector<int>>("Common_fatjet_idxs_new").size() >= 2);}, UNITY);




    // Define selections
    // CommonCut will contain selections that should be common to all categories, starting from this cut, add cuts for this category of the analysis.
    ana.cutflow.getCut("CommonCut");
    // ana.cutflow.addCutToLastActiveCut("Cut_1Lep2fatjet_GenHT", [&]() { return ana.tx.getBranch<float>("Common_genHT"); }, UNITY);
    ana.cutflow.addCutToLastActiveCut("Cut_1Lep2fatjet_Has1Lepton", [&]()
            {
                if (not (ana.tx.getBranchLazy<vector<int>>("Common_lep_idxs").size() == 1))
                    return false;
                ana.tx.setBranch<LorentzVector>("1Lep2fatjet_lep_p4", ana.tx.getBranchLazy<vector<LorentzVector>>("Common_lep_p4")[0]);
                ana.tx.setBranch<int>("1Lep2fatjet_lep_pdgid", ana.tx.getBranchLazy<vector<int>>("Common_lep_pdgid")[0]);
                return true;
            }, UNITY);
    ana.cutflow.addCutToLastActiveCut("Cut_1Lep2fatjet_Has2FatJet", [&]() { return (ana.tx.getBranchLazy<vector<int>>("Common_fatjet_idxs").size() >= 2);}, UNITY);
    ana.cutflow.addCutToLastActiveCut("Cut_1Lep2fatjet_nbVeto", [&]() { return (ana.tx.getBranch<int>("Common_nb_tight") == 0);}, UNITY);
    ana.cutflow.addCutToLastActiveCut("Cut_1Lep2fatjet_FatJetMassPresel", [&]() { return (nt.FatJet_msoftdrop()[ana.tx.getBranch<vector<int>>("Common_fatjet_idxs")[0]] > 20.) and (nt.FatJet_msoftdrop()[ana.tx.getBranch<vector<int>>("Common_fatjet_idxs")[1]] > 20.);}, UNITY);
    ana.cutflow.addCutToLastActiveCut("Cut_1Lep2fatjet_DeepWProd", [&]() { return (nt.FatJet_deepTag_WvsQCD()[ana.tx.getBranch<vector<int>>("Common_fatjet_idxs")[0]] * nt.FatJet_deepTag_WvsQCD()[ana.tx.getBranch<vector<int>>("Common_fatjet_idxs")[1]]) > 0.5;}, UNITY);
    ana.cutflow.addCutToLastActiveCut("Cut_1Lep2fatjet_Preselection", UNITY, UNITY);

    // Create histograms used in this category.
    // Please follow the convention of h_<category>_<varname> structure.
    // N.B. Using nbins of size 180 or 360 can provide flexibility as it can be rebinned easily, as 180, 360 are highly composite numbers.
    RooUtil::Histograms hists_1Lep2fatjet_new;
        
    hists_1Lep2fatjet_new.addHistogram("h_1Lep2fatjet_new_event_id", 100, 0, 100, [&]() { return ana.tx.getBranch<unsigned long long>("Common_evt"); } );
    hists_1Lep2fatjet_new.addHistogram("h_1Lep2fatjet_new_fatjet_msoftdrop_0", 180, 0, 250, [&]() { return nt.FatJet_msoftdrop()[ana.tx.getBranch<vector<int>>("Common_fatjet_idxs_new")[0]]; } );

    RooUtil::Histograms hists_1Lep2fatjet;
//    hists_1Lep2fatjet.addHistogram("h_1Lep2fatjet_lep_pt", 180, 0, 500, [&]() { return ana.tx.getBranch<LorentzVector>("1Lep2fatjet_lep_p4").pt(); } );
    hists_1Lep2fatjet.addHistogram("h_1Lep2fatjet_event_id", 100, 0, 100, [&]() { return ana.tx.getBranch<unsigned long long>("Common_evt"); } );
//    hists_1Lep2fatjet.addHistogram("h_1Lep2fatjet_lep_pdgid", 4, 10, 14, [&]() { return abs(ana.tx.getBranch<int>("1Lep2fatjet_lep_pdgid")); } );
//    hists_1Lep2fatjet.addHistogram("h_1Lep2fatjet_lep_eta", 180, -3, 3, [&]() { return ana.tx.getBranch<LorentzVector>("1Lep2fatjet_lep_p4").eta(); } );
//    hists_1Lep2fatjet.addHistogram("h_1Lep2fatjet_mt", 180, 0, 250, [&]() { return RooUtil::Calc::mT(ana.tx.getBranch<LorentzVector>("1Lep2fatjet_lep_p4"), ana.tx.getBranch<LorentzVector>("Common_met_p4")); } );
//    hists_1Lep2fatjet.addHistogram("h_1Lep2fatjet_pt_lv", 180, 0, 1000, [&]() { return (ana.tx.getBranch<LorentzVector>("1Lep2fatjet_lep_p4") + ana.tx.getBranch<LorentzVector>("Common_met_p4")).pt(); } );
    hists_1Lep2fatjet.addHistogram("h_1Lep2fatjet_met_pt", 180, 0, 500, [&]() { return ana.tx.getBranch<LorentzVector>("Common_met_p4").pt(); } );
    hists_1Lep2fatjet.addHistogram("h_1Lep2fatjet_met_phi", 180, -3.1416, 3.1416, [&]() { return ana.tx.getBranch<LorentzVector>("Common_met_p4").phi(); } );
    hists_1Lep2fatjet.addHistogram("h_1Lep2fatjet_nb_loose", 8, 0, 8, [&]() { return ana.tx.getBranch<int>("Common_nb_loose"); } );
    hists_1Lep2fatjet.addHistogram("h_1Lep2fatjet_nb_medium", 8, 0, 8, [&]() { return ana.tx.getBranch<int>("Common_nb_medium"); } );
    hists_1Lep2fatjet.addHistogram("h_1Lep2fatjet_nb_tight", 8, 0, 8, [&]() { return ana.tx.getBranch<int>("Common_nb_tight"); } );
    hists_1Lep2fatjet.addHistogram("h_1Lep2fatjet_njet", 15, 0, 15, [&]() { return ana.tx.getBranch<vector<int>>("Common_jet_idxs").size(); } );
    hists_1Lep2fatjet.addHistogram("h_1Lep2fatjet_nfatjet", 5, 0, 5, [&]() { return ana.tx.getBranch<vector<int>>("Common_fatjet_idxs").size(); } );
    hists_1Lep2fatjet.addHistogram("h_1Lep2fatjet_difatjet_mass", 180, 0, 3000, [&]() { return (ana.tx.getBranch<vector<LorentzVector>>("Common_fatjet_p4")[0] + ana.tx.getBranch<vector<LorentzVector>>("Common_fatjet_p4")[1]).mass(); } );
    hists_1Lep2fatjet.addHistogram("h_1Lep2fatjet_difatjet_pt", 180, 0, 1000, [&]() { return (ana.tx.getBranch<vector<LorentzVector>>("Common_fatjet_p4")[0] + ana.tx.getBranch<vector<LorentzVector>>("Common_fatjet_p4")[1]).pt(); } );
    hists_1Lep2fatjet.addHistogram("h_1Lep2fatjet_fatjet_eta_0", 180, -5, 5, [&]() { return ana.tx.getBranch<vector<LorentzVector>>("Common_fatjet_p4")[0].eta(); } );
    hists_1Lep2fatjet.addHistogram("h_1Lep2fatjet_fatjet_eta_1", 180, -5, 5, [&]() { return ana.tx.getBranch<vector<LorentzVector>>("Common_fatjet_p4")[1].eta(); } );
    hists_1Lep2fatjet.addHistogram("h_1Lep2fatjet_fatjet_pt_0", 180, 0, 1000, [&]() { return ana.tx.getBranch<vector<LorentzVector>>("Common_fatjet_p4")[0].pt(); } );
    hists_1Lep2fatjet.addHistogram("h_1Lep2fatjet_fatjet_pt_1", 180, 0, 750, [&]() { return ana.tx.getBranch<vector<LorentzVector>>("Common_fatjet_p4")[1].pt(); } );
    hists_1Lep2fatjet.addHistogram("h_1Lep2fatjet_fatjet_msoftdrop_0", 180, 0, 250, [&]() { return nt.FatJet_msoftdrop()[ana.tx.getBranch<vector<int>>("Common_fatjet_idxs")[0]]; } );
    hists_1Lep2fatjet.addHistogram("h_1Lep2fatjet_fatjet_msoftdrop_1", 180, 0, 250, [&]() { return nt.FatJet_msoftdrop()[ana.tx.getBranch<vector<int>>("Common_fatjet_idxs")[1]]; } );
    // hists_1Lep2fatjet.addHistogram("h_1Lep2fatjet_fatjet_genjetmass_0", 180, 0, 250, [&]() { return nt.GenJetAK8_mass()[nt.FatJet_genJetAK8Idx()[ana.tx.getBranch<vector<int>>("Common_fatjet_idxs")[0]]]; } );
    // hists_1Lep2fatjet.addHistogram("h_1Lep2fatjet_fatjet_genjetmass_1", 180, 0, 250, [&]() { return nt.GenJetAK8_mass()[nt.FatJet_genJetAK8Idx()[ana.tx.getBranch<vector<int>>("Common_fatjet_idxs")[1]]]; } );
    hists_1Lep2fatjet.addHistogram("h_1Lep2fatjet_fatjet_MDdeepW_0", 180, -1, 1, [&]() { return nt.FatJet_deepTagMD_WvsQCD()[ana.tx.getBranch<vector<int>>("Common_fatjet_idxs")[0]]; } );
    hists_1Lep2fatjet.addHistogram("h_1Lep2fatjet_fatjet_MDdeepW_1", 180, -1, 1, [&]() { return nt.FatJet_deepTagMD_WvsQCD()[ana.tx.getBranch<vector<int>>("Common_fatjet_idxs")[1]]; } );
    hists_1Lep2fatjet.addHistogram("h_1Lep2fatjet_fatjet_MDdeepW_sum", 180, -1, 2, [&]() { return nt.FatJet_deepTagMD_WvsQCD()[ana.tx.getBranch<vector<int>>("Common_fatjet_idxs")[0]] + nt.FatJet_deepTagMD_WvsQCD()[ana.tx.getBranch<vector<int>>("Common_fatjet_idxs")[1]]; } );
    hists_1Lep2fatjet.addHistogram("h_1Lep2fatjet_fatjet_MDdeepW_prod", 180, -1, 2, [&]() { return nt.FatJet_deepTagMD_WvsQCD()[ana.tx.getBranch<vector<int>>("Common_fatjet_idxs")[0]] * nt.FatJet_deepTagMD_WvsQCD()[ana.tx.getBranch<vector<int>>("Common_fatjet_idxs")[1]]; } );
    hists_1Lep2fatjet.addHistogram("h_1Lep2fatjet_fatjet_deepW_0", 180, -1, 1, [&]() { return nt.FatJet_deepTag_WvsQCD()[ana.tx.getBranch<vector<int>>("Common_fatjet_idxs")[0]]; } );
    hists_1Lep2fatjet.addHistogram("h_1Lep2fatjet_fatjet_deepW_1", 180, -1, 1, [&]() { return nt.FatJet_deepTag_WvsQCD()[ana.tx.getBranch<vector<int>>("Common_fatjet_idxs")[1]]; } );
    hists_1Lep2fatjet.addHistogram("h_1Lep2fatjet_fatjet_deepW_sum", 180, -1, 2, [&]() { return nt.FatJet_deepTag_WvsQCD()[ana.tx.getBranch<vector<int>>("Common_fatjet_idxs")[0]] + nt.FatJet_deepTag_WvsQCD()[ana.tx.getBranch<vector<int>>("Common_fatjet_idxs")[1]]; } );
    hists_1Lep2fatjet.addHistogram("h_1Lep2fatjet_fatjet_deepW_prod", 180, -1, 2, [&]() { return nt.FatJet_deepTag_WvsQCD()[ana.tx.getBranch<vector<int>>("Common_fatjet_idxs")[0]] * nt.FatJet_deepTag_WvsQCD()[ana.tx.getBranch<vector<int>>("Common_fatjet_idxs")[1]]; } );
    hists_1Lep2fatjet.addHistogram("h_1Lep2fatjet_fatjet_deepW_prod_zoom", 180, 0.8, 1, [&]() { return nt.FatJet_deepTag_WvsQCD()[ana.tx.getBranch<vector<int>>("Common_fatjet_idxs")[0]] * nt.FatJet_deepTag_WvsQCD()[ana.tx.getBranch<vector<int>>("Common_fatjet_idxs")[1]]; } );
    hists_1Lep2fatjet.addHistogram("h_1Lep2fatjet_fatjet_MDdeepT_0", 180, -1, 1, [&]() { return nt.FatJet_deepTagMD_TvsQCD()[ana.tx.getBranch<vector<int>>("Common_fatjet_idxs")[0]]; } );
    hists_1Lep2fatjet.addHistogram("h_1Lep2fatjet_fatjet_MDdeepT_1", 180, -1, 1, [&]() { return nt.FatJet_deepTagMD_TvsQCD()[ana.tx.getBranch<vector<int>>("Common_fatjet_idxs")[1]]; } );
    hists_1Lep2fatjet.addHistogram("h_1Lep2fatjet_fatjet_deepT_0", 180, -1, 1, [&]() { return nt.FatJet_deepTag_TvsQCD()[ana.tx.getBranch<vector<int>>("Common_fatjet_idxs")[0]]; } );
    hists_1Lep2fatjet.addHistogram("h_1Lep2fatjet_fatjet_deepT_1", 180, -1, 1, [&]() { return nt.FatJet_deepTag_TvsQCD()[ana.tx.getBranch<vector<int>>("Common_fatjet_idxs")[1]]; } );
    hists_1Lep2fatjet.addHistogram("h_1Lep2fatjet_fatjet_MDdeepbbvsL_0", 180, -1, 1, [&]() { return nt.FatJet_deepTagMD_bbvsLight()[ana.tx.getBranch<vector<int>>("Common_fatjet_idxs")[0]]; } );
    hists_1Lep2fatjet.addHistogram("h_1Lep2fatjet_fatjet_MDdeepbbvsL_1", 180, -1, 1, [&]() { return nt.FatJet_deepTagMD_bbvsLight()[ana.tx.getBranch<vector<int>>("Common_fatjet_idxs")[1]]; } );
    // hists_1Lep2fatjet.addHistogram("h_1Lep2fatjet_fatjet_lsf3_0", 180, -1, 1, [&]() { return nt.FatJet_lsf3()[ana.tx.getBranch<vector<int>>("Common_fatjet_idxs")[0]]; } );
    // hists_1Lep2fatjet.addHistogram("h_1Lep2fatjet_fatjet_lsf3_1", 180, -1, 1, [&]() { return nt.FatJet_lsf3()[ana.tx.getBranch<vector<int>>("Common_fatjet_idxs")[1]]; } );
    // hists_1Lep2fatjet.addHistogram("h_1Lep2fatjet_fatjet_nBHadrons_0", 5, 0, 5, [&]() { return nt.FatJet_nBHadrons()[ana.tx.getBranch<vector<int>>("Common_fatjet_idxs")[0]]; } );
    // hists_1Lep2fatjet.addHistogram("h_1Lep2fatjet_fatjet_nBHadrons_1", 5, 0, 5, [&]() { return nt.FatJet_nBHadrons()[ana.tx.getBranch<vector<int>>("Common_fatjet_idxs")[1]]; } );
    // hists_1Lep2fatjet.addHistogram("h_1Lep2fatjet_fatjet_nCHadrons_0", 5, 0, 5, [&]() { return nt.FatJet_nCHadrons()[ana.tx.getBranch<vector<int>>("Common_fatjet_idxs")[0]]; } );
    // hists_1Lep2fatjet.addHistogram("h_1Lep2fatjet_fatjet_nCHadrons_1", 5, 0, 5, [&]() { return nt.FatJet_nCHadrons()[ana.tx.getBranch<vector<int>>("Common_fatjet_idxs")[1]]; } );
    hists_1Lep2fatjet.addHistogram("h_1Lep2fatjet_fatjet_tau_21_0", 180, 0, 1, [&]() { return nt.FatJet_tau2()[ana.tx.getBranch<vector<int>>("Common_fatjet_idxs")[0]] / nt.FatJet_tau1()[ana.tx.getBranch<vector<int>>("Common_fatjet_idxs")[0]]; } );
    hists_1Lep2fatjet.addHistogram("h_1Lep2fatjet_fatjet_tau_21_1", 180, 0, 1, [&]() { return nt.FatJet_tau2()[ana.tx.getBranch<vector<int>>("Common_fatjet_idxs")[1]] / nt.FatJet_tau1()[ana.tx.getBranch<vector<int>>("Common_fatjet_idxs")[1]]; } );
    hists_1Lep2fatjet.addHistogram("h_1Lep2fatjet_fatjet_tau_21_sum", 180, 0, 2, [&]() { return nt.FatJet_tau2()[ana.tx.getBranch<vector<int>>("Common_fatjet_idxs")[0]] / nt.FatJet_tau1()[ana.tx.getBranch<vector<int>>("Common_fatjet_idxs")[0]] + nt.FatJet_tau2()[ana.tx.getBranch<vector<int>>("Common_fatjet_idxs")[1]] / nt.FatJet_tau1()[ana.tx.getBranch<vector<int>>("Common_fatjet_idxs")[1]]; } );
    hists_1Lep2fatjet.addHistogram("h_1Lep2fatjet_fatjet_tau_21_prod", 180, 0, 1, [&]() { return nt.FatJet_tau2()[ana.tx.getBranch<vector<int>>("Common_fatjet_idxs")[0]] / nt.FatJet_tau1()[ana.tx.getBranch<vector<int>>("Common_fatjet_idxs")[0]] * nt.FatJet_tau2()[ana.tx.getBranch<vector<int>>("Common_fatjet_idxs")[1]] / nt.FatJet_tau1()[ana.tx.getBranch<vector<int>>("Common_fatjet_idxs")[1]]; } );
    hists_1Lep2fatjet.addHistogram("h_1Lep2fatjet_fatjet_subjet_pt_0_0", 180, -20, 1000, [&]() { if (nt.FatJet_subJetIdx1()[ana.tx.getBranch<vector<int>>("Common_fatjet_idxs")[0]] >= 0) return nt.Jet_pt()[nt.FatJet_subJetIdx1()[ana.tx.getBranch<vector<int>>("Common_fatjet_idxs")[0]]]; else return -999.f; } );
    hists_1Lep2fatjet.addHistogram("h_1Lep2fatjet_fatjet_subjet_pt_1_0", 180, -20, 1000, [&]() { if (nt.FatJet_subJetIdx2()[ana.tx.getBranch<vector<int>>("Common_fatjet_idxs")[0]] >= 0) return nt.Jet_pt()[nt.FatJet_subJetIdx2()[ana.tx.getBranch<vector<int>>("Common_fatjet_idxs")[0]]]; else return -999.f; } );
    hists_1Lep2fatjet.addHistogram("h_1Lep2fatjet_fatjet_subjet_pt_0_1", 180, -20, 1000, [&]() { if (nt.FatJet_subJetIdx1()[ana.tx.getBranch<vector<int>>("Common_fatjet_idxs")[1]] >= 0) return nt.Jet_pt()[nt.FatJet_subJetIdx1()[ana.tx.getBranch<vector<int>>("Common_fatjet_idxs")[1]]]; else return -999.f; } );
    hists_1Lep2fatjet.addHistogram("h_1Lep2fatjet_fatjet_subjet_pt_1_1", 180, -20, 1000, [&]() { if (nt.FatJet_subJetIdx2()[ana.tx.getBranch<vector<int>>("Common_fatjet_idxs")[1]] >= 0) return nt.Jet_pt()[nt.FatJet_subJetIdx2()[ana.tx.getBranch<vector<int>>("Common_fatjet_idxs")[1]]]; else return -999.f; } );
    hists_1Lep2fatjet.addHistogram("h_1Lep2fatjet_fatjet_subjet_ptfrac_0_0", 180, 0, 1, [&]() { if (nt.FatJet_subJetIdx1()[ana.tx.getBranch<vector<int>>("Common_fatjet_idxs")[0]] >= 0) return nt.Jet_pt()[nt.FatJet_subJetIdx1()[ana.tx.getBranch<vector<int>>("Common_fatjet_idxs")[0]]] / ana.tx.getBranch<vector<LorentzVector>>("Common_fatjet_p4")[0].pt(); else return -999.f; } );
    hists_1Lep2fatjet.addHistogram("h_1Lep2fatjet_fatjet_subjet_ptfrac_1_0", 180, 0, 1, [&]() { if (nt.FatJet_subJetIdx2()[ana.tx.getBranch<vector<int>>("Common_fatjet_idxs")[0]] >= 0) return nt.Jet_pt()[nt.FatJet_subJetIdx2()[ana.tx.getBranch<vector<int>>("Common_fatjet_idxs")[0]]] / ana.tx.getBranch<vector<LorentzVector>>("Common_fatjet_p4")[0].pt(); else return -999.f; } );
    hists_1Lep2fatjet.addHistogram("h_1Lep2fatjet_fatjet_subjet_ptfrac_0_1", 180, 0, 1, [&]() { if (nt.FatJet_subJetIdx1()[ana.tx.getBranch<vector<int>>("Common_fatjet_idxs")[1]] >= 0) return nt.Jet_pt()[nt.FatJet_subJetIdx1()[ana.tx.getBranch<vector<int>>("Common_fatjet_idxs")[1]]] / ana.tx.getBranch<vector<LorentzVector>>("Common_fatjet_p4")[1].pt(); else return -999.f; } );
    hists_1Lep2fatjet.addHistogram("h_1Lep2fatjet_fatjet_subjet_ptfrac_1_1", 180, 0, 1, [&]() { if (nt.FatJet_subJetIdx2()[ana.tx.getBranch<vector<int>>("Common_fatjet_idxs")[1]] >= 0) return nt.Jet_pt()[nt.FatJet_subJetIdx2()[ana.tx.getBranch<vector<int>>("Common_fatjet_idxs")[1]]] / ana.tx.getBranch<vector<LorentzVector>>("Common_fatjet_p4")[1].pt(); else return -999.f; } );
    hists_1Lep2fatjet.addHistogram("h_1Lep2fatjet_fatjet_subjet_pt_ratio_0", 180, -20, 10, [&]() { if (nt.FatJet_subJetIdx1()[ana.tx.getBranch<vector<int>>("Common_fatjet_idxs")[0]] >= 0 and nt.FatJet_subJetIdx2()[ana.tx.getBranch<vector<int>>("Common_fatjet_idxs")[0]] >= 0) return nt.Jet_pt()[nt.FatJet_subJetIdx1()[ana.tx.getBranch<vector<int>>("Common_fatjet_idxs")[0]]] / nt.Jet_pt()[nt.FatJet_subJetIdx2()[ana.tx.getBranch<vector<int>>("Common_fatjet_idxs")[0]]]; else return -999.f; } );
    hists_1Lep2fatjet.addHistogram("h_1Lep2fatjet_fatjet_subjet_pt_ratio_1", 180, -20, 10, [&]() { if (nt.FatJet_subJetIdx1()[ana.tx.getBranch<vector<int>>("Common_fatjet_idxs")[1]] >= 0 and nt.FatJet_subJetIdx2()[ana.tx.getBranch<vector<int>>("Common_fatjet_idxs")[1]] >= 0) return nt.Jet_pt()[nt.FatJet_subJetIdx1()[ana.tx.getBranch<vector<int>>("Common_fatjet_idxs")[1]]] / nt.Jet_pt()[nt.FatJet_subJetIdx2()[ana.tx.getBranch<vector<int>>("Common_fatjet_idxs")[1]]]; else return -999.f; } );
/*    hists_1Lep2fatjet.addHistogram("h_1Lep2fatjet_htSum", 180, 0, 4000, [&]() { return ana.tx.getBranch<vector<LorentzVector>>("Common_fatjet_p4")[0].pt() + ana.tx.getBranch<vector<LorentzVector>>("Common_fatjet_p4")[1].pt() + ana.tx.getBranch<LorentzVector>("1Lep2fatjet_lep_p4").pt() + ana.tx.getBranch<LorentzVector>("Common_met_p4").pt(); } );
    hists_1Lep2fatjet.addHistogram("h_1Lep2fatjet_mass_lvJJ", 180, 0, 4000, [&]() { return (ana.tx.getBranch<vector<LorentzVector>>("Common_fatjet_p4")[0] + ana.tx.getBranch<vector<LorentzVector>>("Common_fatjet_p4")[1] + ana.tx.getBranch<LorentzVector>("1Lep2fatjet_lep_p4") + ana.tx.getBranch<LorentzVector>("Common_met_p4")).mass(); } );
    hists_1Lep2fatjet.addHistogram("h_1Lep2fatjet_mt_lvJJ", 180, 0, 4000, [&]() { return (ana.tx.getBranch<vector<LorentzVector>>("Common_fatjet_p4")[0] + ana.tx.getBranch<vector<LorentzVector>>("Common_fatjet_p4")[1] + ana.tx.getBranch<LorentzVector>("1Lep2fatjet_lep_p4") + ana.tx.getBranch<LorentzVector>("Common_met_p4")).mt(); } );
    hists_1Lep2fatjet.addHistogram("h_1Lep2fatjet_pt_lvJJ", 180, 0, 4000, [&]() { return (ana.tx.getBranch<vector<LorentzVector>>("Common_fatjet_p4")[0] + ana.tx.getBranch<vector<LorentzVector>>("Common_fatjet_p4")[1] + ana.tx.getBranch<LorentzVector>("1Lep2fatjet_lep_p4") + ana.tx.getBranch<LorentzVector>("Common_met_p4")).mass(); } );
    hists_1Lep2fatjet.addHistogram("h_1Lep2fatjet_mass_lvJ_0", 180, 0, 2000, [&]() { return (ana.tx.getBranch<vector<LorentzVector>>("Common_fatjet_p4")[0] + ana.tx.getBranch<LorentzVector>("1Lep2fatjet_lep_p4") + ana.tx.getBranch<LorentzVector>("Common_met_p4")).mass(); } );
    hists_1Lep2fatjet.addHistogram("h_1Lep2fatjet_mass_lvJ_1", 180, 0, 2000, [&]() { return (ana.tx.getBranch<vector<LorentzVector>>("Common_fatjet_p4")[1] + ana.tx.getBranch<LorentzVector>("1Lep2fatjet_lep_p4") + ana.tx.getBranch<LorentzVector>("Common_met_p4")).mass(); } );
    hists_1Lep2fatjet.addHistogram("h_1Lep2fatjet_mt_lvJ_0", 180, 0, 2000, [&]() { return (ana.tx.getBranch<vector<LorentzVector>>("Common_fatjet_p4")[0] + ana.tx.getBranch<LorentzVector>("1Lep2fatjet_lep_p4") + ana.tx.getBranch<LorentzVector>("Common_met_p4")).mt(); } );
    hists_1Lep2fatjet.addHistogram("h_1Lep2fatjet_mt_lvJ_1", 180, 0, 2000, [&]() { return (ana.tx.getBranch<vector<LorentzVector>>("Common_fatjet_p4")[1] + ana.tx.getBranch<LorentzVector>("1Lep2fatjet_lep_p4") + ana.tx.getBranch<LorentzVector>("Common_met_p4")).mt(); } );
*/

    // Now book cutflow histogram (could be commented out if user does not want.)
    // N.B. Cutflow histogramming can be CPU consuming.
    ana.cutflow.bookCutflows();

    // Book histograms to cuts that user wants for this category.
    ana.cutflow.bookHistogramsForCut(hists_1Lep2fatjet, "Cut_1Lep2fatjet_Has2FatJet");
//    ana.cutflow.bookHistogramsForCut(hists_1Lep2fatjet, "Cut_1Lep2fatjet_Has2FatJet_new");
//    ana.cutflow.bookHistogramsForCut(hists_1Lep2fatjet_new, "Cut_1Lep2fatjet_Has1Lepton_new");
    ana.cutflow.bookHistogramsForCut(hists_1Lep2fatjet_new, "Cut_1Lep2fatjet_Has2FatJet_new");
//    ana.cutflow.bookHistogramsForCut(hists_1Lep2fatjet, "Cut_1Lep2fatjet_Has1Lepton_new");
    ana.cutflow.bookHistogramsForCut(hists_1Lep2fatjet, "Cut_1Lep2fatjet_FatJetMassPresel");
    ana.cutflow.bookHistogramsForCut(hists_1Lep2fatjet, "Cut_1Lep2fatjet_Preselection");
    
    ana.tx.fill();
}
