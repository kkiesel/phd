#include "HistogramProducer.h"

void HistogramProducer::initTriggerStudies() {
  effMap["eff_hlt_pt"] = TEfficiency("", ";#it{p}_{T} (GeV);#varepsilon", 250, 0, 1000);
  effMap["eff_hlt_eta"] = TEfficiency("", ";|#eta|;#varepsilon", 15, 0, 1.5);
  effMap["eff_hlt_ht"] = TEfficiency("", ";#it{EMH}_{T} (GeV);#varepsilon", 200, 0, 2000);
  effMap["eff_hlt_ht_ptMin"] = TEfficiency("", ";#it{EMH}_{T} (GeV);#varepsilon", 200, 0, 2000);
  effMap["eff_hlt_ht_etaMax"] = TEfficiency("", ";#it{EMH}_{T} (GeV);#varepsilon", 310, 0, 3.1);
  effMap["eff_hlt_ht_ct"] = TEfficiency("", ";#it{EMH}_{T} (GeV);#varepsilon", 200, 0, 2000);
  effMap["eff_hlt_ht_ct_preScaled"] = TEfficiency("", ";#it{EMH}_{T} (GeV);#varepsilon (prescaled)", 200, 0, 2000);
  effMap["eff_hlt_ht_ct2"] = TEfficiency("", ";#it{EMH}_{T} (GeV);#varepsilon", 200, 0, 2000);
  effMap["eff_hlt_ht_ct2_preScaled"] = TEfficiency("", ";#it{EMH}_{T} (GeV);#varepsilon (prescaled)", 200, 0, 2000);
  effMap["eff_hlt_nVertex"] = TEfficiency("", ";vertex multiplicity", 41, -0.5, 40.5);
  effMap["eff_hlt_sie"] = TEfficiency("", ";#sigma_{i#etai#eta}", 400, 0, 0.02);
  effMap["eff_hlt_hoe"] = TEfficiency("", ";H/E", 100, 0, 0.15);
  effMap["eff_hlt_r9"] = TEfficiency("", ";r9", 110, 0, 1.1);
  effMap["eff_hlt_cIso"] = TEfficiency("", ";I_{#pi} (GeV)", 100, 0, 10);
  effMap["eff_hlt_nIso"] = TEfficiency("", ";I_{n} (GeV)", 100, 0, 20);
  effMap["eff_hlt_pIso"] = TEfficiency("", ";I_{#gamma} (GeV)", 100, 0, 20);
  effMap["eff_hlt_nJet"] = TEfficiency("", ";uncleaned jet multiplicity", 15, -0.5, 14.5);
  effMap["eff_hlt_met"] = TEfficiency("", ";#it{E}_{T}^{miss} (GeV)", 150, 0, 150);
  effMap["eff_hlt_met_ct"] = TEfficiency("", ";#it{E}_{T}^{miss} (GeV)", 150, 0, 150);
  effMap["eff_hlt_met_ct2"] = TEfficiency("", ";#it{E}_{T}^{miss} (GeV)", 150, 0, 150);
  effMap["eff_hlt_ele_pt"] = TEfficiency("", ";#it{p}_{T} (GeV);#varepsilon", 100, 0, 100);
  effMap["eff_hlt_pt_endcap"] = TEfficiency("", ";#it{p}_{T} (GeV);#varepsilon", 250, 0, 1000);
  effMap["eff_hlt_pt_endcap_ps"] = TEfficiency("", ";#it{p}_{T} (GeV);#varepsilon", 250, 0, 1000);
  effMap["eff_hlt_eta_endcap"] = TEfficiency("", ";|#eta|;#varepsilon", 15, 1.5, 3);
  effMap["eff_hlt_eta_endcap_ps"] = TEfficiency("", ";|#eta|;#varepsilon", 15, 1.5, 3);

}

void HistogramProducer::fillTriggerStudies() {

  tree::Photon* thisPhoton = 0;
  for (auto& photon : *photons) {
    if (photon.p.Pt() > 15 && fabs(photon.p.Eta()) < photonsEtaMaxBarrel && photon.isLoose && !photon.hasPixelSeed) {
      thisPhoton = &photon;
      break; // take leading(first) photon
    }
  }

  float ht = 0;
  if (thisPhoton) ht += thisPhoton->p.Pt();
  tree::Jet *minPtJet = 0;
  tree::Jet *maxEtaJet = 0;
  for (auto& jet : *jets) {
    if (jet.p.Pt() > 40 && fabs(jet.p.Eta()) < 3) {
      if (!minPtJet || minPtJet->p.Pt() > jet.p.Pt()) minPtJet = &jet;
      if (!maxEtaJet || fabs(minPtJet->p.Eta()) < fabs(jet.p.Eta())) maxEtaJet = &jet;
      if (!thisPhoton || jet.p.DeltaR(thisPhoton->p) > 0.3) {
        ht += jet.p.Pt();
      }
    }
  }

  if (thisPhoton) {

    // get trigger efficiencies for each run
    if (thisPhoton->p.Pt() > 100  && ht > 700 && *hlt_ht600) {
      if (!rawEff_vs_run.count(*runNo)) rawEff_vs_run[*runNo] = make_pair(0,0);
      if (*hlt_photon90_ht600) rawEff_vs_run.at(*runNo).first += 1;
      rawEff_vs_run.at(*runNo).second += 1;
    }

    // main trigger plots
    if (thisPhoton->p.Pt() > 100 && *hlt_photon90) {
      effMap.at("eff_hlt_ht").Fill(*hlt_photon90_ht600, ht);
      effMap.at("eff_hlt_ht_ptMin").Fill(*hlt_photon90_ht600, minPtJet->p.Pt());
      if (ht>600) effMap.at("eff_hlt_ht_etaMax").Fill(*hlt_photon90_ht600, fabs(maxEtaJet->p.Eta()));
    }
    if (thisPhoton->p.Pt() > 100 && *hlt_photon90_ht600) {
      effMap.at("eff_hlt_ht_ct").Fill(*hlt_ht600, ht);
      effMap.at("eff_hlt_ht_ct_preScaled").FillWeighted(*hlt_ht600, *hlt_ht600_pre, ht);
      if (ht > 700) effMap.at("eff_hlt_met_ct").Fill(*hlt_ht600, met->p.Pt());
    }
    if (thisPhoton->p.Pt() > 100 && *hlt_photon90) {
      effMap.at("eff_hlt_ht_ct2").Fill(*hlt_ht600, ht);
      effMap.at("eff_hlt_ht_ct2_preScaled").FillWeighted(*hlt_ht600, *hlt_ht600_pre, ht);
      if (ht > 700) effMap.at("eff_hlt_met_ct2").Fill(*hlt_ht600, met->p.Pt());
    }
    if (ht > 700 && *hlt_ht600) {
      effMap.at("eff_hlt_pt").Fill(*hlt_photon90_ht600, thisPhoton->p.Pt());
      if (thisPhoton->p.Pt()>100) {
        effMap.at("eff_hlt_eta").Fill(*hlt_photon90_ht600, fabs(thisPhoton->p.Eta()));
        effMap.at("eff_hlt_eta").Fill(*hlt_photon90_ht600, fabs(thisPhoton->p.Eta()));
      }
    }
  }

  if (ht > 700 && *hlt_ht600) {
    for (auto& el : *electrons) {
      if (fabs(el.p.Eta())>2.1 || !el.isTight) continue;
      effMap.at("eff_hlt_ele_pt").Fill(*hlt_el27, el.p.Pt());
      break; // only leading electron
    }
    for (auto& photon : *photons) {
      auto eta = fabs(photon.p.Eta());
      if (photon.p.Pt() > 100 && eta < photonsEtaMaxBarrel) {
        if (looseCutFlowPhoton.check(photon)) {
          effMap.at("eff_hlt_r9").Fill(*hlt_photon90_ht600, photon.r9);
          effMap.at("eff_hlt_nVertex").Fill(*hlt_photon90_ht600, *nGoodVertices);
          effMap.at("eff_hlt_nJet").Fill(*hlt_photon90_ht600, jets->size());
          effMap.at("eff_hlt_met").Fill(*hlt_photon90_ht600, met->p.Pt());
        }
        if (
          looseCutFlowPhoton.passHoe()
          && looseCutFlowPhoton.passSie()
          && looseCutFlowPhoton.passNIso()
          && looseCutFlowPhoton.passPIso()
         ) effMap.at("eff_hlt_cIso").Fill(*hlt_photon90_ht600, photon.cIso);
        if (
          looseCutFlowPhoton.passHoe()
          && looseCutFlowPhoton.passSie()
          && looseCutFlowPhoton.passCIso()
          && looseCutFlowPhoton.passPIso()
         ) effMap.at("eff_hlt_nIso").Fill(*hlt_photon90_ht600, photon.nIso);
        if (
          looseCutFlowPhoton.passHoe()
          && looseCutFlowPhoton.passSie()
          && looseCutFlowPhoton.passCIso()
          && looseCutFlowPhoton.passNIso()
         ) effMap.at("eff_hlt_pIso").Fill(*hlt_photon90_ht600, photon.pIso);
        if (
          looseCutFlowPhoton.passHoe()
          && looseCutFlowPhoton.passIso()
       ) effMap.at("eff_hlt_sie").Fill(*hlt_photon90_ht600, photon.sigmaIetaIeta);
        if (
          looseCutFlowPhoton.passSie()
          && looseCutFlowPhoton.passIso()
       ) effMap.at("eff_hlt_hoe").Fill(*hlt_photon90_ht600, photon.hOverE);
      }
      if (photon.isLoose && photonsEtaMinEndcap<eta && eta<photonsEtaMaxEndcap) {
        effMap.at("eff_hlt_pt_endcap").Fill(*hlt_photon90_ht600, photon.p.Pt());
        if (!photon.hasPixelSeed) effMap.at("eff_hlt_pt_endcap_ps").Fill(*hlt_photon90_ht600, photon.p.Pt());
        if (photon.p.Pt()>100) effMap.at("eff_hlt_eta_endcap").Fill(*hlt_photon90_ht600, fabs(photon.p.Eta()));
        if (photon.p.Pt()>100 && !photon.hasPixelSeed) effMap.at("eff_hlt_eta_endcap_ps").Fill(*hlt_photon90_ht600, fabs(photon.p.Eta()));
      }
    } // photon loop
  } // ht cross trigger
}

void HistogramProducer::initUncut() {
  map<string,TH1F> h;
  h["genHt"] = TH1F("", ";#it{H}_{T}^{gen}", 3000, 0, 3000);
  h["ht600_prescale"] = TH1F("", ";Prescales for HLT_PFHT600", 65, 0.5, 65.5);
  h1Maps["uncut"] = h;
  map<string,TH2F> h2;
  h2["dr_vs_relpt"] = TH2F("", ";#DeltaR;#it{p}_{T}^{jet}/#it{p}_{T}^{#gamma}", 100, 0, 0.5, 300, 0, 3);
  h2["dr_vs_relpt_genG"] = TH2F("", ";#DeltaR;#it{p}_{T}^{jet}/#it{p}_{T}^{#gamma}", 100, 0, 0.5, 300, 0, 3);
  h2Maps["uncut"] = h2;
}

void HistogramProducer::fillUncut() {
  h1Maps["uncut"]["genHt"].Fill(*genHt);
  h1Maps["uncut"]["ht600_prescale"].Fill(*hlt_ht600_pre);

  for (auto& g : *photons) {
    if (g.isLoose && !g.hasPixelSeed && g.p.Pt() > 100 && fabs(g.p.Eta()) < photonsEtaMaxBarrel) {

      // search gen match
      bool genMatch = false;
      for (auto& gen : *genParticles) {
        if (gen.pdgId == 22 && gen.isPrompt && gen.fromHardProcess && gen.p.DeltaR(g.p) < 0.3) {
          genMatch = true;
          break;
        }
      }

      for (auto& j : *jets) {
        auto dr = g.p.DeltaR(j.p);
        auto relPt = j.p.Pt()/g.p.Pt();
        h2Maps["uncut"]["dr_vs_relpt"].Fill(dr, relPt);
        if (genMatch) {
          h2Maps["uncut"]["dr_vs_relpt_genG"].Fill(dr, relPt);
        }
      }
    }
  }
}

map<string,TH2F> initHistograms2() {
  map<string,TH2F> hMap;

  hMap["n_jets_vs_photonPosition"] = TH2F("",";jet multiplicity;#gamma position",10, -0.5, 9.5, 10, -0.5, 9.5);
  hMap["g_eta_vs_g_phi"] = TH2F("",";|#eta|;|#phi|", 260, 2.6, 2.6, 100, -3.1, 3.1);
  hMap["met_vs_emht"] = TH2F("", ";#it{E}_{T}^{miss} (GeV);#it{EMH}_{T} (GeV)", 300, 0, 3000, 450, 500, 5000);
  hMap["metPar_vs_emht"] = TH2F("", ";#it{E}_{T}^{miss} #parallel (GeV);#it{EMH}_{T} (GeV)", 600, -3000, 3000, 450, 500, 5000);
  hMap["metPer_vs_emht"] = TH2F("", ";#it{E}_{T}^{miss #perp  } (GeV);#it{EMH}_{T} (GeV)", 300, 0, 3000, 450, 500, 5000);
  hMap["metRaw_vs_emht"] = TH2F("", ";uncorrected #it{E}_{T}^{miss} (GeV);#it{EMH}_{T} (GeV)", 300, 0, 3000, 450, 500, 5000);
  hMap["metParRaw_vs_emht"] = TH2F("", ";uncorrected #it{E}_{T}^{miss} #parallel (GeV);#it{EMH}_{T} (GeV)", 600, -3000, 3000, 450, 500, 5000);
  hMap["metPerRaw_vs_emht"] = TH2F("", ";uncorrected #it{E}_{T}^{miss #perp  } (GeV);#it{EMH}_{T} (GeV)", 300, 0, 3000, 450, 500, 5000);
  hMap["met_vs_n_jet"] = TH2F("", ";#it{E}_{T}^{miss} (GeV);N_{jet}", 300, 0, 3000, 15, -.5, 14.5);
  hMap["met_vs_n_obj"] = TH2F("", ";#it{E}_{T}^{miss} (GeV);N_{jet}", 300, 0, 3000, 15, -.5, 14.5);
  hMap["metRaw_vs_n_jet"] = TH2F("", ";uncorrected #it{E}_{T}^{miss} (GeV);N_{jet}", 300, 0, 3000, 15, -.5, 14.5);
  hMap["metRaw_vs_n_obj"] = TH2F("", ";uncorrected #it{E}_{T}^{miss} (GeV);N_{jet}", 300, 0, 3000, 15, -.5, 14.5);
  hMap["met_vs_g_pt"] = TH2F("", ";#it{E}_{T}^{miss} (GeV);#it{p}_{T} (GeV)", 300, 0, 3000, 100, 0, 1000);
  hMap["metRaw_vs_g_pt"] = TH2F("", ";#it{E}_{T}^{miss} (GeV);#it{p}_{T} (GeV)", 300, 0, 3000, 100, 0, 1000);
  hMap["met_vs_g_e"] = TH2F("", ";#it{E}_{T}^{miss} (GeV);#it{E} (GeV)", 300, 0, 3000, 100, 0, 1000);
  hMap["met_vs_j1_pt"] = TH2F("", ";#it{E}_{T}^{miss} (GeV);#it{p}_{T}^{jet} (GeV)", 300, 0, 3000, 200, 0, 2000);
  hMap["met_vs_mht"] = TH2F("", ";#it{E}_{T}^{miss} (GeV);|#vec{H}_{T}| (GeV)", 300, 0, 3000, 1500, 0, 1500);
  hMap["memht_vs_emht"] = TH2F("", ";#it{EMH}_{T}^{miss} (GeV);#it{EMH}_{T} (GeV)", 300, 0, 3000, 450, 500, 5000);
  hMap["mht_vs_emht"] = TH2F("", ";#it{H}_{T}^{miss} (GeV);#it{EMH}_{T} (GeV)", 300, 0, 3000, 450, 500, 5000);
  hMap["razor"] = TH2F("",";M_{R} (GeV);R^{2}", 500, 0, 5000, 200, 0, 2);
  hMap["st_vs_n_jet"] = TH2F("",";S_{T} (GeV);jet multiplicity", 300, 0, 3000, 10, -.5, 9.5);

  return hMap;
}

map<string,TH1F> initHistograms() {
  map<string,TH1F> hMap;

  // mets
  hMap["met"] = TH1F("", ";#it{E}_{T}^{miss} (GeV)", 200, 0, 2000);
  hMap["metPhotonPtUp"] = TH1F("", ";#it{E}_{T}^{miss} (GeV)", 200, 0, 2000);
  hMap["metPhotonPtDn"] = TH1F("", ";#it{E}_{T}^{miss} (GeV)", 200, 0, 2000);
  hMap["metSmearedPhotonByJERUp"] = TH1F("", ";smeared #it{E}_{T}^{miss} (GeV)", 200, 0, 2000);
  hMap["metSmearedPhotonByJERDn"] = TH1F("", ";smeared #it{E}_{T}^{miss} (GeV)", 200, 0, 2000);
  hMap["metSmearedPhotonByJER"] = TH1F("", ";smeared #it{E}_{T}^{miss} (GeV)", 200, 0, 2000);
  hMap["metSmearedPhotonByJERUp_fixedRho"] = TH1F("", ";smeared #it{E}_{T}^{miss} (GeV)", 200, 0, 2000);
  hMap["metSmearedPhotonByJERDn_fixedRho"] = TH1F("", ";smeared #it{E}_{T}^{miss} (GeV)", 200, 0, 2000);
  hMap["metSmearedPhotonByJER_fixedRho"] = TH1F("", ";smeared #it{E}_{T}^{miss} (GeV)", 200, 0, 2000);
  hMap["metPhotonResUp"] = TH1F("", ";smeared #it{E}_{T}^{miss} (GeV)", 200, 0, 2000);
  hMap["metPhotonResDn"] = TH1F("", ";smeared #it{E}_{T}^{miss} (GeV)", 200, 0, 2000);
  hMap["metUncertUp"] = TH1F("", ";#it{E}_{T}^{miss} (GeV)", 200, 0, 2000);
  hMap["metUncertDn"] = TH1F("", ";#it{E}_{T}^{miss} (GeV)", 200, 0, 2000);
  hMap["metJet1ResPhi"] = TH1F("", ";#it{E}_{T}^{miss} (GeV)", 200, 0, 2000);
  hMap["metJet1Res"] = TH1F("", ";#it{E}_{T}^{miss} (GeV)", 200, 0, 2000);
  hMap["metJet2Res"] = TH1F("", ";#it{E}_{T}^{miss} (GeV)", 200, 0, 2000);
  hMap["metJet3Res"] = TH1F("", ";#it{E}_{T}^{miss} (GeV)", 200, 0, 2000);
  hMap["metJet4Res"] = TH1F("", ";#it{E}_{T}^{miss} (GeV)", 200, 0, 2000);
  hMap["metJet5Res"] = TH1F("", ";#it{E}_{T}^{miss} (GeV)", 200, 0, 2000);


//  int maxNjets = 8;
//  for (int i=0; i<maxNjets; i++) {
//    for (int j=0; j<i; j++) {
//      string hname = "metJetSmeared_"+to_string(j)+"_"+to_string(i);
//      hMap[hname+"_up"] = TH1F("", ";#it{E}_{T}^{miss} (GeV)", 200, 0, 2000);
//      hMap[hname+"_dn"] = TH1F("", ";#it{E}_{T}^{miss} (GeV)", 200, 0, 2000);
//    }
//  }

  // met projections
  hMap["metPar"] = TH1F("", ";#it{E}_{T}^{miss} #parallel (GeV)", 400, -2000, 2000);
  hMap["metParUp"] = TH1F("", ";#it{E}_{T}^{miss} #parallel (GeV)", 400, -2000, 2000);
  hMap["metParDn"] = TH1F("", ";#it{E}_{T}^{miss} #parallel (GeV)", 400, -2000, 2000);
  hMap["metPer"] = TH1F("", ";#it{E}_{T}^{miss #perp  } (GeV)", 200, 0, 2000);
  hMap["metRaw"] = TH1F("", ";uncorrected #it{E}_{T}^{miss} (GeV)", 200, 0, 2000);
  hMap["metParRaw"] = TH1F("", ";uncorrected #it{E}_{T}^{miss} #parallel (GeV)", 400, -2000, 2000);
  hMap["metPerRaw"] = TH1F("", ";uncorrected #it{E}_{T}^{miss #perp  } (GeV)", 200, 0, 2000);
  hMap["mt_g_met"] = TH1F("", ";m(#gamma,#it{E}_{T}^{miss}_{T}) (GeV)", 150, 0, 1500);

  hMap["metSig"] = TH1F("", ";#it{S}", 3000, 0, 3000);
  hMap["tremht"] = TH1F("", ";#it{EMH}_{T}^{trigger-like} (GeV)", 300, 0, 3000);
  hMap["emht"] = TH1F("", ";#it{EMH}_{T} (GeV)", 300, 0, 3000);
  hMap["emhtNoLep"] = TH1F("", ";lepton cleaned #it{EMH}_{T} (GeV)", 300, 0, 3000);
  hMap["emhtStar"] = TH1F("", ";#it{EMH}_{T}* (GeV)", 300, 0, 3000);
  hMap["ht"] = TH1F("", ";#it{H}_{T} (GeV)", 300, 0, 3000);
  hMap["st"] = TH1F("", ";#it{S}_{T} (GeV)", 200, 0, 4000);
  hMap["emrecoilt"] = TH1F("", ";#vec{#it{EMH}}_{T} (GeV)", 150, 0, 1500);
  hMap["recoilt"] = TH1F("", ";#vec{#it{H}}_{T} (GeV)", 150, 0, 1500);

  // photon
  hMap["g_pt"] = TH1F("", ";#it{p}_{T} (GeV)", 150, 0, 1500);
  hMap["g_e"] = TH1F("", ";#it{E} (GeV)", 150, 0, 1500);
  hMap["g_ptStar"] = TH1F("", ";#it{p}_{T}* (GeV)", 150, 0, 1500);
  hMap["g_eta"] = TH1F("", ";|#eta|", 2600, 0, 2.6);

  // he-jet
  hMap["j1_pt"] = TH1F("", ";#it{p}_{T}^{1.jet} (GeV)", 150, 0, 1500);
  hMap["j1_eta"] = TH1F("", ";|#eta^{1.jet}|", 150, 0, 3);
  hMap["j2_pt"] = TH1F("", ";#it{p}_{T}^{2.jet} (GeV)", 150, 0, 1500);
  hMap["j2_eta"] = TH1F("", ";|#eta^{2.jet}|", 150, 0, 3);
  hMap["j3_pt"] = TH1F("", ";#it{p}_{T}^{3.jet} (GeV)", 150, 0, 1500);
  hMap["j3_eta"] = TH1F("", ";|#eta^{3.jet}|", 150, 0, 3);

  // angles
  hMap["dphi_met_g"] = TH1F("", ";|#Delta#phi(#it{E}_{T}^{miss},#gamma)|", 70, 0, 3.5);
  hMap["dphi_met_j1"] = TH1F("", ";|#Delta#phi(#it{E}_{T}^{miss},1.jet)|", 70, 0, 3.5);
  hMap["dphi_met_j2"] = TH1F("", ";|#Delta#phi(#it{E}_{T}^{miss},2.jet)|", 70, 0, 3.5);
  hMap["dphi_met_j3"] = TH1F("", ";|#Delta#phi(#it{E}_{T}^{miss},3.jet)|", 70, 0, 3.5);
  hMap["dphi_met_recoil"] = TH1F("", ";|#Delta#phi(#it{E}_{T}^{miss},#vec{#it{H}}_{T})|", 70, 0, 3.5);
  hMap["dphi_g_j1"] = TH1F("", ";|#Delta#phi(#gamma,1.jet)|", 70, 0, 3.5);
  hMap["dphi_g_j2"] = TH1F("", ";|#Delta#phi(#gamma,2.jet)|", 70, 0, 3.5);

  // multiplicities
  hMap["n_vertex"] = TH1F("", ";vertex multiplicity", 61, -0.5, 60.5);
  hMap["n_vertex_unw"] = TH1F("", ";vertex multiplicity", 61, -0.5, 60.5);
  hMap["n_photon"] = TH1F("", ";photon multiplicity", 6, -0.5, 5.5);
  hMap["n_jet"] = TH1F("", ";jet multiplicity", 16, -0.5, 15.5);
  hMap["n_bjet"] = TH1F("", ";b-jet multiplicity", 16, -0.5, 15.5);
  hMap["n_hejet"] = TH1F("", ";(EB,#it{p}_{T}>100) jet multiplicity", 16, -0.5, 15.5);
  hMap["n_electron"] = TH1F("", ";electron multiplicity", 4, -0.5, 3.5);
  hMap["n_muon"] = TH1F("", ";muon multiplicity", 4, -0.5, 3.5);
  hMap["n_heJet"] = TH1F("", ";photon-like jet multiplicity", 11, -0.5, 10.5);
  hMap["n_tracksPV"] = TH1F("", ";PV track multiplicity", 51, -0.5, 50.5);

  hMap["genMatch"] = TH1F("", ";pdg id for gen match", 47, -23.5, 23.5);
  hMap["genHt"] = TH1F("", ";#it{H}_{T}^{gen}", 3000, 0, 3000);

  return hMap;
}

void HistogramProducer::fillSelection(string const& s, bool fillTree=false, float addWeight=1.) {
  auto weight = selW*addWeight;
  //auto signalPointName = getSignalPointName(*signal_nBinos, *signal_m1, *signal_m2);
  if (std::isnan(met->p.X()) and std::isnan(met->p.Y())) return;

  float tree_m, tree_mRaw, tree_w, tree_emht, tree_pt, tree_emrecoilt, tree_topWeight, tree_dPhi, tree_tremht, tree_ht;
  UInt_t tree_njet, tree_run;
  if (!h1Maps.count(s)) {
    h1Maps[s] = initHistograms();
    h2Maps[s] = initHistograms2();
    treeMap[s] = new TTree("simpleTree", "");
    treeMap[s]->Branch("met", &tree_m);
    treeMap[s]->Branch("emht", &tree_emht);
    treeMap[s]->Branch("tremht", &tree_tremht);
    treeMap[s]->Branch("weight", &tree_w);
 //   treeMap[s]->Branch("run", &tree_run);
    treeMap[s]->Branch("pt", &tree_pt);
    treeMap[s]->Branch("ht", &tree_ht);
//    treeMap[s]->Branch("dPhi", &tree_dPhi);
  }
  auto m1 = &h1Maps[s];
  auto m2 = &h2Maps[s];

  vector<float> dPhis;
  int nJets = 0;
  for (auto& j : *jets ) {
    if (nJets < 3 && j.p.Pt()>100 && fabs(j.p.Eta()<3)) {
      dPhis.push_back(fabs(met->p.DeltaPhi(j.p)));
      nJets ++;
    }
  }
  tree_dPhi = dPhis.size() ? *std::min_element(dPhis.begin(),dPhis.end()) : 10;
  tree_m = met->p.Pt();
  tree_mRaw = metRaw->p.Pt();
  tree_w = weight * sampleW;
  tree_pt = 0;
  tree_topWeight = topPtReweighting(*genParticles);
  tree_run = *runNo;
  tree_tremht = 0;
  for (auto& j : *jets) {
    auto pt = j.p.Pt();
    if (fabs(j.p.Eta())<3 and pt>30) tree_tremht += pt;
  }
  tree_ht = 0;
  for (auto& j: selJets) {
    auto pt = j->p.Pt();
    if (fabs(j->p.Eta())<3 and pt>30) tree_ht += pt;
  }

  if (selPhotons.size()) {
    tree_pt = selPhotons.at(0)->p.Pt();
  } else if (selJets.size()) {
    tree_pt = selJets.at(0)->p.Pt();
  }

  // calculate variables
  float ht = 0;
  float emhtNoLep = 0;
  TVector3 recoil(0,0,0);
  for (auto& j : selJets) {
    ht += j->p.Pt();
    recoil += j->p;
    if (!j->hasPhotonMatch && !j->hasElectronMatch && !j->hasMuonMatch) {
      emhtNoLep += j->p.Pt();
    }
  }

  TVector3 emrecoil = recoil;
  float emht = ht;
  float emhtStar = ht;
  for (auto& p : selPhotons) {
    emrecoil += p->p;
    emht += p->p.Pt();
    emhtNoLep += p->p.Pt();
    auto mJet = matchedJet(*p);
    if (mJet) emhtStar += mJet->p.Pt();
    else emhtStar = p->p.Pt();
  }
  tree_emht = emht;
  tree_emrecoilt = emrecoil.Pt();
  tree_njet = selPhotons.size() ? selJets.size() : selJets.size() - 1;
  if (fillTree) treeMap[s]->Fill();
  float st = emht + met->p.Pt();

  float tremht = 0;
  for (auto& jet : *jets) {
    auto pt = jet.p.Pt();
    if ( pt > 30 && fabs(jet.p.Eta()) < 3) {
      tremht += pt;
    }
  }

  int nHEJets = selHEJets.size();
//  for (int i=0; i<nHEJets && nHEJets<8; i++) { // histograms only created for up to 7 jets
//    string hname = "metJetSmeared_"+to_string(i)+"_"+to_string(nHEJets);
//    m1->at(hname+"_up").Fill((met->p+selHEJets.at(i)->p*selHEJets.at(i)->ptRes).Pt(), weight);
//    m1->at(hname+"_dn").Fill((met->p-selHEJets.at(i)->p*selHEJets.at(i)->ptRes).Pt(), weight);
//  }

  m1->at("metSig").Fill(met->sig, weight);
  m1->at("emhtNoLep").Fill(emhtNoLep, weight);
  m1->at("tremht").Fill(tremht, weight);
  m1->at("emht").Fill(emht, weight);
  m1->at("emhtStar").Fill(emhtStar, weight);
  m1->at("ht").Fill(ht, weight);
  m1->at("st").Fill(st, weight);
  m1->at("recoilt").Fill(recoil.Pt(), weight);
  m1->at("emrecoilt").Fill(emrecoil.Pt(), weight);

  m1->at("met").Fill(met->p.Pt(), weight);
  m1->at("metRaw").Fill(metRaw->p.Pt(), weight);
  m1->at("metUncertUp").Fill(met->p.Pt()+met->uncertainty, weight);
  m1->at("metUncertDn").Fill(met->p.Pt()-met->uncertainty, weight);

  if (selJets.size()>4) m1->at("metJet5Res").Fill((met->p+selJets.at(4)->p*rand.Gaus(0, selJets.at(4)->ptRes)).Pt(), weight);
  if (selJets.size()>3) m1->at("metJet4Res").Fill((met->p+selJets.at(3)->p*rand.Gaus(0, selJets.at(3)->ptRes)).Pt(), weight);
  if (selJets.size()>2) m1->at("metJet3Res").Fill((met->p+selJets.at(2)->p*rand.Gaus(0, selJets.at(2)->ptRes)).Pt(), weight);
  if (selJets.size()>1) m1->at("metJet2Res").Fill((met->p+selJets.at(1)->p*rand.Gaus(0, selJets.at(1)->ptRes)).Pt(), weight);
  if (selJets.size()>0) m1->at("metJet1Res").Fill((met->p+selJets.at(0)->p*rand.Gaus(0, selJets.at(0)->ptRes)).Pt(), weight);
  if (selJets.size()>0) {
    auto thisJet = selJets.at(0)->p;
    thisJet.SetPhi(thisJet.Phi()+rand.Gaus(0,selJets.at(0)->phiRes));
    m1->at("metJet1ResPhi").Fill((met->p+thisJet-selJets.at(0)->p).Pt(), weight);
  }

  if (selPhotons.size() > 0) {
    auto g = selPhotons.at(0);
    m1->at("genMatch").Fill(genMatchNegativePrompt(*g, *genParticles), weight);
    m1->at("metPhotonPtUp").Fill((met->p+g->sigmaPt/g->p.Pt()*g->p).Pt(), weight);
    m1->at("metPhotonPtDn").Fill((met->p-g->sigmaPt/g->p.Pt()*g->p).Pt(), weight);
    auto mJet = matchedJet(*g);
    if (mJet) {
      m1->at("g_ptStar").Fill(mJet->p.Pt(), weight);
      m1->at("metSmearedPhotonByJERUp").Fill((met->p+g->p*mJet->ptRes).Pt(), weight);
      m1->at("metSmearedPhotonByJERDn").Fill((met->p-g->p*mJet->ptRes).Pt(), weight);
      for (int i=0;i<1000;i++) m1->at("metSmearedPhotonByJER").Fill((met->p-g->p*rand.Gaus(0, mJet->ptRes)).Pt(), weight);
    }
    float res = resolution.get(g->p.Pt(), g->p.Eta(), 15.4);
    m1->at("metSmearedPhotonByJERUp_fixedRho").Fill((met->p+res*g->p).Pt(), weight);
    m1->at("metSmearedPhotonByJERDn_fixedRho").Fill((met->p-res*g->p).Pt(), weight);
    for (int i=0;i<1000;i++) m1->at("metSmearedPhotonByJER_fixedRho").Fill((met->p-rand.Gaus(0,res)*g->p).Pt(), weight);

    m1->at("mt_g_met").Fill((g->p + met->p).Pt(), weight);
    m1->at("g_pt").Fill(g->p.Pt(), weight);
    m1->at("g_e").Fill(g->p.Mag(), weight);
    m1->at("g_eta").Fill(fabs(g->p.Eta()), weight);
    float dphi_met_g = fabs(met->p.DeltaPhi(g->p));
    m1->at("dphi_met_g").Fill(dphi_met_g, weight);
    m1->at("metPar").Fill(met->p.Pt()*cos(dphi_met_g), weight);
    m1->at("metParUp").Fill(met->p.Pt()*cos(dphi_met_g)+g->sigmaPt, weight);
    m1->at("metParDn").Fill(met->p.Pt()*cos(dphi_met_g)-g->sigmaPt, weight);
    m1->at("metPer").Fill(fabs(met->p.Pt()*sin(dphi_met_g)), weight);
    m1->at("metParRaw").Fill(metRaw->p.Pt()*cos(met->p.DeltaPhi(g->p)), weight);
    m1->at("metPerRaw").Fill(fabs(metRaw->p.Pt()*sin(met->p.DeltaPhi(g->p))), weight);
    m2->at("metPar_vs_emht").Fill(met->p.Pt()*cos(dphi_met_g), emht, weight);
    m2->at("metPer_vs_emht").Fill(met->p.Pt()*sin(dphi_met_g), emht, weight);
    m2->at("metParRaw_vs_emht").Fill(metRaw->p.Pt()*cos(dphi_met_g), emht, weight);
    m2->at("metPerRaw_vs_emht").Fill(metRaw->p.Pt()*sin(dphi_met_g), emht, weight);
    unsigned photonPosition=0;
    for (;photonPosition<selJets.size() && selJets.at(photonPosition)->p.Pt() > g->p.Pt();photonPosition++);
    m2->at("n_jets_vs_photonPosition").Fill(selJets.size(), photonPosition, weight);
    m2->at("g_eta_vs_g_phi").Fill(g->p.Eta(), g->p.Phi(), weight);
  }

  if (selJets.size() > 2) {
    m1->at("j3_pt").Fill(selJets.at(2)->p.Pt(), weight);
    m1->at("j3_eta").Fill(fabs(selJets.at(2)->p.Eta()), weight);
    m1->at("dphi_met_j3").Fill(fabs(met->p.DeltaPhi(selJets.at(2)->p)), weight);
  }
  if (selJets.size() > 1) {
    m1->at("j2_pt").Fill(selJets.at(1)->p.Pt(), weight);
    m1->at("j2_eta").Fill(fabs(selJets.at(1)->p.Eta()), weight);
    m1->at("dphi_met_j2").Fill(fabs(met->p.DeltaPhi(selJets.at(1)->p)), weight);
    if (selPhotons.size() > 0) m1->at("dphi_g_j2").Fill(fabs(selPhotons.at(0)->p.DeltaPhi(selJets.at(1)->p)), weight);
  }
  if (selJets.size() > 0) {
    m1->at("j1_pt").Fill(selJets.at(0)->p.Pt(), weight);
    m1->at("j1_eta").Fill(fabs(selJets.at(0)->p.Eta()), weight);
    m1->at("dphi_met_j1").Fill(fabs(met->p.DeltaPhi(selJets.at(0)->p)), weight);
    if (selPhotons.size() > 0) m1->at("dphi_g_j1").Fill(fabs(selPhotons.at(0)->p.DeltaPhi(selJets.at(0)->p)), weight);
  }

  m1->at("dphi_met_recoil").Fill(fabs(met->p.DeltaPhi(recoil)), weight);

  m1->at("n_vertex").Fill(*nGoodVertices, weight);
  m1->at("n_vertex_unw").Fill(*nGoodVertices);
  m1->at("n_photon").Fill(selPhotons.size(), weight);
  m1->at("n_jet").Fill(selJets.size(), weight);
  m1->at("n_bjet").Fill(selBJets.size(), weight);
  m1->at("n_hejet").Fill(nHEJets, weight);
  m1->at("n_electron").Fill(selElectrons.size(), weight);
  m1->at("n_muon").Fill(selMuons.size(), weight);
  m1->at("n_heJet").Fill(selJets.size(), weight);
  m1->at("n_tracksPV").Fill(*nTracksPV, weight);
  m1->at("genHt").Fill(*genHt, weight);

  m2->at("met_vs_n_jet").Fill(met->p.Pt(), selJets.size(), weight);
  m2->at("met_vs_n_obj").Fill(met->p.Pt(), selJets.size()+selPhotons.size(), weight);
  m2->at("metRaw_vs_n_jet").Fill(metRaw->p.Pt(), selJets.size(), weight);
  m2->at("metRaw_vs_n_obj").Fill(metRaw->p.Pt(), selJets.size()+selPhotons.size(), weight);
  if (selJets.size()) m2->at("met_vs_j1_pt").Fill(met->p.Pt(), selJets.at(0)->p.Pt() , weight);
  if (selPhotons.size()) {
    m2->at("met_vs_g_pt").Fill(met->p.Pt(), selPhotons.at(0)->p.Pt() , weight);
    m2->at("metRaw_vs_g_pt").Fill(metRaw->p.Pt(), selPhotons.at(0)->p.Pt() , weight);
    m2->at("met_vs_g_e").Fill(met->p.Pt(), selPhotons.at(0)->p.Pt() , weight);
  } else {
    m2->at("metRaw_vs_g_pt").Fill(metRaw->p.Pt(), 0., weight);
  }
  m2->at("met_vs_emht").Fill(met->p.Pt(), emht, weight);
  m2->at("metRaw_vs_emht").Fill(metRaw->p.Pt(), emht, weight);
  m2->at("met_vs_mht").Fill(met->p.Pt(), recoil.Pt(), weight);
  m2->at("memht_vs_emht").Fill(emrecoil.Pt(), emht, weight);
  m2->at("mht_vs_emht").Fill(recoil.Pt(), emht, weight);

} // end fill


HistogramProducer::HistogramProducer():
  photons(fReader, "photons"),
  jets(fReader, "jets"),
  electrons(fReader, "electrons"),
  muons(fReader, "muons"),
  genJets(fReader, "genJets"),
  genParticles(fReader, "genParticles"),
  intermediateGenParticles(fReader, "intermediateGenParticles"),
  met(fReader, "met"),
  metRaw(fReader, "met_raw"),
  nGoodVertices(fReader, "nGoodVertices"),
  nTracksPV(fReader, "nTracksPV"),
  pu_weight(fReader, "pu_weight"),
  mc_weight(fReader, "mc_weight"),
  genHt(fReader, "genHt"),
  //rho(fReader, "rho"),
  runNo(fReader, "runNo"),
  hlt_photon90_ht600(fReader, "HLT_Photon90_CaloIdL_PFHT600_v"),
  hlt_photon90(fReader, "HLT_Photon90_v"),
  hlt_ht600(fReader, "HLT_PFHT600_v"),
  hlt_ht800(fReader, "HLT_PFHT800_v"),
  hlt_el27(fReader, "HLT_Ele27_eta2p1_WPTight_Gsf_v"),
  hlt_ht600_pre(fReader, "HLT_PFHT600_v_pre"),
  //signal_nBinos(fReader, "signal_nBinos"),
  //signal_m1(fReader, "signal_m1"),
  //signal_m2(fReader, "signal_m2"),
  looseCutFlowPhoton({{"sigmaIetaIeta_eb",0.0102}, {"cIso_eb",3.32}, {"nIso1_eb",1.92}, {"nIso2_eb",0.014}, {"nIso3_eb",0.000019}, {"pIso1_eb",0.81}, {"pIso2_eb",0.0053},
    {"sigmaIetaIeta_ee",0.0274}, {"cIso_ee",1.97}, {"nIso1_ee",11.86}, {"nIso2_ee",0.0139}, {"nIso3_ee",0.000025}, {"pIso1_ee",0.83}, {"pIso2_ee",0.0034} }),
  startTime(time(NULL)),
  rand()
{
}

void HistogramProducer::Init(TTree *tree)
{
  fReader.SetTree(tree);
  inputName = fReader.GetTree()->GetCurrentFile()->GetName();
  isData = inputName.find("Run201") != string::npos;
  resolution = Resolution(isData? "Spring16_25nsV6_DATA_PtResolution_AK4PFchs.txt": "Spring16_25nsV6_MC_PtResolution_AK4PFchs.txt");
  weighters["fakeRate_eta"] = Weighter("../plotter/weights.root", isData?"fakeRate__data_40pt_eta":"fakeRate__sim_40pt_eta");
  weighters["fakeRate_pt"] = Weighter("../plotter/weights.root", isData?"fakeRate__data_x17_EB_40pt_pt":"fakeRate__sim_x17_EB_40pt_pt");
  weighters["sf_photon_id_loose"] = Weighter("../plotter/data/egammaEffi.txt_SF2DLoose.root", "EGamma_SF2D");
  weighters["sf_photon_pixel"] = Weighter("../plotter/data/EleVeto_SFs_80X.root", "Scaling Factors_HasPix_InclusiveR9");
  weighters.at("sf_photon_id_loose").fillOverflow2d();
  weighters.at("sf_photon_pixel").fillOverflow2d();

  float lumi = 36.53e3; // pb^{-1}
  cutFlow = *((TH1F*)fReader.GetTree()->GetCurrentFile()->Get("TreeWriter/hCutFlow"));
  fReader.GetEntries(true); // jumps to last file
  string lastInputName = fReader.GetTree()->GetCurrentFile()->GetName();
  if (inputName!=lastInputName) {
    // adds cut flow of last file. This makes only scence if there is for 1 or two files
    cutFlow.Add((TH1F*)fReader.GetTree()->GetCurrentFile()->Get("TreeWriter/hCutFlow"));
  }
  float nGen = cutFlow.GetBinContent(2);
  sampleW = isData ? 1. : lumi * sampleCrossSection(inputName) / nGen;

  genPt130 = inputName.find("0to130") != string::npos;
  genHt600 = inputName.find("HT-0to600") != string::npos;

  initTriggerStudies();
  initUncut();
}

void HistogramProducer::SlaveBegin(TTree *tree)
{
}

void HistogramProducer::defaultSelection()
{
  for (auto& mu : *muons) {
    if (mu.p.Pt() < 15) continue;
    if (indexOfMatchedParticle<tree::Photon*>(mu, selPhotons, .3) >= 0) continue;
    selMuons.push_back(&mu);
  }
  for (auto& el : *electrons) {
    if (!el.isLoose || el.p.Pt() < 15) continue;
    if (indexOfMatchedParticle<tree::Photon*>(el, selPhotons, .3) >= 0) continue;
    selElectrons.push_back(&el);
  }
  for (auto& jet : *jets) {
    if (!jet.isLoose
//      || jet.hasPhotonMatch || jet.hasElectronMatch || jet.hasMuonMatch
      || indexOfMatchedParticle<tree::Photon*>(jet, selPhotons, .3) >= 0
      || jet.p.Pt() < 40 || fabs(jet.p.Eta()) > 3) continue;
    selJets.push_back(&jet);
    if (jet.p.Pt() > 100 && fabs(jet.p.Eta()) < 1.4442) selHEJets.push_back(&jet);
    if (jet.bDiscriminator > bTaggingWorkingPoints.at(CSVv2M) && fabs(jet.p.Eta()) < 2.5)
      selBJets.push_back(&jet);
  }
}

void HistogramProducer::fillSelectedPhotons(const tree::Particle& p) {
  tree::Photon photon;
  photon.p = p.p;
  artificialPhotons.push_back(photon);
  selPhotons.push_back(&artificialPhotons.back());
}

tree::Jet* HistogramProducer::matchedJet(const tree::Particle& p) {
  for (auto& jet : *jets) {
    if (jet.p.DeltaR(p.p) < .2) {
      return &jet;
    }
  }
  return 0;
}

float HistogramProducer::getPhotonWeight(const tree::Photon& p) {
  float weight = 1.;
  if (!isData) {
    auto pt = p.p.Pt();
    auto eta = p.p.Eta();
    // sf for id and electron veto
    weight = weighters.at("sf_photon_id_loose").getWeight(eta, pt) * weighters.at("sf_photon_pixel").getWeight(fabs(eta), pt);
    // trigger efficiency
    weight *= fabs(eta)<photonsEtaMaxBarrel ? 0.977 : 0.953;
  }
  return weight;
}

Bool_t HistogramProducer::Process(Long64_t entry)
{
  resetSelection();
  fReader.SetLocalEntry(entry);
  if (genPt130) {
    float gPt = 0;
    for (auto const& genP : *genParticles) {
      float pt = genP.p.Pt();
      if (fabs(genP.pdgId)==22 && pt>gPt) {
        gPt = pt;
      }
    }
    if (gPt>130) return kTRUE;
  }
  if (genHt600 && *genHt>600) {
    return kTRUE;
  }

  fillUncut();

  if (isData) fillTriggerStudies();

  // selection: nTracksPV >= 2 for e->gamma fake-rate
  if (*nTracksPV<2) return kTRUE;

  // set weight
  selW = *mc_weight * *pu_weight;

  /////////////////////////////////////////////////////////////////////////////
  // signal sample
  /////////////////////////////////////////////////////////////////////////////

  for (auto& photon : *photons) {
    if (photon.isLoose && !photon.hasPixelSeed && photon.p.Pt() > 100 && fabs(photon.p.Eta()) < photonsEtaMaxBarrel) {
      selPhotons.push_back(&photon);
    }
  }
  defaultSelection();

  float myHt=0;
  for (auto& p : selPhotons) myHt += p->p.Pt();
  for (auto& p : selJets) myHt += p->p.Pt();

  if (selPhotons.size() && myHt > 700 && (*hlt_photon90_ht600 || !isData)) {
    auto addWeight = getPhotonWeight(*selPhotons.at(0));
    fillSelection("tr", true, addWeight);
    if (!selElectrons.size() && !selMuons.size()) fillSelection("tr_noLep", true, addWeight);
    if (selPhotons.at(0)->isTrue == MATCHED_FROM_GUDSCB) fillSelection("tr_true_GUDSCB", false, addWeight);
    else if (selPhotons.at(0)->isTrue == MATCHED_FROM_PI0) fillSelection("tr_true_pi0", false, addWeight);
    else if (selPhotons.at(0)->isTrue == MATCHED_FROM_OTHER_SOURCES) fillSelection("tr_true_other", false, addWeight);
    else fillSelection("tr_unmatched", false, addWeight);
    if (selPhotons.at(0)->isTight) fillSelection("tr_tight", false, addWeight);
    if (met->p.Pt() < 100) fillSelection("tr_0met100", false, addWeight);
    else                    fillSelection("tr_100met", false, addWeight);
    fillSelection(string("tr_genWZ")+to_string(genMatchWZDecay(*selPhotons.at(0), *intermediateGenParticles)), false, addWeight);
    fillSelection(string("tr_gen")+to_string(genMatchNegativePrompt(*selPhotons.at(0), *genParticles)), false, addWeight);
    if (fabs(genMatchNegativePrompt(*selPhotons.at(0), *genParticles)) == 11) {
      fillSelection("tr_genE", true, addWeight);
    } else {
      fillSelection("tr_noGenE", true, addWeight);
    }

  }

  if (!selPhotons.size() && myHt > 700 && (*hlt_ht600 || !isData)) {
    fillSelection("tr_jControl", true, *hlt_ht600_pre);
    if (!selElectrons.size() && !selMuons.size()) fillSelection("tr_jControl_noLep", true, *hlt_ht600_pre);
    for (auto& j : selJets) {
      if (j->nef>.9 && j->p.Pt()>100 && fabs(j->p.Eta())<1.4442) {
        fillSelection("tr_jControl_neutralEM9", false, *hlt_ht600_pre);
        break;
      }
    }
    for (auto& j : selJets) {
      if (j->chf>.9 && j->p.Pt()>100 && fabs(j->p.Eta())<1.4442) {
        fillSelection("tr_jControl_neutralCH", false, *hlt_ht600_pre);
        break;
      }
    }
  }


  resetSelection();
  /////////////////////////////////////////////////////////////////////////////
  // ee signal sample
  /////////////////////////////////////////////////////////////////////////////

  for (auto& photon : *photons) {
    auto eta = fabs(photon.p.Eta());
    if (photon.isLoose && !photon.hasPixelSeed && photon.p.Pt() > 100 && photonsEtaMinEndcap < eta && eta < photonsEtaMaxEndcap) {
      selPhotons.push_back(&photon);
    }
  }
  defaultSelection();

  myHt=0;
  for (auto& p : selPhotons) myHt += p->p.Pt();
  for (auto& p : selJets) myHt += p->p.Pt();

  if (selPhotons.size() && myHt > 700 && (*hlt_photon90_ht600 || !isData)) {
    auto addWeight = getPhotonWeight(*selPhotons.at(0));
    fillSelection("tr_ee", true, addWeight);
    if (fabs(genMatchNegativePrompt(*selPhotons.at(0), *genParticles)) == 11) {
      fillSelection("tr_ee_genE", true, addWeight);
    } else {
      fillSelection("tr_noGenE_ee", true, addWeight);
    }
  }

  if (!selPhotons.size() && myHt > 700 && (*hlt_ht600 || !isData)) {
    fillSelection("tr_jControl_ee", true, *hlt_ht600_pre);
    if (!selElectrons.size() && !selMuons.size()) fillSelection("tr_jControl_noLep_ee", true, *hlt_ht600_pre);
  }

  resetSelection();
  /////////////////////////////////////////////////////////////////////////////
  // electron sample
  /////////////////////////////////////////////////////////////////////////////

  for (auto& photon : *photons) {
    if (photon.isLoose && photon.hasPixelSeed && photon.p.Pt() > 100 && fabs(photon.p.Eta()) < photonsEtaMaxBarrel) {
      selPhotons.push_back(&photon);
    }
  }
  defaultSelection();

  myHt=0;
  for (auto& p : selPhotons) myHt += p->p.Pt();
  for (auto& p : selJets) myHt += p->p.Pt();
  if (selPhotons.size() && myHt > 700 && (*hlt_photon90_ht600 || !isData)) {
    auto addWeight = getPhotonWeight(*selPhotons.at(0));
    fillSelection("tr_eControl", true, addWeight);
    fillSelection("tr_eControl_vtxWeighted", true,
        isData ? (1.238+0.0935**nGoodVertices)/100 : (1.051+0.0318**nGoodVertices)/100*addWeight);
    fillSelection("tr_eControl_etaWeighted", true,
        weighters.at("fakeRate_eta").getWeight(fabs(selPhotons.at(0)->p.Eta())) * addWeight);
    fillSelection("tr_eControl_ptWeighted", true,
        weighters.at("fakeRate_pt").getWeight(selPhotons.at(0)->p.Pt()) * addWeight);
  }


  resetSelection();
  /////////////////////////////////////////////////////////////////////////////
  // ee electron sample
  /////////////////////////////////////////////////////////////////////////////

  for (auto& photon : *photons) {
    auto eta = fabs(photon.p.Eta());
    if (photon.isLoose && photon.hasPixelSeed && photon.p.Pt() > 100 && photonsEtaMinEndcap < eta && eta < photonsEtaMaxEndcap) {
      selPhotons.push_back(&photon);
    }
  }
  defaultSelection();

  myHt=0;
  for (auto& p : selPhotons) myHt += p->p.Pt();
  for (auto& p : selJets) myHt += p->p.Pt();
  if (selPhotons.size() && myHt > 700 && (*hlt_photon90_ht600 || !isData)) {
    auto addWeight = getPhotonWeight(*selPhotons.at(0));
    fillSelection("tr_eControl_ee", true, addWeight);
  }

  resetSelection();
  return kTRUE;
}

template<typename T>
void save2File(const map<string,map<string,T>>& hMaps, TFile& file)
{
  for (auto& hMapIt : hMaps) {
    if (!file.Get(hMapIt.first.c_str())) {
      file.mkdir(hMapIt.first.c_str());
    }
    file.cd(hMapIt.first.c_str());
    for (auto& h : hMapIt.second) {
      h.second.Write(h.first.c_str(), TObject::kWriteDelete);
    }
    file.cd();
  }
}

void HistogramProducer::Terminate()
{
  auto outputName = getOutputFilename(inputName);
  TFile file(outputName.c_str(), "RECREATE");
  save2File(h1Maps, file);
  save2File(h2Maps, file);
  for (auto& metsMapIt : treeMap) {
    file.cd(metsMapIt.first.c_str());
    metsMapIt.second->Write("simpleTree", TObject::kWriteDelete);
    file.cd();
  }

  file.mkdir("triggerStudies");
  file.cd("triggerStudies");
  for (auto& effIt : effMap) effIt.second.Write(effIt.first.c_str(), TObject::kWriteDelete);
  file.GetDirectory("triggerStudies")->WriteObject(&rawEff_vs_run, "rawEff_vs_run");
  file.cd();

  cutFlow.Write("hCutFlow");
  file.Close();
  cout << "Created " << outputName << " in " << (time(NULL) - startTime)/60 << " min" << endl;
}

void HistogramProducer::resetSelection() {
  selPhotons.clear();
  selJets.clear();
  selBJets.clear();
  selHEJets.clear();
  selElectrons.clear();
  selMuons.clear();
  artificialPhotons.clear();
}
