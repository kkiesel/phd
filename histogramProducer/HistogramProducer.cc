#include "HistogramProducer.h"

void HistogramProducer::initTriggerStudies() {
  effMap["eff_pt__p90ht600__ht600"] = TEfficiency("", ";#it{p}_{T} (GeV);#varepsilon", 250, 0, 1000);
  effMap["eff_pt_ee__p90ht600__ht600"] = TEfficiency("", ";#it{p}_{T} (GeV);#varepsilon", 250, 0, 1000);
  effMap["eff_eta__p90ht600__ht600"] = TEfficiency("", ";|#eta|;#varepsilon", 300, 0, 3);
  effMap["eff_eta_ee__p90ht600__ht600"] = TEfficiency("", ";|#eta|;#varepsilon", 300, 0, 3);
  effMap["eff_nVertex__p90ht600__ht600"] = TEfficiency("", ";vertex multiplicity", 41, -0.5, 40.5);
  effMap["eff_sie__p90ht600__ht600"] = TEfficiency("", ";#sigma_{i#etai#eta}", 400, 0, 0.02);
  effMap["eff_hoe__p90ht600__ht600"] = TEfficiency("", ";H/E", 100, 0, 0.15);
  effMap["eff_r9__p90ht600__ht600"] = TEfficiency("", ";r9", 110, 0, 1.1);
  effMap["eff_cIso__p90ht600__ht600"] = TEfficiency("", ";I_{#pi} (GeV)", 100, 0, 10);
  effMap["eff_nIso__p90ht600__ht600"] = TEfficiency("", ";I_{n} (GeV)", 100, 0, 20);
  effMap["eff_pIso__p90ht600__ht600"] = TEfficiency("", ";I_{#gamma} (GeV)", 100, 0, 20);
  effMap["eff_hasPixelSeeds__p90ht600__ht600"] = TEfficiency("", ";has pixel seeds", 2, 0, 2);
  effMap["eff_nJet__p90ht600__ht600"] = TEfficiency("", ";uncleaned jet multiplicity", 15, -0.5, 14.5);
  effMap["eff_met__p90ht600__ht600"] = TEfficiency("", ";#it{E}_{T}^{miss} (GeV)", 250, 0, 250);
  effMap["eff_pt__ele27__ht600"] = TEfficiency("", ";#it{p}_{T} (GeV);#varepsilon", 100, 0, 100);
  effMap["eff_emht__ht800__ht600"] = TEfficiency("", ";#it{EMH}_{T} (GeV);#varepsilon", 200, 0, 2000);
  effMap["eff_emht__p90ht600__p90"] = TEfficiency("", ";#it{EMH}_{T} (GeV);#varepsilon", 200, 0, 2000);
  effMap["eff_emht__ht800__p90"] = TEfficiency("", ";#it{EMH}_{T} (GeV);#varepsilon", 200, 0, 2000);
  effMap["eff_emht__ht600__p90"] = TEfficiency("", ";#it{EMH}_{T} (GeV);#varepsilon", 200, 0, 2000);
  effMap["eff_emht__ht600__p90_ps"] = TEfficiency("", ";#it{EMH}_{T} (GeV);#varepsilon (prescaled)", 200, 0, 2000);
  effMap["eff_met__ht600__p90"] = TEfficiency("", ";#it{E}_{T}^{miss} (GeV)", 250, 0, 250);
  effMap["eff_met__ht600__p90_ps"] = TEfficiency("", ";#it{E}_{T}^{miss} (GeV) (prescaled)", 250, 0, 250);
  effMap.at("eff_emht__ht600__p90_ps").SetUseWeightedEvents();
  effMap.at("eff_met__ht600__p90_ps").SetUseWeightedEvents();
}

void HistogramProducer::fillTriggerStudies() {
  // Barrel selection
  tree::Photon* selPhoton = 0;
  for (auto& photon : *photons) {
    if (photon.p.Pt() > 50 && fabs(photon.p.Eta()) < photonsEtaMaxBarrel && photon.isLoose && !photon.hasPixelSeed) {
      selPhoton = &photon;
      break; // take leading(first) photon
    }
  }
  float emht = 0;
  if (selPhoton) emht += selPhoton->p.Pt();
  for (auto& jet : *jets) {
    if (jet.p.Pt() > 30 && fabs(jet.p.Eta()) < 3) {
      if (!selPhoton || jet.p.DeltaR(selPhoton->p) > 0.3) {
        emht += jet.p.Pt();
      }
    }
  }
  if (selPhoton) {
    if (selPhoton->p.Pt() > 100  && emht > 700 && *hlt_ht600) {
      if (!rawEff_vs_run.count(*runNo)) rawEff_vs_run[*runNo] = make_pair(0,0);
      if (*hlt_photon90_ht600) rawEff_vs_run.at(*runNo).first += 1;
      rawEff_vs_run.at(*runNo).second += 1;

      effMap.at("eff_eta__p90ht600__ht600").Fill(*hlt_photon90_ht600, fabs(selPhoton->p.Eta()));
      effMap.at("eff_met__p90ht600__ht600").Fill(*hlt_photon90_ht600, met->p.Pt());
      effMap.at("eff_nVertex__p90ht600__ht600").Fill(*hlt_photon90_ht600, *nGoodVertices);
      effMap.at("eff_nJet__p90ht600__ht600").Fill(*hlt_photon90_ht600, *nGoodVertices);
      effMap.at("eff_r9__p90ht600__ht600").Fill(*hlt_photon90_ht600, selPhoton->r9);
      effMap.at("eff_hasPixelSeeds__p90ht600__ht600").Fill(*hlt_photon90_ht600, selPhoton->hasPixelSeed);
      effMap.at("eff_cIso__p90ht600__ht600").Fill(*hlt_photon90_ht600, selPhoton->cIso);
      effMap.at("eff_nIso__p90ht600__ht600").Fill(*hlt_photon90_ht600, selPhoton->nIso);
      effMap.at("eff_pIso__p90ht600__ht600").Fill(*hlt_photon90_ht600, selPhoton->pIso);
      effMap.at("eff_sie__p90ht600__ht600").Fill(*hlt_photon90_ht600, selPhoton->sigmaIetaIeta);
      effMap.at("eff_hoe__p90ht600__ht600").Fill(*hlt_photon90_ht600, selPhoton->hOverE);
    }
    if (emht>700 && *hlt_ht600) {
      effMap.at("eff_pt__p90ht600__ht600").Fill(*hlt_photon90_ht600, selPhoton->p.Pt());
    }

    if (selPhoton->p.Pt() > 100 && *hlt_photon90) {
      effMap.at("eff_emht__p90ht600__p90").Fill(*hlt_photon90_ht600, emht);
      effMap.at("eff_emht__ht800__p90").Fill(*hlt_ht800, emht);
      effMap.at("eff_emht__ht600__p90").Fill(*hlt_ht600, emht);
      effMap.at("eff_emht__ht600__p90_ps").FillWeighted(*hlt_ht600, *hlt_ht600_pre, emht);
      if (emht > 700) {
        effMap.at("eff_met__ht600__p90").Fill(*hlt_ht600, met->p.Pt());
        effMap.at("eff_met__ht600__p90_ps").FillWeighted(*hlt_ht600, *hlt_ht600_pre, met->p.Pt());
      }
    }
  }
  // Endcap selection
  selPhoton = 0;
  for (auto& photon : *photons) {
    auto eta = fabs(photon.p.Eta());
    if (photon.p.Pt() > 50 && photonsEtaMinEndcap < eta && eta < photonsEtaMaxEndcap)
      selPhoton = &photon;
      break; // take leading(first) photon
  }
  emht = 0;
  if (selPhoton) emht += selPhoton->p.Pt();
  for (auto& jet : *jets) {
    if (jet.p.Pt() > 30 && fabs(jet.p.Eta()) < 3) {
      if (!selPhoton || jet.p.DeltaR(selPhoton->p) > 0.3) {
        emht += jet.p.Pt();
      }
    }
  }
  if (selPhoton) {
     if (emht>700 && *hlt_ht600) {
      effMap.at("eff_pt_ee__p90ht600__ht600").Fill(*hlt_photon90_ht600, selPhoton->p.Pt());
      effMap.at("eff_eta_ee__p90ht600__ht600").Fill(*hlt_photon90_ht600, fabs(selPhoton->p.Eta()));
    }
  }

  // Selection without photon
  emht = 0;
  for (auto& jet : *jets) {
    if (jet.p.Pt() > 30 && fabs(jet.p.Eta()) < 3) {
        emht += jet.p.Pt();
      }
  }

  if (emht > 700 && *hlt_ht600) {
    effMap.at("eff_emht__ht800__ht600").Fill(*hlt_ht800, emht);
    for (auto& el : *electrons) {
      if (fabs(el.p.Eta())>2.1 || !el.isTight) continue;
      effMap.at("eff_pt__ele27__ht600").Fill(*hlt_el27, el.p.Pt());
      break; // only leading electron
    }
  }
}

void HistogramProducer::initUncut() {
  map<string,TH1F> h;
  h["genHt"] = TH1F("", ";#it{H}_{T}^{gen}", 3000, 0, 3000);
  h["ht600_prescale"] = TH1F("", ";Prescales for HLT_PFHT600", 65, 0.5, 65.5);
  h["johannesSt"] = TH1F("", ";S_{T}^{#gamma} (GeV)", 10, 600, 1600);
  h1Maps["uncut"] = h;

  map<string,TH2F> h2;
  h2["dr_vs_relpt"] = TH2F("", ";#DeltaR;#it{p}_{T}^{jet}/#it{p}_{T}^{#gamma}", 100, 0, 0.5, 300, 0, 3);
  h2["dr_vs_relpt_genG"] = TH2F("", ";#DeltaR;#it{p}_{T}^{jet}/#it{p}_{T}^{#gamma}", 100, 0, 0.5, 300, 0, 3);
  h2Maps["uncut"] = h2;

  effMap["noPromptEvaluation"] = TEfficiency("", ";;#varepsilon(isPrompt)", 1, 0, 1);
}

void HistogramProducer::fillUncut() {
  h1Maps["uncut"].at("genHt").Fill(*genHt);
  h1Maps["uncut"].at("ht600_prescale").Fill(*hlt_ht600_pre);

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
        h2Maps["uncut"].at("dr_vs_relpt").Fill(dr, relPt);
        if (genMatch) {
          h2Maps["uncut"].at("dr_vs_relpt_genG").Fill(dr, relPt);
        }
      }
    }
  }
  // Johannes selection:
  vector<TVector3> selPhotonsJo;
  for (auto& g : *photons) {
    if (g.isLoose && !g.hasPixelSeed && g.p.Pt() > 180 && fabs(g.p.Eta()) < photonsEtaMaxBarrel) {
      bool foundJet = false;
      for (auto& j: *jets) {
        if (j.isLoose && !j.hasPhotonMatch && j.p.DeltaR(g.p)<.5) {
          foundJet = true;
        }
      }
      if (!foundJet) {
        selPhotonsJo.push_back(g.p);
      }
    }
  }
  if (selPhotonsJo.size() && selPhotonsJo.at(0).Pt()>180 && met->p.Pt()>300 && transverseMass(selPhotonsJo.at(0), met->p)>300) {
    auto st = met->p.Pt();
    for (auto& g : selPhotonsJo) {
      st += g.Pt();
    }
    h1Maps["uncut"].at("johannesSt").Fill(st, *mc_weight * *pu_weight);
  }
}

map<string,TH2F> initHistograms2() {
  map<string,TH2F> hMap;

  hMap["n_jets_vs_photonPosition"] = TH2F("",";jet multiplicity;#gamma position",10, -0.5, 9.5, 10, -0.5, 9.5);
  hMap["g_eta_vs_g_phi"] = TH2F("",";|#eta|;|#phi|", 260, 2.6, 2.6, 100, -3.1, 3.1);
  hMap["met_vs_emht"] = TH2F("", ";#it{E}_{T}^{miss} (GeV);#it{EMH}_{T} (GeV)", 300, 0, 3000, 450, 500, 5000);
  hMap["met_vs_emht_JESu"] = TH2F("", ";#it{E}_{T}^{miss} (GeV);#it{EMH}_{T} (GeV)", 300, 0, 3000, 450, 500, 5000);
  hMap["met_vs_emht_JESd"] = TH2F("", ";#it{E}_{T}^{miss} (GeV);#it{EMH}_{T} (GeV)", 300, 0, 3000, 450, 500, 5000);
  hMap["met_vs_emht_JERu"] = TH2F("", ";#it{E}_{T}^{miss} (GeV);#it{EMH}_{T} (GeV)", 300, 0, 3000, 450, 500, 5000);
  hMap["met_vs_emht_JERd"] = TH2F("", ";#it{E}_{T}^{miss} (GeV);#it{EMH}_{T} (GeV)", 300, 0, 3000, 450, 500, 5000);
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
  hMap["metCrystalSeedCorrected"] = TH1F("", ";#it{E}_{T}^{miss} (GeV)", 200, 0, 2000);
  hMap["metCrystalSeedCorrected2"] = TH1F("", ";#it{E}_{T}^{miss} (GeV)", 200, 0, 2000);
  hMap["metCrystalSeedCorrected3"] = TH1F("", ";#it{E}_{T}^{miss} (GeV)", 200, 0, 2000);
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
  hMap["metTimesPt"] = TH1F("", ";#vec{#it{E}}_{T}^{miss}#upoint #vec{#it{p}}_{T} (GeV^{2})", 400, -200000, 200000);

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
  hMap["g_sie"] = TH1F("", ";#sigma_{i#etai#eta}", 400, 0, 0.04);
  hMap["g_sip"] = TH1F("", ";#sigma_{i#phii#phi}", 500, 0, 0.1);
  hMap["g_hoe"] = TH1F("", ";H/E", 500, 0, 0.5);
  hMap["g_r9"] = TH1F("", ";r9", 100, 0, 1);
  hMap["g_cIso"] = TH1F("", ";I_{#pm} (GeV)", 350, 0, 3.5);
  hMap["g_nIso"] = TH1F("", ";I_{n} (GeV)", 400, 0, 40);
  hMap["g_pIso"] = TH1F("", ";I_{p} (GeV)", 400, 0, 4);
  hMap["g_cIsoWorst"] = TH1F("", ";worst I_{#pm} (GeV)", 600, 0, 600);

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

void HistogramProducer::fillSelection(string const& s, float addWeight=1., bool fillTree=false) {
  auto weight = selW*addWeight;
  if (std::isnan(met->p.X()) and std::isnan(met->p.Y())) return;

  float tree_met, tree_metRaw, tree_emht, tree_w;
  float tree_dPhi, tree_pt, tree_eta, tree_eCrystal;
  bool tree_hasPromptPhoton = false;
  bool tree_genE = false;

  if (!h1Maps.count(s)) {
    h1Maps[s] = initHistograms();
    h2Maps[s] = initHistograms2();
    treeMap[s] = new TTree("simpleTree", "");
    treeMap[s]->Branch("met", &tree_met);
    treeMap[s]->Branch("metRaw", &tree_metRaw);
    treeMap[s]->Branch("emht", &tree_emht);
    treeMap[s]->Branch("weight", &tree_w);
    treeMap[s]->Branch("dPhi", &tree_dPhi);
    treeMap[s]->Branch("pt", &tree_pt);
    treeMap[s]->Branch("eta", &tree_eta);
    treeMap[s]->Branch("eCrystal", &tree_eCrystal);
    treeMap[s]->Branch("hasPromptPhoton", &tree_hasPromptPhoton, "hasPromptPhoton/O");
    treeMap[s]->Branch("genE", &tree_genE, "genE/O");
  }
  auto m1 = &h1Maps[s];
  auto m2 = &h2Maps[s];

  tree_met = met->p.Pt();
  tree_metRaw = metRaw->p.Pt();
  tree_w = weight * sampleW;
  if (selPhotons.size()) {
    auto g = selPhotons.at(0);
    tree_hasPromptPhoton = noPromptPhotons && count_if(genParticles->begin(), genParticles->end(), [] (const tree::GenParticle& p) { return p.pdgId==22 && p.promptStatus == DIRECTPROMPT;});
    tree_genE = fabs(genMatchNegativePrompt(*g, *genParticles)) == 11;
    tree_dPhi = fabs(met->p.DeltaPhi(g->p));
    tree_eta = g->p.Eta();
    tree_eCrystal = g->seedCrystalE;
    tree_pt = g->p.Pt();
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

    // correcting the barrel crystal seed energy
    float corr = isData && fabs(g->p.Eta()) < 1.4442 ? crystalResponse(g->seedCrystalE) : 1;
    auto corrMet1 = (met->p+g->p*(1-corr)).Pt();
    auto corrMet2 = (met->p-g->p*(1-corr)).Pt();
    auto corrMet3 = min(corrMet1, corrMet2);
    m1->at("metCrystalSeedCorrected").Fill((met->p+g->p*(1-corr)).Pt(), weight);
    m1->at("metCrystalSeedCorrected2").Fill((met->p-g->p*(1-corr)).Pt(), weight);
    m1->at("metCrystalSeedCorrected3").Fill(corrMet3, weight);

    m1->at("mt_g_met").Fill(transverseMass(g->p, met->p), weight);
    m1->at("metTimesPt").Fill(g->p * met->p, weight);
    m1->at("g_pt").Fill(g->p.Pt(), weight);
    m1->at("g_e").Fill(g->p.Mag(), weight);
    m1->at("g_eta").Fill(fabs(g->p.Eta()), weight);
    m1->at("g_sie").Fill(g->sigmaIetaIeta, weight);
    m1->at("g_sip").Fill(g->sigmaIphiIphi, weight);
    m1->at("g_hoe").Fill(g->hOverE, weight);
    m1->at("g_r9").Fill(g->r9, weight);
    m1->at("g_cIso").Fill(g->cIso, weight);
    m1->at("g_nIso").Fill(g->nIso, weight);
    m1->at("g_pIso").Fill(g->pIso, weight);
    m1->at("g_cIsoWorst").Fill(g->cIsoWorst, weight);

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
  m2->at("met_vs_emht_JESu").Fill(met_JESu->p.Pt(), emht, weight);
  m2->at("met_vs_emht_JESd").Fill(met_JESd->p.Pt(), emht, weight);
  m2->at("met_vs_emht_JERu").Fill(met_JERu->p.Pt(), emht, weight);
  m2->at("met_vs_emht_JERd").Fill(met_JERd->p.Pt(), emht, weight);
  m2->at("metRaw_vs_emht").Fill(metRaw->p.Pt(), emht, weight);
  m2->at("met_vs_mht").Fill(met->p.Pt(), recoil.Pt(), weight);
  m2->at("memht_vs_emht").Fill(emrecoil.Pt(), emht, weight);
  m2->at("mht_vs_emht").Fill(recoil.Pt(), emht, weight);

} // end fill

map<string,TH1F> initSignalHistograms(unsigned pdfWeightSize=0) {
  map<string,TH1F> hMap;
  hMap["met"] = TH1F("", ";#it{E}_{T}^{miss} (GeV)", 100, 0, 1000);
  for (unsigned i=0; i<pdfWeightSize; i++) {
    string hname = "met_weight_"+to_string(i);
    hMap[hname] = TH1F("", ";#it{E}_{T}^{miss} (GeV)", 100, 0, 1000);
  }
  hMap["met_puUp"] = TH1F("", ";#it{E}_{T}^{miss} (GeV)", 100, 0, 1000);
  hMap["met_puDn"] = TH1F("", ";#it{E}_{T}^{miss} (GeV)", 100, 0, 1000);
  hMap["met_jesUp"] = TH1F("", ";#it{E}_{T}^{miss} (GeV)", 100, 0, 1000);
  hMap["met_jesDn"] = TH1F("", ";#it{E}_{T}^{miss} (GeV)", 100, 0, 1000);
  hMap["met_jerUp"] = TH1F("", ";#it{E}_{T}^{miss} (GeV)", 100, 0, 1000);
  hMap["met_jerDn"] = TH1F("", ";#it{E}_{T}^{miss} (GeV)", 100, 0, 1000);

  return hMap;
}

void HistogramProducer::fillSignalSelection(const string& s, float addWeight=1.)
{
  auto weight = selW*addWeight;
  if (std::isnan(met->p.X()) and std::isnan(met->p.Y())) return;

  if (!h1Maps.count(s)) {
    h1Maps[s] = initSignalHistograms(pdf_weights->size());
  }
  auto m1 = &h1Maps[s];
  auto _met = met->p.Pt();
  m1->at("met").Fill(_met, weight);
  for (unsigned i=0; i<pdf_weights->size(); i++) {
    string hname = "met_weight_"+to_string(i);
    m1->at(hname).Fill(_met, weight*pdf_weights->at(i));
  }
  if (!isData) {
    m1->at("met_puUp").Fill(_met, weight*weighters.at("puWeightUp").getWeight(*nTruePV)/ *pu_weight);
    m1->at("met_puDn").Fill(_met, weight*weighters.at("puWeightDn").getWeight(*nTruePV)/ *pu_weight);
    m1->at("met_jesUp").Fill(met_JESu->p.Pt(), weight);
    m1->at("met_jesDn").Fill(met_JESd->p.Pt(), weight);
    m1->at("met_jerUp").Fill(met_JERu->p.Pt(), weight);
    m1->at("met_jerDn").Fill(met_JERd->p.Pt(), weight);
  }
}


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
  met_JESu(fReader, "met_JESu"),
  met_JESd(fReader, "met_JESd"),
  met_JERu(fReader, "met_JERu"),
  met_JERd(fReader, "met_JERd"),
  nGoodVertices(fReader, "nGoodVertices"),
  nTracksPV(fReader, "nTracksPV"),
  pu_weight(fReader, "pu_weight"),
  mc_weight(fReader, "mc_weight"),
  pdf_weights(fReader, "pdf_weights"),
  genHt(fReader, "genHt"),
  nTruePV(fReader, "true_nPV"),
  runNo(fReader, "runNo"),
  lumNo(fReader, "lumNo"),
  evtNo(fReader, "evtNo"),
  hlt_photon90_ht600(fReader, "HLT_Photon90_CaloIdL_PFHT600_v"),
  hlt_photon90(fReader, "HLT_Photon90_v"),
  hlt_ht600(fReader, "HLT_PFHT600_v"),
  hlt_ht800(fReader, "HLT_PFHT800_v"),
  hlt_el27(fReader, "HLT_Ele27_eta2p1_WPTight_Gsf_v"),
  hlt_ht600_pre(fReader, "HLT_PFHT600_v_pre"),
  looseCutFlowPhoton({{"hoe_eb",0.0597}, {"sigmaIetaIeta_eb",0.01031}, {"cIso_eb",1.295}, {"nIso1_eb",10.910}, {"nIso2_eb",0.0148}, {"nIso3_eb",0.000017}, {"pIso1_eb",3.630}, {"pIso2_eb",0.0047},
    {"hoe_ee",0.0481}, {"sigmaIetaIeta_ee",0.03013}, {"cIso_ee",1.011}, {"nIso1_ee",5.931}, {"nIso2_ee",0.0163}, {"nIso3_ee",0.000014}, {"pIso1_ee",6.641}, {"pIso2_ee",0.0034} }),
  startTime(time(NULL)),
  rand()
{
}

void HistogramProducer::Init(TTree *tree)
{
  fReader.SetTree(tree);
  inputName = fReader.GetTree()->GetCurrentFile()->GetName();
  isData = inputName.find("Run201") != string::npos;
  isSignal = inputName.find("SMS") != string::npos;
  resolution = Resolution(isData? "Spring16_25nsV6_DATA_PtResolution_AK4PFchs.txt": "Spring16_25nsV6_MC_PtResolution_AK4PFchs.txt");
  weighters["sf_photon_id_loose"] = Weighter("../plotter/data/dataMcScaleFactors_80X.root", "EGamma_SF2D");
  weighters["sf_photon_pixel"] = Weighter("../plotter/data/ScalingFactors_80X_Summer16.root", "Scaling_Factors_HasPix_R9 Inclusive");
  string puUp = "pileupWeightUp_mix_2016_25ns_Moriond17MC_PoissonOOTPU";
  string puDn = "pileupWeightDown_mix_2016_25ns_Moriond17MC_PoissonOOTPU";
  if (isSignal) {
    puUp = "pileupWeightUp_mix_2016_25ns_SpringMC_PUScenarioV1_PoissonOOTPU";
    puDn = "pileupWeightDown_mix_2016_25ns_SpringMC_PUScenarioV1_PoissonOOTPU";
  }
  weighters["puWeightUp"] = Weighter("../../CMSSW/treewriter/CMSSW_8_0_25/src/TreeWriter/PUreweighting/data/puWeights.root", puUp);
  weighters["puWeightDn"] = Weighter("../../CMSSW/treewriter/CMSSW_8_0_25/src/TreeWriter/PUreweighting/data/puWeights.root", puDn);
  weighters.at("sf_photon_id_loose").fillOverflow2d();
  weighters.at("sf_photon_pixel").fillOverflow2d();

  float lumi = 35.867e3; // pb^{-1}
  cutFlow = *((TH1F*)fReader.GetTree()->GetCurrentFile()->Get("TreeWriter/hCutFlow"));
  fReader.GetEntries(true); // jumps to last file
  string lastInputName = fReader.GetTree()->GetCurrentFile()->GetName();
  if (inputName!=lastInputName) {
    // adds cut flow of last file. This makes only scence if there is for 1 or two files
    cutFlow.Add((TH1F*)fReader.GetTree()->GetCurrentFile()->Get("TreeWriter/hCutFlow"));
  }
  float nGen = cutFlow.GetBinContent(2);
  sampleW = isData ? 1. : lumi * sampleCrossSection(inputName) / nGen;

  genHt600 = inputName.find("HT-0to600") != string::npos;
  noPromptPhotons = inputName.find("QCD_HT") != string::npos
    || inputName.find("TTJets") != string::npos
    || inputName.find("WJetsToLNu") != string::npos
    || inputName.find("ZJetsToNuNu") != string::npos;

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
      || indexOfMatchedParticle<tree::Photon*>(jet, selPhotons, .4) >= 0
      || jet.p.Pt() < 30 || fabs(jet.p.Eta()) > 3) continue;
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
    weight *= fabs(eta)<photonsEtaMaxBarrel ? 0.964 : 0.94;
  }
  return weight;
}

Bool_t HistogramProducer::Process(Long64_t entry)
{
  resetSelection();
  fReader.SetLocalEntry(entry);

  if (genHt600 && *genHt>600) {
    return kTRUE;
  }
  if (isSignal && unMatchedSuspiciousJet(*jets, *genJets)) {
    for (int i=0;i<cutFlow.GetNbinsX()+2;i++) {
      cutFlow.AddBinContent(i, -1); // this is not considered in nGen for tree weights
    }
    return kTRUE;
  }

  bool cutPrompt = noPromptPhotons && count_if(genParticles->begin(), genParticles->end(), [] (const tree::GenParticle& p) { return p.pdgId==22 && p.promptStatus == DIRECTPROMPT;});
  effMap.at("noPromptEvaluation").Fill(cutPrompt, 0);

  fillUncut();
  fillTriggerStudies();

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

  float emht=0;
  for (auto& p : selPhotons) emht += p->p.Pt();
  for (auto& p : selJets) emht += p->p.Pt();

  if (selPhotons.size() && emht > 700 && (*hlt_photon90_ht600 || !isData)) {
    auto g = selPhotons.at(0);
    auto addWeight = getPhotonWeight(*g);
    auto dPhi = fabs(met->p.DeltaPhi(g->p));
    bool orthogonal = .3<dPhi && dPhi<2.84;
    bool genE = fabs(genMatchNegativePrompt(*g, *genParticles)) == 11;
    bool genE2 = genMatchWZDecay(*g, *intermediateGenParticles);
    bool isGenEclean = isData || isSignal || !genE;

    if (!cutPrompt && orthogonal && isGenEclean && emht<2000) fillSignalSelection("signal_lowEMHT", addWeight);
    if (!cutPrompt && orthogonal && isGenEclean && 2000<emht) fillSignalSelection("signal_highEMHT", addWeight);

    if (!cutPrompt && orthogonal && genE && emht<2000) fillSignalSelection("signal_lowEMHT_genE", addWeight);
    if (!cutPrompt && orthogonal && genE && 2000<emht) fillSignalSelection("signal_highEMHT_genE", addWeight);

    if (!cutPrompt && isGenEclean && emht<2000) fillSignalSelection("signal_lowEMHT2", addWeight);
    if (!cutPrompt && isGenEclean && 2000<emht) fillSignalSelection("signal_highEMHT2", addWeight);

    if (!cutPrompt && genE && emht<2000) fillSignalSelection("signal_lowEMHT2_genE", addWeight);
    if (!cutPrompt && genE && 2000<emht) fillSignalSelection("signal_highEMHT2_genE", addWeight);

    if (g->p.Pt()>300 && met->p.Pt()>200) cout << *runNo << ":" << *lumNo << ":" << *evtNo << endl;

    fillSelection("tr", addWeight, true);
    if (genE) fillSelection("tr_genE", addWeight, true);
    if (!selElectrons.size() && !selMuons.size()) fillSelection("tr_noLep", addWeight);
    if (g->isTight) fillSelection("tr_tight", addWeight);
    if (met->p.Pt() < 100) fillSelection("tr_0met100", addWeight);
    else fillSelection("tr_100met", addWeight);
  }

  if (!selPhotons.size() && emht > 700 && (*hlt_ht600 || !isData)) {
    fillSelection("tr_jControl", *hlt_ht600_pre, true);
    if (!selElectrons.size() && !selMuons.size()) fillSelection("tr_jControl_noLep", *hlt_ht600_pre);
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

  emht=0;
  for (auto& p : selPhotons) emht += p->p.Pt();
  for (auto& p : selJets) emht += p->p.Pt();

  if (selPhotons.size() && emht > 700 && (*hlt_photon90_ht600 || !isData)) {
    auto g = selPhotons.at(0);
    auto addWeight = getPhotonWeight(*g);
    auto dPhi = fabs(met->p.DeltaPhi(g->p));
    bool orthogonal = .3<dPhi && dPhi<2.84;
    bool genE = fabs(genMatchNegativePrompt(*g, *genParticles)) == 11;

    bool isGenEclean = isData || isSignal || !genE;

    if (!cutPrompt && orthogonal && isGenEclean && emht<2000) fillSignalSelection("signal_lowEMHT_ee", addWeight);
    if (!cutPrompt && orthogonal && isGenEclean && 2000<emht) fillSignalSelection("signal_highEMHT_ee", addWeight);

    if (!cutPrompt && isGenEclean && emht<2000) fillSignalSelection("signal_lowEMHT2_ee", addWeight);
    if (!cutPrompt && isGenEclean && 2000<emht) fillSignalSelection("signal_highEMHT2_ee", addWeight);


    fillSelection("tr_ee", addWeight);
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

  emht=0;
  for (auto& p : selPhotons) emht += p->p.Pt();
  for (auto& p : selJets) emht += p->p.Pt();
  if (selPhotons.size() && emht > 700 && (*hlt_photon90_ht600 || !isData)) {
    auto g = selPhotons.at(0);
    auto addWeight = getPhotonWeight(*g);
    auto dPhi = fabs(met->p.DeltaPhi(g->p));
    bool orthogonal = .3<dPhi && dPhi<2.84;

    if (!cutPrompt && orthogonal && emht<2000) fillSignalSelection("signal_lowEMHT_eControl", addWeight);
    if (!cutPrompt && orthogonal && 2000<emht) fillSignalSelection("signal_highEMHT_eControl", addWeight);

    if (!cutPrompt && emht<2000) fillSignalSelection("signal_lowEMHT2_eControl", addWeight);
    if (!cutPrompt && 2000<emht) fillSignalSelection("signal_highEMHT2_eControl", addWeight);


    fillSelection("tr_eControl", addWeight, true);
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

  emht=0;
  for (auto& p : selPhotons) emht += p->p.Pt();
  for (auto& p : selJets) emht += p->p.Pt();
  if (selPhotons.size() && emht > 700 && (*hlt_photon90_ht600 || !isData)) {
    auto g = selPhotons.at(0);
    auto addWeight = getPhotonWeight(*g);
    auto dPhi = fabs(met->p.DeltaPhi(g->p));
    bool orthogonal = .3<dPhi && dPhi<2.84;

    if (!cutPrompt && orthogonal && emht<2000) fillSignalSelection("signal_lowEMHT_ee_eControl", addWeight);
    if (!cutPrompt && orthogonal && 2000<emht) fillSignalSelection("signal_highEMHT_ee_eControl", addWeight);

    if (!cutPrompt && emht<2000) fillSignalSelection("signal_lowEMHT2_ee_eControl", addWeight);
    if (!cutPrompt && 2000<emht) fillSignalSelection("signal_highEMHT2_ee_eControl", addWeight);


    fillSelection("tr_eControl_ee", addWeight);
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
