#include "TreeParticles.hpp"
#include <map>

class CutFlowPhoton {
  public:

    CutFlowPhoton( std::map<std::string,float> s ) {
      cut_sigmaIetaIeta_eb = s.at("sigmaIetaIeta_eb");
      cut_cIso_eb = s.at("cIso_eb");
      cut_nIso1_eb = s.at("nIso1_eb");
      cut_nIso2_eb = s.at("nIso2_eb");
      cut_nIso3_eb = s.at("nIso3_eb");
      cut_pIso1_eb = s.at("pIso1_eb");
      cut_pIso2_eb = s.at("pIso2_eb");

      cut_sigmaIetaIeta_ee = s.at("sigmaIetaIeta_ee");
      cut_cIso_ee = s.at("cIso_ee");
      cut_nIso1_ee = s.at("nIso1_ee");
      cut_nIso2_ee = s.at("nIso2_ee");
      cut_nIso3_ee = s.at("nIso3_ee");
      cut_pIso1_ee = s.at("pIso1_ee");
      cut_pIso2_ee = s.at("pIso2_ee");
    }

    bool pass(){
      return pass_pt
        && pass_eta
        && pass_hoe
        && pass_sie
        && pass_cIso
        && pass_nIso
        && pass_pIso;
    }


    bool check( const tree::Photon& photon ){

      float pt = photon.p.Pt();
      float aEta = fabs( photon.p.Eta() );

      // kinematic criteria
      pass_pt = pt > 15;
      pass_EB = aEta < 1.4442;
      pass_EE = 1.566 < aEta && aEta < 2.5;
      pass_eta = pass_EB || pass_EE;

      // identification & isolation
      pass_hoe = photon.hOverE < 0.05;
      pass_sie = false;
      pass_cIso = false;
      pass_nIso = false;
      pass_pIso = false;

      if( pass_EB ) {
        pass_sie = photon.sigmaIetaIeta < cut_sigmaIetaIeta_eb;
        pass_cIso = photon.cIso < cut_cIso_eb;
        pass_nIso = photon.nIso < cut_nIso1_eb + cut_nIso2_eb * pt + cut_nIso3_eb * pow(pt,2);
        pass_pIso = photon.pIso < cut_pIso1_eb + cut_pIso2_eb * pt;
      } else if( pass_EE ) {
        pass_sie = photon.sigmaIetaIeta < cut_sigmaIetaIeta_ee;
        pass_cIso = photon.cIso < cut_cIso_ee;
        pass_nIso = photon.nIso < cut_nIso1_ee + cut_nIso2_ee * pt + cut_nIso3_ee * pow(pt,2);
        pass_pIso = photon.pIso < cut_pIso1_ee + cut_pIso2_ee * pt;
      }

      return pass();
    };

    bool passPt(){ return pass_pt; };
    bool passEta(){ return pass_eta; };
    bool passEB(){ return pass_EB; };
    bool passEE(){ return pass_EE; };
    bool passHoe(){ return pass_hoe; };
    bool passSie(){ return pass_sie; };
    bool passCIso(){ return pass_cIso; };
    bool passNIso(){ return pass_nIso; };
    bool passPIso(){ return pass_pIso; };
    bool passIso(){ return pass_nIso && pass_cIso && pass_pIso; };


 private:
    bool pass_pt;
    bool pass_eta;
    bool pass_EB;
    bool pass_EE;
    bool pass_hoe;
    bool pass_sie;
    bool pass_cIso;
    bool pass_nIso;
    bool pass_pIso;

    float cut_sigmaIetaIeta_eb;
    float cut_cIso_eb;
    float cut_nIso1_eb;
    float cut_nIso2_eb;
    float cut_nIso3_eb;
    float cut_pIso1_eb;
    float cut_pIso2_eb;

    float cut_sigmaIetaIeta_ee;
    float cut_cIso_ee;
    float cut_nIso1_ee;
    float cut_nIso2_ee;
    float cut_nIso3_ee;
    float cut_pIso1_ee;
    float cut_pIso2_ee;


};

//loose = CutFlowPhoton( 0.0103, 2.44, 2.57, 0.0044, 0.5809, 1.92, 0.0043, 0.0277, 1.84, 4.00, 0.0040, 0.9402, 2.15, 0.0041 );
//tight = CutFlowPhoton( 0.0100, 0.91, 0.33, 0.0044, 0.5809, 0.61, 0.0043, 0.0267, 0.65, 0.93, 0.0040, 0.9402, 0.54, 0.0041 );
