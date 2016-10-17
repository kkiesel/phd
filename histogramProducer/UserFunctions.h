#include <algorithm>
#include <regex>
#include <TH2F.h>

// https://twiki.cern.ch/twiki/bin/view/CMS/BtagRecommendation74X
enum bTaggingEnum { CSVv2L, CSVv2M, CSVv2T };
map<bTaggingEnum,float> bTaggingWorkingPoints = {
  { CSVv2L, 0.605 },
  { CSVv2M, 0.89 },
  { CSVv2T, 0.97}
};

pair<TVector3,TVector3> megajets( const vector<TVector3>& jets ) {
  // code from https://twiki.cern.ch/twiki/bin/view/CMSPublic/RazorLikelihoodHowTo

  TVector3 j1, j2;
  int N_comb = (int) pow( 2, jets.size() );

  double M_min = numeric_limits<double>::max();
  int j_count;
  for(int i = 1; i < N_comb-1; i++){
    TVector3 j_temp1, j_temp2;
    int itemp = i;
    j_count = N_comb/2;
    int count = 0;
    while(j_count > 0){
      if(itemp/j_count == 1){
        j_temp1 += jets[count];
      } else {
        j_temp2 += jets[count];
      }
      itemp -= j_count*(itemp/j_count);
      j_count /= 2;
      count++;
    }
    double M_temp = j_temp1.Mag2()+j_temp2.Mag2();
    // smallest mass
    if(M_temp < M_min){
      M_min = M_temp;
      j1 = j_temp1;
      j2 = j_temp2;
    }
  }

  if(j2.Pt() > j1.Pt()) return pair<TVector3,TVector3>(j2,j1);
  else                  return pair<TVector3,TVector3>(j1,j2);
}


pair<float,float> razorVariables( const pair<TVector3,TVector3>& megajets, const TVector3& met ) {
  float mr2 = pow( megajets.first.Mag() + megajets.second.Mag(), 2 ) - pow( megajets.first.Z() + megajets.second.Z(), 2 );
  float mrt2 = 0.5 * ( met.Mag()*( megajets.first.Pt() + megajets.second.Pt() ) - met.XYvector() * ( megajets.first.XYvector() + megajets.second.XYvector() ) );
  float r2 = mrt2 / mr2;
  return pair<float,float>(sqrt(mr2),r2);
}

std::ostream& operator << ( std::ostream& os, const TVector3& p ) {
  os << p.Pt() << "\t" << p.Eta() << "\t" << p.Phi();
  return os;
}


string getOutputFilename( const string& inputFileName, const string& appendix="hists" ) {

  // Converts "/path/to/ntuple/QCD_nTuple.root" to "QCD_hists.root"

  auto startPos = inputFileName.rfind("/");
  auto endPos = inputFileName.find("_nTuple.root");
  string outputName = "out.root";
  if( endPos != string::npos ) {
    outputName = inputFileName.substr( startPos+1, endPos-startPos-1 ) + "_"+appendix+".root";
  }

  // for signal scans
  if( inputFileName.find("/T5") != string::npos ) {
    endPos = inputFileName.find(".root");
    if( endPos != string::npos ) {
      outputName = inputFileName.substr( startPos+1, endPos-startPos-1 ) + "_"+appendix+".root";
    }
  }
  return outputName;

}

template <typename VectorClass>
int indexOfMatchedParticle( const tree::Particle& tag, const std::vector<VectorClass>& particles, float deltaR=.1, float relPt=-1 ) {
  for( int i=0; i<(int)particles.size(); ++i ) {
    auto p = particles.at(i)->p;
    if(    ( deltaR<0 || p.DeltaR( tag.p ) < deltaR )
        && ( relPt<0  || fabs(p.Pt()-tag.p.Pt())/tag.p.Pt() < relPt ) ) {
      return i;
    }
  }
  return -1;
}

vector<int> getRunList( const string& filename="/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions15/13TeV/Cert_246908-260627_13TeV_PromptReco_Collisions15_25ns_JSON.txt" ) {
  ifstream t(filename);
  string s((std::istreambuf_iterator<char>(t)),std::istreambuf_iterator<char>());
  smatch m;
  regex e ("\"([\\d]+)\"");
  vector<int> runList;
  while (regex_search (s,m,e)) {
    runList.push_back( stoi( m[1] ) );
    s = m.suffix().str();
  }
  return runList;
}

std::ostream &operator<<(std::ostream &os, TVector3 &p) {
      return os << p.Pt() << "\t" << p.Eta() << "\t" << p.Phi();
}

inline int pseudoRandom( float f ) {
  //return *(int*)&f; // Interprets float as int
  return ((int)(f*1e6)) % 1000;
}

inline bool pseudoRandomSort(const tree::Particle p1, const tree::Particle p2) {
  return pseudoRandom(p1.p.Phi()) > pseudoRandom(p2.p.Phi());
}

const float photonsEtaMaxBarrel = 1.4442;

class JetSelector {
  public:
    JetSelector( const string& filename, const string& histname ) {
      TFile f( filename.c_str() );
      TH2F* h2 = (TH2F*)f.Get( histname.c_str() );
      if(h2) {
        for(unsigned i=1; i<(unsigned)h2->GetNbinsX()+2; i++) {
          auto h = h2->ProjectionY(("proj"+to_string(i)).c_str(), i, i );
          if(h->Integral()) {
            h->Scale(1./h->Integral());
            hMap[i-1] = *h;
          }
        }
      }
    }

    ~JetSelector(){
    }

    unsigned getJetN( unsigned nJets ){
      if(!hMap.count(nJets)) return 0;
      return (unsigned)std::round(hMap.at(nJets).GetRandom());
    }

  private:
    map<unsigned,TH1D> hMap;
};

float sampleCrossSection(const string& inputFileName) {
  auto startPos = inputFileName.rfind("/");
  auto endPos = inputFileName.find("_nTuple.root");
  string sampleName = "";
  if (endPos != string::npos) {
    sampleName = inputFileName.substr(startPos+1, endPos-startPos-1);
  }

  map<string,float> xs = {
    {"GJets_HT-40To100", 20790},
    {"GJets_HT-100To200", 9238},
    {"GJets_HT-200To400", 2305},
    {"GJets_HT-400To600", 274.4},
    {"GJets_HT-600ToInf", 93.46},
    {"QCD_HT100to200", 27990000},
    {"QCD_HT200to300", 1712000},
    {"QCD_HT300to500", 347700},
    {"QCD_HT500to700", 32100},
    {"QCD_HT700to1000", 6831},
    {"QCD_HT1000to1500", 1207},
    {"QCD_HT1500to2000", 119.9},
    {"QCD_HT2000toInf", 25.24},
    {"TTJets", 670.3},
    {"WJetsToLNu_HT-100To200", 1345.*1.21},
    {"WJetsToLNu_HT-200To400", 359.7*1.21},
    {"WJetsToLNu_HT-400To600", 48.91*1.21},
    {"WJetsToLNu_HT-600To800", 12.05*1.21},
    {"WJetsToLNu_HT-800To1200", 5.501*1.21},
    {"WJetsToLNu_HT-1200To2500", 1.329*1.21},
    {"WJetsToLNu_HT-2500ToInf", 0.03216*1.21},
    {"WJetsToLNu_HT-600ToInf", 18.77*1.21},
    {"TTGJets", 3.697},
    {"TGJets_amcatnlo_madspin", 2.967},
    {"WGToLNuG-amcatnloFXFX", 489.},
    {"WGToLNuG-madgraphMLM", 405.271},
    {"WGToLNuG-madgraphMLM_PtG-0to130", 405.271},
    {"WGToLNuG_PtG-500", 0.0117887},
    {"ZJetsToNuNu_HT-100To200", 280.47*1.23},
    {"ZJetsToNuNu_HT-200To400", 78.36*1.23},
    {"ZJetsToNuNu_HT-400To600", 10.94*1.23},
    {"ZJetsToNuNu_HT-600ToInf", 4.20*1.23},
    {"ZJetsToNuNu_HT-600To800", 0.853*1.23},
    {"ZJetsToNuNu_HT-800To1200", 0.3942*1.23},
    {"ZJetsToNuNu_HT-1200To2500", 0.0974*1.23},
    {"ZJetsToNuNu_HT-2500ToInf", 0.002308*1.23},
    {"ZNuNuGJets_MonoPhoton_PtG-40to130", 2.816},
    {"ZNuNuGJets_MonoPhoton_PtG-130", 0.223},
    {"ZGTo2LG", 117.864},
    {"ZGTo2LGmod", 117.864},
    {"DYJetsToLL_M-50", 6025.2},
    {"WGJets_MonoPhoton_PtG-130", 0.834} // source: jlange
  };
  if (xs.count(sampleName)) {
    return xs.at(sampleName);
  } else {
    cout << "No cross section defined for " << inputFileName << endl;
    return -1;
  }
}

std::pair<float, float> jerScales(float eta) {
  // Copied on 2016-09-14 from https://twiki.cern.ch/twiki/bin/viewauth/CMS/JetResolution
  // for 80X
  eta = std::abs(eta);
  float sf=0, e=0;
  if (eta<.5) {
    sf = 1.122; e = 0.026;
  } else if (eta<.8) {
    sf = 1.167; e = 0.048;
  } else if (eta<1.1) {
    sf = 1.168; e = 0.046;
  } else if (eta<1.3) {
    sf = 1.029; e = 0.066;
  } else if (eta<1.7) {
    sf = 1.115; e = 0.03;
  } else if (eta<1.9) {
    sf = 1.041; e = 0.062;
  } else if (eta<2.1) {
    sf = 1.167; e = 0.086;
  } else if (eta<2.3) {
    sf = 1.094; e = 0.093;
  } else if (eta<2.5) {
    sf = 1.168; e = 0.120;
  } else if (eta<2.8) {
    sf = 1.266; e = 0.132;
  } else if (eta<3.0) {
    sf = 1.595; e = 0.175;
  } else if (eta<3.2) {
    sf = 0.998; e = 0.066;
  } else {
    sf = 1.226; e = 0.145;
  }
  return std::pair<float, float> (sf, e);
}

float smearedPtDataMC(const TVector3& p, const vector<tree::Particle>& genJets, TRandom2& rand) {
  auto sfe = jerScales(p.Eta());
  float oldPt = p.Pt();
  float genPt = 0;
  for (const auto& genJ : genJets) {
    if (p.DeltaR(genJ.p) < .2 && std::abs(oldPt-genJ.p.Pt()) < 3*sfe.second*oldPt) {
      genPt = genJ.p.Pt();
      break;
    }
  }
  float newPt = genPt > 10 ?
      std::max((float)0., genPt + sfe.first*(oldPt-genPt))
    : rand.Gaus(oldPt, sqrt(pow(sfe.first,2)-1) * sfe.second*oldPt);
  return newPt;
}

int genMatchNegativePrompt(const tree::Particle& p, const std::vector<tree::GenParticle>& particles) {
  for (auto const& genP : particles) {
    auto id = fabs(genP.pdgId);
    auto dr = p.p.DeltaR(genP.p);
    auto dpt = p.p.Pt()/genP.p.Pt();
    if (dr < 0.15 && fabs(dpt-1) < 0.15) {
      if (genP.isPrompt) return id;
      else return -id;
    }
  }
  return 0;
}

int genMatchWZDecay(const tree::Particle& p, const std::vector<tree::IntermediateGenParticle>& particles) {
  // match to daghters of massive particles
  for (auto const& genP : particles) {
    for (auto const & d : genP.daughters) {
      auto id = fabs(d.pdgId);
      auto dr = p.p.DeltaR(d.p);
      auto dpt = p.p.Pt()/d.p.Pt();
      if (dr < 0.4) return id;
    }
  }
  return 0;
}
