#include<regex>
#include<TH2F.h>

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
    {"WGToLNuG_PtG-500", 0.0117887},
    {"ZJetsToNuNu_HT-100To200", 280.47*1.23},
    {"ZJetsToNuNu_HT-200To400", 78.36*1.23},
    {"ZJetsToNuNu_HT-400To600", 10.94*1.23},
    {"ZJetsToNuNu_HT-600ToInf", 4.20*1.23},
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


