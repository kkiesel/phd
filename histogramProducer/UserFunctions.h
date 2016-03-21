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

