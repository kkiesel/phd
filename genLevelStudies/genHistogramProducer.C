#include<MyGenParticle.h>
#include<TH3F.h>
#include<TH2F.h>
#include<TChain.h>
#include<TFile.h>
#include<TVector3.h>
#include<TCanvas.h>
#include<TSystem.h>
#include<iostream>
#include <algorithm> // sort

using namespace std;

ostream& operator<< (ostream& os, const TVector3& p) {
  os << p.Pt() << "\t" << p.Eta() << "\t" << p.Phi();
  return os;
}

string getProcessString( const TVector3& g, const vector<gen::MyGenParticle>& jets, int idMother1, int idMother2 ) {
  unsigned nGluon = 0;
  unsigned nQuark = 0;
  for( auto& j : jets ){
    if(j.id==21) nGluon++;
    else if( fabs(j.id)<6.1 ) nQuark ++;
    else cout << "No idea what to do with " << j.id << endl;
  }
  string out;
  if( idMother1 == 21 && idMother2 == 21 ) out = "gg";
  else if( idMother1 ==21 || idMother2 == 21 ) out = "gq";
  else out = "qq";
  out.append("TO");
  for( unsigned i=0;i<nGluon;i++) out.append("g");
  for( unsigned i=0;i<nQuark;i++) out.append("q");
  if( g.Pt()>10 ) out.append("G");
  return out;
}

string getOutputName( const string& inputName ) {
  return "hists"+inputName.substr(inputName.find("_"),inputName.back());
}

bool isHeJet( const gen::MyGenParticle& jet ) {
  return jet.p.Pt()>100 && fabs(jet.p.Eta())<1.4442;
}

bool isNotHeJet( const gen::MyGenParticle& jet ) {
  return !isHeJet(jet);
}


bool isAnJet( const gen::MyGenParticle& jet ) {
  return jet.p.Pt()>40 && fabs(jet.p.Eta())<3;
}

class PhotonPicker{
 public:
  PhotonPicker(int nMax=5) :
   hG("",";nPart;gPos", nMax, .5, nMax+.5, nMax, .5, nMax+.5 ),
   hQ("",";nPart", nMax, .5, nMax+.5 )
  {}

  int getGPos( const TVector3& g, const vector<gen::MyGenParticle>& jets ) {
    int pos=1;
    for( auto& j:jets ) {
      if(j.p.Pt()>g.Pt() ) pos++;
    }
    return pos;
  }

  void fillG( const TVector3& g, const vector<gen::MyGenParticle>& jets, const float xsec=0 ) {
    int pos = getGPos( g, jets );
    hG.Fill( jets.size()+1, pos, xsec );
  }

  void fillQ( unsigned nJets, const float xsec=0 ) {
    hQ.Fill( nJets, xsec );
  }

  void divide(){
    ratio = *hG.ProjectionX();
    ratio.Scale(1./ratio.Integral(0,-1));
    hQ.Scale(1./hQ.Integral(0,-1));
    ratio.Divide( &hQ );
    for( int i=1;i<hG.GetNbinsY()+1;i++) {
      auto htmp = hG.ProjectionY(("_py"+to_string(i)).c_str(),i,i);
      htmp->Scale(1./htmp->Integral());
      vecGpos.push_back(*htmp);
    }
  }

  float getWeight( unsigned nJets ) {
    auto bin = ratio.FindFixBin( nJets );
    return ratio.GetBinContent(bin);
  };

  float getPos( unsigned nJet ) {
    return vecGpos.at(nJet-1).GetRandom();
  }
  void write(const string& fname="picker.root") {
    TFile* f = TFile::Open(fname.c_str(), "recreate");
    hG.Write("hG");
    hQ.Write("hQ");
    f->Close();
  }
  void read(const string& fname="picker.root") {
    TFile* f = TFile::Open(fname.c_str());
    hG = *(TH2F*)f->Get("hG");
    hQ = *(TH1D*)f->Get("hQ");
  }
 private:
  TH2F hG;
  TH1D hQ;
  vector<TH1D> vecGpos;
  TH1D ratio;
};


class HistogramBooker {
 public:
  HistogramBooker(){
    defaultH["ptO1"] = TH1F("",";p_{T}", 100, 0, 1000 );
    defaultH["ptO1G"] = TH1F("",";p_{T}", 100, 0, 1000 );
    defaultH["ptO1J"] = TH1F("",";p_{T}", 100, 0, 1000 );
    defaultH["ptO2"] = TH1F("",";p_{T}", 100, 0, 1000 );
    defaultH["ptO2G"] = TH1F("",";p_{T}", 100, 0, 1000 );
    defaultH["ptO2J"] = TH1F("",";p_{T}", 100, 0, 1000 );
    defaultH["ptO3"] = TH1F("",";p_{T}", 100, 0, 1000 );
    defaultH["ptO3G"] = TH1F("",";p_{T}", 100, 0, 1000 );
    defaultH["ptO3J"] = TH1F("",";p_{T}", 100, 0, 1000 );
    defaultH["ptO4"] = TH1F("",";p_{T}", 100, 0, 1000 );
    defaultH["ptO4G"] = TH1F("",";p_{T}", 100, 0, 1000 );
    defaultH["ptO4J"] = TH1F("",";p_{T}", 100, 0, 1000 );
    defaultH["ptO5"] = TH1F("",";p_{T}", 100, 0, 1000 );
    defaultH["ptO5G"] = TH1F("",";p_{T}", 100, 0, 1000 );
    defaultH["ptO5J"] = TH1F("",";p_{T}", 100, 0, 1000 );
    defaultH["ptG"] = TH1F("",";p_{T}", 100, 0, 1000 );
    defaultH["ptG0"] = TH1F("",";p_{T}", 100, 0, 1000 );
    defaultH["ptG1"] = TH1F("",";p_{T}", 100, 0, 1000 );
    defaultH["ptG2"] = TH1F("",";p_{T}", 100, 0, 1000 );
    defaultH["ptG3"] = TH1F("",";p_{T}", 100, 0, 1000 );
    defaultH["ptG4"] = TH1F("",";p_{T}", 100, 0, 1000 );
    defaultH["ptJ1"] = TH1F("",";p_{T}", 100, 0, 1000 );
    defaultH["ptJ2"] = TH1F("",";p_{T}", 100, 0, 1000 );
    defaultH["ptJ3"] = TH1F("",";p_{T}", 100, 0, 1000 );
    defaultH["ptJ4"] = TH1F("",";p_{T}", 100, 0, 1000 );
    defaultH["etaG"] = TH1F("",";|#eta|", 5, 0, 1.4442 );
    defaultH["etaJ1"] = TH1F("",";|#eta|", 5, 0, 1.4442 );
    defaultH["etaJ2"] = TH1F("",";|#eta|", 5, 0, 1.4442 );
    defaultH["etaJ3"] = TH1F("",";|#eta|", 5, 0, 1.4442 );
    defaultH["etaJ4"] = TH1F("",";|#eta|", 5, 0, 1.4442 );
    defaultH["nHeJet"] = TH1F("",";n_{he jet}", 5, 0.5, 5.5 );
    defaultH["nJet"] = TH1F("",";n_{jet}", 5, 0.5, 5.5 );
    defaultH["nAll"] = TH1F("",";n_{jet}", 5, 0.5, 5.5 );
    defaultH["nPart"] = TH1F("",";n_{he jet+#gamma}", 5, 0.5, 5.5 );
    defaultH["posG"] = TH1F("",";#gamma position", 5, 0.5, 5.5 );
    defaultH["ht"] = TH1F("",";H_{T}", 100, 0, 2000 );
    defaultH["emht"] = TH1F("",";EMH_{T}", 100, 0, 2000 );
    defaultH["heht"] = TH1F("",";he H_{T}", 100, 0, 2000 );
    defaultH["heemht"] = TH1F("",";he EMH_{T}", 100, 0, 2000 );
    defaultH["met"] = TH1F("","EMH_{T}^{miss}", 100, 0, 1000 );
    defaultH["hemet"] = TH1F("","he EMH_{T}^{miss}", 100, 0, 1000 );

  } // end definition

  void fill( const string& name, const TVector3& g, const vector<gen::MyGenParticle>& jets, float xsec ) {
    if(!hMapMap.count(name)) hMapMap[name] = defaultH;

    float ht=0,heht=0;
    int njet=0, nhejet=0,higherJets=0;
    TVector3   met(0,0,0);
    TVector3 hemet(0,0,0);

    for(unsigned i=0;i<jets.size();i++) {
      float pt = jets.at(i).p.Pt();
      if( isHeJet(jets.at(i)) ) {
        nhejet++;
        heht += pt;
        hMapMap[name]["ptJ"+to_string(i+1)].Fill( pt, xsec );
        hMapMap[name]["etaJ"+to_string(i+1)].Fill( fabs(jets.at(i).p.Eta()), xsec );
        hemet -= jets.at(i).p;
        if(pt>g.Pt()) higherJets++;
      }
      if( isAnJet(jets.at(i)) ) {
        njet++;
        ht += pt;
        met -= jets.at(i).p;
      }
    }
    hMapMap[name]["nJet"].Fill( njet, xsec );
    hMapMap[name]["nAllJet"].Fill( jets.size(), xsec );
    hMapMap[name]["nHeJet"].Fill( nhejet, xsec );
    hMapMap[name]["ht"].Fill( ht, xsec );
    hMapMap[name]["heht"].Fill( heht, xsec );

    if( g.Pt()>10 ) { // only for gjet
      hMapMap[name]["ptG"+to_string(higherJets)].Fill(g.Pt(), xsec );

      auto objects = jets;
      objects.erase( remove_if(objects.begin(), objects.end(), isNotHeJet), objects.end());
      gen::MyGenParticle photonJet;
      photonJet.p = g;
      photonJet.id = 22;
      objects.push_back(photonJet);
      sort(objects.begin(),objects.end(),gen::PtGreater);

      for( unsigned i=0;i<objects.size();i++) {
        hMapMap[name]["ptO"+to_string(i+1)].Fill( objects.at(i).p.Pt(), xsec );
        if(objects.at(i).id == 22 ) hMapMap[name]["ptO"+to_string(i+1)+"G"].Fill( objects.at(i).p.Pt(), xsec );
        else hMapMap[name]["ptO"+to_string(i+1)+"J"].Fill( objects.at(i).p.Pt(), xsec );
      }
      unsigned photonPosition = 1; //starting at 1
      for(auto j:jets) {
        if( j.p.Pt() > g.Pt() ) photonPosition++;
      }
      hMapMap[name]["posG"].Fill( photonPosition, xsec );
      hMapMap[name]["ptG"].Fill( g.Pt(), xsec );
      hMapMap[name]["etaG"].Fill( fabs(g.Eta()), xsec );
      ht += g.Pt();
      heht += g.Pt();
      met -= g;
      hemet -= g;
      nhejet++;
    }
    hMapMap[name]["emht"].Fill( ht, xsec );
    hMapMap[name]["heemht"].Fill( heht, xsec );
    hMapMap[name]["nPart"].Fill( nhejet, xsec );
    hMapMap[name]["met"].Fill( met.Pt(), xsec );
    hMapMap[name]["hemet"].Fill( hemet.Pt(), xsec );

  } // end fill

  void fillProcessString( const string& name, const string& processStr, float xsec ) {
    if(!processMapMap.count(name)) processMapMap[name] = map<string,float>();
    if(!processMapMap.at(name).count(processStr)) processMapMap.at(name)[processStr] = 0;
    processMapMap.at(name).at(processStr) += xsec;
  } // end fillProcessString

  void save(const string& filename="outHists.root"){
    TFile f(filename.c_str(), "recreate" );
    f.cd();
    for(auto& mapIt : hMapMap ) {
      for(auto& mapItIt : mapIt.second ) {
        auto integral = mapItIt.second.Integral(0,-1);
        if( integral < 1e-6 ) continue;
        mapItIt.second.SetDrawOption("hist");
        //mapItIt.second.Scale( 1./integral );
        //mapItIt.second.Scale( 1./integral );
        mapItIt.second.SetLineColor( mapIt.first.find("qcd")!=string::npos ? 2 : 1 );
        mapItIt.second.Write( (mapIt.first+"_"+mapItIt.first).c_str(), TObject::kWriteDelete );

      }
    }
    for(auto& mapIt : processMapMap ) {
      TH1F h1("","",mapIt.second.size(),0,mapIt.second.size());
      unsigned bin=1;
      for(auto& mapItIt : mapIt.second ) {
        h1.GetXaxis()->SetBinLabel(bin,mapItIt.first.c_str());
        h1.SetBinContent(bin,mapItIt.second);
        bin++;
      }
      h1.Write( ("processes_"+mapIt.first).c_str(), TObject::kWriteDelete );
    }
    cout << "Written to file " << filename << endl;
  }

 private:
  map<string,map<string,TH1F>> hMapMap;
  map<string,TH1F> defaultH;
  map<string,map<string,float>> processMapMap;

};

void doThings(){
  string gjetName = "/user/kiesel/nTuples/genTreesV03/GJets_nTuple.root";
  string qcdName = "/user/kiesel/nTuples/genTreesV03/QCD_nTuple.root";
  const char* treeName = "GenWriter/genParticleTree";

  treeName = "genParticleTree";
  // all100
  gjetName = "GJets_nTuple_all100.root";
  qcdName = "QCD_nTuple_all100.root";

  // an like
  gjetName = "GJets_slimmedTreeAN.root";
  qcdName = "QCD_slimmedTreeAN.root";

  TChain gChain(treeName);
  gChain.AddFile(gjetName.c_str());

  TChain qChain(treeName);
  qChain.AddFile(qcdName.c_str());


  PhotonPicker pp;

  // qcd tree
  double xsec;
  TVector3* g=0;
  vector<gen::MyGenParticle>* jets=0;
  int idMother1;
  int idMother2;

  qChain.SetBranchAddress("xsec",&xsec);
  gChain.SetBranchAddress("xsec",&xsec);
  qChain.SetBranchAddress("g",&g);
  gChain.SetBranchAddress("g",&g);
  qChain.SetBranchAddress("jets",&jets);
  gChain.SetBranchAddress("jets",&jets);
  gChain.SetBranchAddress("idMother1",&idMother1);
  qChain.SetBranchAddress("idMother1",&idMother1);
  gChain.SetBranchAddress("idMother2",&idMother2);
  qChain.SetBranchAddress("idMother2",&idMother2);

  /*
  // compute weights
  for( int i=0;i<gChain.GetEntries();i++ ){
    gChain.GetEntry(i);
    pp.fillG( *g, *jets, xsec );
  }
  for( int i=0;i<qChain.GetEntries();i++ ){
    qChain.GetEntry(i);
    pp.fillQ( jets->size(), xsec );
  }
  cout << "write" << endl;
  pp.write();
  pp.read();

  pp.divide();

  */
  // plotting
  HistogramBooker h;

  for( int i=0;i<gChain.GetEntries();i++ ){
    if(!(i%20000)) { cout.flush(); cout << "GJets: " << 100.*i/gChain.GetEntries() << " %\r"; }
    gChain.GetEntry(i);
    sort(jets->begin(),jets->end(),gen::PtGreater);
    string processStr = getProcessString(*g,*jets,idMother1,idMother2);
    h.fillProcessString( "gjet", processStr, xsec );
    int nHeJets = count_if( jets->begin(),jets->end(), isHeJet );
    h.fill( "gjet", *g, *jets, xsec );
    h.fill( "gjet_"+to_string(nHeJets+1), *g, *jets, xsec );
    h.fill( "gjet_"+processStr, *g, *jets, xsec );
  }
  cout << endl;

  for( int i=0;i<qChain.GetEntries();i++ ){
    if(!(i%20000)) { cout.flush(); cout << "QCD:  " << 100.*i/qChain.GetEntries() << " %\r"; }
    qChain.GetEntry(i);
    sort(jets->begin(),jets->end(),gen::PtGreater);
    string processStr = getProcessString(*g,*jets,idMother1,idMother2);
    h.fillProcessString( "qcd", processStr, xsec );
    int nHeJets = count_if( jets->begin(),jets->end(), isHeJet );
    h.fill( "qcd", *g, *jets, xsec );
    h.fill( "qcd_"+to_string(nHeJets), *g, *jets, xsec );
    h.fill( "qcd_"+processStr, *g, *jets, xsec );
  }
  cout << endl;
  h.save(getOutputName(gjetName));
}



void genHistogramProducer() {
  gSystem->Load("pluginGenWriterGenWriterAuto");
  doThings();
}
