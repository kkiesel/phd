
class Weighter {
  public:
    Weighter( const string& filename, const string& histname ) {
      TFile f( filename.c_str() );
      h = (TH1F*)f.Get( histname.c_str() );
      if(!h) {
        cerr << "Error in <Weighter>:  Could not find histogram "
             << histname << " in file " << filename << endl;
        throw 5;
      }
      h->SetDirectory(0);
    }

    ~Weighter(){
      delete h;
    }

    float getWeight( float value ){
      return h->GetBinContent( h->FindFixBin( value ) );
    }

  private:
    TH1F* h;
};
