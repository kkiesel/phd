
class Weighter {
  public:
    Weighter( const string& filename, const string& histname ) {
      TFile f( filename.c_str() );
      h = (TH1F*)f.Get( histname.c_str() );
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
