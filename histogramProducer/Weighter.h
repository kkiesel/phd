
class Weighter {
  public:
    Weighter() {}
    Weighter(const string& filename, const string& histname) {
      TFile f(filename.c_str());
      h = (TH1F*) f.Get(histname.c_str());
      if (h) {
        h->SetDirectory(0);
      } else {
        cerr << "Error in <Weighter::Weighter>:  Could not find histogram "
             << histname << " in file " << filename << endl;
      }
      f.Close();
    }

    // deleting the histogram pointer leads to segmentation violations
    ~Weighter() {}

    float getWeight(float value) {
      return h ? h->GetBinContent( h->FindFixBin(value) ) : 1.;
    }

  private:
    TH1F* h = 0;
};
