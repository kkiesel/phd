#ifndef EXPGAUSSEXP
#define EXPGAUSSEXP

#include <Riostream.h>
#include <math.h>

#include "RooAbsCategory.h"
#include "RooAbsPdf.h"
#include "RooAbsReal.h"
#include "RooCategoryProxy.h"
#include "RooRealProxy.h"
#include "TMath.h"

class ExpGaussExp : public RooAbsPdf {
public:
  ExpGaussExp() {}; 
  ExpGaussExp(const char *name, const char *title,
   RooAbsReal& _x,
   RooAbsReal& _p0,
   RooAbsReal& _p1,
   RooAbsReal& _p2,
   RooAbsReal& _p3);
  ExpGaussExp(const ExpGaussExp& other, const char* name=0);
  virtual TObject* clone(const char* newname) const { return new ExpGaussExp(*this,newname); }
  inline virtual ~ExpGaussExp() { }

protected:
  RooRealProxy x ;
  RooRealProxy p0 ;
  RooRealProxy p1 ;
  RooRealProxy p2 ;
  RooRealProxy p3 ;
  
  Double_t evaluate() const ;

private:
  ClassDef(ExpGaussExp,1) // Your description goes here...
};
 
#endif
