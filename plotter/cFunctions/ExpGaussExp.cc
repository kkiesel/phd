#include "ExpGaussExp.h"

ClassImp(ExpGaussExp)

ExpGaussExp::ExpGaussExp(const char *name, const char *title,
                        RooAbsReal& _x,
                        RooAbsReal& _p0,
                        RooAbsReal& _p1,
                        RooAbsReal& _p2,
                        RooAbsReal& _p3) :
 RooAbsPdf(name,title),
 x("x","x",this,_x),
 p0("p0","p0",this,_p0),
 p1("p1","p1",this,_p1),
 p2("p2","p2",this,_p2),
 p3("p3","p3",this,_p3)
{ }


ExpGaussExp::ExpGaussExp(const ExpGaussExp& other, const char* name) :
 RooAbsPdf(other,name),
 x("x",this,other.x),
 p0("p0",this,other.p0),
 p1("p1",this,other.p1),
 p2("p2",this,other.p2),
 p3("p3",this,other.p3)
{ }

Double_t ExpGaussExp::evaluate() const
{
  Double_t std = (x-p0)/p1;
  Double_t result = 0;
  if (std>p2) {
    result = exp( p2*p2/2. - p2*std );
  } else if (std<=p2 && std>-p3) {
    result = exp( -0.5*pow(std, 2) );
  } else if (std<=-p3) {
    result = exp( p3*p3/2. + p3*std );
  }
  return result;
}
