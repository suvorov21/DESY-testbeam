#ifndef LINE
#define LINE

#include <iostream>

#include "TVector3.h"

/// Class describing the line gometry. Used in crossing reconstruction
class TLine_att {
public:
  /// a constructor given two dots
  TLine_att(const TVector3& dot1, const TVector3& dot2, bool vertical = false) {
    double k, a, b, c;
    if (!vertical) {
      k = dot2.Y() - dot1.Y();

      a = (dot2.X() - dot1.X()) / k;
      b = 1.;
      c = (dot2.Z() - dot1.Z()) / k;
    } else {
      k = dot2.X() - dot1.X();

      a  = 1.;
      b = (dot2.Y() - dot1.Y()) / k;
      c = (dot2.Z() - dot1.Z()) / k;
    }

    dir = TVector3(a, b, c);
    dir = dir.Unit();

    pos = dot1;
  }

  double GetDistance(const TVector3& dot) {
    TVector3 hyp = dot - pos;

    float cos = hyp.Dot(dir);
    cos /= hyp.Mag();
    cos /= dir.Mag();

    return hyp.Mag() * sqrt(1 - cos * cos);
  }

  TVector3 GetDistVec(const TVector3& dot) {
    TVector3 hyp = dot - pos;

    TVector3 norm_point = pos + dir * (dir.Dot(hyp));

    return dot - norm_point;
  }

  Double_t GetDistX(const TVector3& dot) {
    return GetDistVec(dot).X();
  }

  Double_t GetDistY(const TVector3& dot) {
    return GetDistVec(dot).Y();
  }

  Double_t GetDistZ(const TVector3& dot) {
    return GetDistVec(dot).Z();
  }

  TVector3 EvalY(double y) {
    double k = (y - pos.Y()) / dir.Y();
    return {pos.X() + k * dir.X(),
            pos.Y() + k * dir.Y(),
            pos.Z() + k * dir.Z()};
  }

  TVector3 EvalX(double x) {
    double k = (x - pos.X()) / dir.X();
    return {pos.X() + k * dir.X(),
            pos.Y() + k * dir.Y(),
            pos.Z() + k * dir.Z())};
  }

  TVector3 GetPos() {
    return pos;
  }

  TVector3 GetDir() {
    return dir;
  }

  void Print() {
    std::cout << "TLine_att: pos:dir" << std::endl;
    pos.Print();
    dir.Print();
  }
private:

  TVector3 pos;
  TVector3 dir;
};

#endif