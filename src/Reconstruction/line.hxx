#ifndef LINE
#define LINE

#include "TVector3.h"

/// Class describing the line gometry. Used in crossing reconstruction
class TLine_att {
public:
  /// a constructor given two dots
  TLine_att(const TVector3& dot1, const TVector3& dot2, bool vertical = false) {
    float k, a, b, c;
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

  float GetDistance(const TVector3& dot) {
    TVector3 hyp = dot - pos;

    float cos = hyp.Dot(dir);
    cos /= hyp.Mag();
    cos /= dir.Mag();

    return hyp.Mag() * sqrt(1 - cos * cos);
  }

  TVector3 GetDistVec(const TVector3& dot) {
    TVector3 hyp = dot - pos;

    TVector3 norm_point = pos + dir * (dir.Dot(hyp));

    return hyp - norm_point;
  }

  float GetDistX(const TVector3& dot) {
    return GetDistVec(dot).X();
  }

  float GetDistY(const TVector3& dot) {
    return GetDistVec(dot).Y();
  }

  float GetDistZ(const TVector3& dot) {
    return GetDistVec(dot).Z();
  }

  TVector3 EvalY(float y) {
    double k = (y - pos.Y()) / dir.Y();
    return TVector3(pos.X() + k * dir.X(),
                    pos.Y() + k * dir.Y(),
                    pos.Z() + k * dir.Z());
  }

  TVector3 EvalX(float x) {
    double k = (x - pos.X()) / dir.X();
    return TVector3(pos.X() + k * dir.X(),
                    pos.Y() + k * dir.Y(),
                    pos.Z() + k * dir.Z());
  }

  TVector3 GetPos() {
    return pos;
  }

  TVector3 GetDir() {
    return dir;
  }

  /// a default destructor
  virtual ~TLine_att() {};
private:

  TVector3 pos;
  TVector3 dir;
};

#endif