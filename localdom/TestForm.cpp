#include "IForm.h"
#include <iostream>
using namespace std;

int main()
{
  ImaginaryForm Form;
  Form.init(1.,69.,.02,0.,0.,4,1,1.65,60.);

  double energy;
  cin >> energy;
  double y1 = Form.DeltaAsymmetricHole(energy-.5);
  double y2 = Form.DeltaAsymmetricHole(energy+.5);
  double dy = y2-y1;
  cout << dy << " " << Form.DerDeltaAsymmetricHole(energy);
  return 0;
}
