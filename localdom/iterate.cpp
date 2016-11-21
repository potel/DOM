// to run this program 
// chisq X  to fit data
// chisq a  to not fit

#include "fit.h"


int main(int Narg, char* pszArg[])
{

  if (Narg > 3) return 0;
  string title(pszArg[1]);
  string fileName = title+".inp";

  //read in input file to cound how many fitted parameters
  ifstream file(fileName.c_str());
  if (file.fail())
    {
      cout << "could not open file " << fileName << endl;
      abort();
    }
  char line[80];
  //skip two lines
  file.getline(line,80);
  string variable;
  double parameter,varied,squared,scaling;
  int Ndim = 0;
  for (int i=0;i<fit::TotPara;i++)
    {
      file >> parameter >> varied >> squared>>scaling >> variable;
      if (varied) Ndim++;
    }  
  

  file.close();
  file.clear();




  fit Fit(&title,Ndim);
  
  for (int i=0;i<Fit.Nreact;i++)
    {
      cout << "here " << endl;
      Fit.React[i]->IterateFermi(); 
      cout << "there " << endl;
      Fit.React[i]->OpenRootFile();
      cout << "PlotFit" << endl;
      Fit.React[i]->PlotFit();
      cout << "plotPot" << endl;
      Fit.React[i]->PlotPotentialEcm();
      Fit.React[i]->WriteAHF();
    }

 
  return 0;
}

