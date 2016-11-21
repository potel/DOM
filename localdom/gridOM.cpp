// to run this program 
// chisq X  to fit data
// chisq a  to not fit

#include "fitOM.h"


int main(int Narg, char* pszArg[])
{

  if (Narg > 5) return 0;
  string title(pszArg[1]);
  string reactionTitle(pszArg[2]);
  string fileName = title+".inp";
  int jdata = atoi(pszArg[3]);

  //read in input file to cound how many fitted parameters
  ifstream file(fileName.c_str());
  if (file.fail())
    {
      cout << "could not open file " << fileName << endl;
      abort();
    }
  char line[80];
  //skip two lines
  //file.getline(line,80);
  string variable;
  double parameter,varied,squared,scaling;
  int Ndim = 0;
  for (int i=0;i<fitOM::TotPara;i++)
    {
      file >> parameter >> varied >> squared>>scaling >> variable;
      if (varied) Ndim++;
    }  
  

  file.close();
  file.clear();




  fitOM Fit(&title,reactionTitle,jdata,Ndim);



  
  double * para;
  para = new double [Ndim];

  for (int i=0;i<Ndim;i++)
    {
      int j= Fit.map2[i];
      para[i] = Fit.allPara[j]*Fit.scaling[j];
    }

  double* xi[Ndim];
  for (int i=0;i<Fit.ND;i++) 
    {
      xi[i] = new double [Ndim];
      for (int j=0;j<Fit.ND;j++)
	{
	  if (i == j) xi[i][j] = 1.;
	  else xi[i][j] = 0.;
	}
    }
  double const ftol = 1.e-2;
  double chiMin=0.;

  string toFit(pszArg[4]);  
  //string toFit("X");
  cout << toFit << endl;
  if (toFit == "X" || toFit ==  "x")chiMin=Fit.powell(para,xi,ftol);

  cout << "minimum chisq= " << chiMin << endl;

  Fit.SavePara(para);

  Fit.WriteFit(chiMin);

  for (int i=0;i<Fit.Nreact;i++)
    {
     Fit.React[i]->OpenRootFile();
     cout << "PlotFit" << endl;
      Fit.React[i]->PlotFit();
     Fit.React[i]->CloseRootFile();
     file.close();
     file.clear();
    }

  cout << " for Ecm = "<< Fit.React[0]->data[0].energyCM << " MeV" << endl;
  cout << " for Elab = "<< Fit.React[0]->data[0].energyLab << " MeV" << endl;
  if (Fit.React[0]->Zp == 0) 
    cout << " predicted= " << Fit.React[0]->scatter.TotXsection() << 
    "exp = " << Fit.React[0]->TotXdata[0].xsec << " mb "<< endl;

  Fit.React[0]->scatter.VIntegrals();
  cout << " JReal= " << Fit.React[0]->scatter.JReal 
       << " JImag= " << Fit.React[0]->scatter.JImag 
       << " RrmsReal= " << Fit.React[0]->scatter.RrmsReal 
       << " RrmsImag= " << Fit.React[0]->scatter.RrmsImag << endl;

  return 0;
}

