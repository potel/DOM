{

  ifstream fileIn ("nca40_48.data");
  ofstream fileOut ("new.data");
 
  string out;
  //  getline(fileIn,out);
  //fileOut << out <<endl;
  getline(fileIn,out);
  fileOut << out <<endl;
  int N;
  fileIn >> N;
  fileOut << N << endl;

  double E,y,dy;
  for (int i=0;i<N;i++)
    {
      fileIn >> E >> y >> dy;
      //y*= 1000.;
      //dy*= 1000.;
      y -=.08;
      fileOut << E << " " << y <<" " << dy << endl;
    }
  fileIn.close();
  fileOut.close();
}
