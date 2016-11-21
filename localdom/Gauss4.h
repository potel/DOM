//*** 4 point Gaussian quadrature
const double GaussIntegx4[]={
  -0.861136311594053,
  -0.339981043584856,
   0.339981043584856,
   0.861136311594053
}; 
const  double GaussIntegw4[] = {
  0.347854845137454,
  0.652145154862546,
  0.652145154862546,
  0.347854845137454
};

const GaussInteg Gauss4(4,GaussIntegx4,GaussIntegw4);
