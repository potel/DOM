{
  TF1 * fun = new TF1("funct","51./pow(.821*sqrt(3.14159),3)*exp(-pow(x/.821,2))+20.48/pow(1.195*sqrt(3.14159),3)*exp(-pow(x/1.195,2))",-3,3);
  fun->Draw();


  TF1 * fun2 = new TF1("funct2","1.13*51./pow(.821*sqrt(3.14159),3)*exp(-pow(x/.821,2))",-3,3);
  fun2->SetLineColor(2);
  fun2->Draw("same");

  TF1 * fun3 = new TF1("funct3","20.48/pow(1.195*sqrt(3.14159),3)*exp(-pow(x/1.195,2))",-3,3);
  fun3->SetLineColor(4);
  fun3->Draw("Same");


}
