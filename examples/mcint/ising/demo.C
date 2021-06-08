{
  gROOT->LoadMacro("randwalk_ising.C+");
  gROOT->LoadMacro("coolheat.C");

  coolheat c(5);

  c.setT(0.0);
  c.rampT(5.0, 200);

}
