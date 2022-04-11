#include "TFile.h"
#include "TH2.h"
#include "TH2Poly.h"
#include "TGraphAsymmErrors.h"
#include "TRandom3.h"

#include <iostream>
#include <vector>
#include <algorithm>
#include <numeric>
#include <map>


float mjj_bins[4] = {250,350,700,13500};
float tag_bins[3] = {-0.5, 0.5, 1.5};
float D0minus_bins[4] = {0, 0.7, 0.9, 1};
TH1F *mjj_bins_histo = new TH1F("mjj_bins_histo","",3,mjj_bins);
TH1F *tag_bins_histo = new TH1F("tag_bins_histo","",2,tag_bins);
TH1F *D0minus_bins_histo = new TH1F("dnn_bins_histo","",3,D0minus_bins);

int mela_catIndex(float mjj, float dipho_bdt, float p_vbf, float p_ggH, float DOminus) {
  int mjj_bin = mjj_bins_histo->FindBin(mjj);
  if (mjj_bin==0) return -1; // cut these events 
  int tag_bin = -1;
  if (mjj_bin==1) { // SM cat for low mjj, high ptHjj
    if (dipho_bdt > 0.826 && p_vbf > 0.200 && p_ggH < 0.928) tag_bin = 0;
    else if (dipho_bdt > 0.650 && p_vbf > 0.137 && p_ggH < 0.876) tag_bin = 1;
  }
  else if (mjj_bin==2) { // SM cat
    if (dipho_bdt > 0.800 && p_vbf > 0.3 && p_ggH < 0.565) tag_bin = 0;
    else if (dipho_bdt > 0.779 && p_vbf > 0.2 && p_ggH < 0.726) tag_bin = 1;
    // if (dipho_bdt > 0.800 && p_vbf > 0.379 && p_ggH < 0.565) tag_bin = 0;
    // else if (dipho_bdt > 0.779 && p_vbf > 0.279 && p_ggH < 0.726) tag_bin = 1;
  } else { // SM cat
    if (dipho_bdt > 0.790 && p_vbf > 0.310 && p_ggH < 0.583) tag_bin = 0;
    else if (dipho_bdt > 0.606 && p_vbf > 0.160 && p_ggH < 0.919) tag_bin = 1;
    // if (dipho_bdt > 0.790 && p_vbf > 0.410 && p_ggH < 0.583) tag_bin = 0;
    // else if (dipho_bdt > 0.606 && p_vbf > 0.260 && p_ggH < 0.919) tag_bin = 1;
  }
  if (tag_bin<0) return -1; // cut these events
  int mela_bin = D0minus_bins_histo->FindBin( DOminus );
  // std::cout << mjj << "\t" << dipho_bdt << "\t"  << p_vbf << "\t" << p_ggH << "\t" << DOminus << "\t===> " << tag_bin << std::endl; 
  // std::cout << "mjj_bin = " << mjj_bin << std::endl;
  // std::cout << "tag_bin = " << tag_bin << std::endl;
  // std::cout << "mela_bin = " << mela_bin << std::endl;
  int cat = D0minus_bins_histo->GetNbinsX() * tag_bins_histo->GetNbinsX() * (mjj_bin-1) + 
    tag_bins_histo->GetNbinsX() * (mela_bin - 1) + 
    tag_bin + 1; 
  // std::cout << "category = " << cat << std::endl;
  return cat;
}

float dnn_pbkg_bins[3] = {0., 0.3, 1.0};
float dnn_pbsm_bins[4] = {0., 0.3, 0.7, 1.0};
float D0minus2_bins[3] = {0., 0.2, 1.0};
TH1F *dnn_pbkg_bins_histo = new TH1F("dnn_pbkg_bins_histo","",2,dnn_pbkg_bins);
TH1F *dnn_pbsm_bins_histo = new TH1F("dnn_pbsm_bins_histo","",3,dnn_pbsm_bins);
TH1F *D0minus2_bins_histo = new TH1F("D0minus2_bins_histo","",2,D0minus2_bins);

int dnn_catIndex(float p_bkg, float p_bsm, float DOminus) {
  int dnn_pbkg_bin = dnn_pbkg_bins_histo->FindBin(p_bkg);
  if (dnn_pbkg_bin==0 || dnn_pbkg_bin==dnn_pbkg_bins_histo->GetNbinsX()+1) return -1;
  int dnn_pbsm_bin = dnn_pbsm_bins_histo->FindBin(p_bsm);
  if (dnn_pbsm_bin==0 || dnn_pbsm_bin==dnn_pbsm_bins_histo->GetNbinsX()+1) return -1;
  int mela_bin = D0minus2_bins_histo->FindBin( DOminus );
  // std::cout << "dnn_pbkg_bin = " << dnn_pbkg_bin << std::endl;
  // std::cout << "dnn_pbsm_bin = " << dnn_pbsm_bin << std::endl;
  // std::cout << "mela_bin = " << mela_bin << std::endl;

  // for bkg-dominated events, just divide in D0 minus to separate ggH from pp->gg better
  int cat=-1;
  if (dnn_pbkg_bin==dnn_pbkg_bins_histo->GetNbinsX()) cat = D0minus2_bins_histo->GetNbinsX() * dnn_pbsm_bins_histo->GetNbinsX() + mela_bin;
  else cat = D0minus2_bins_histo->GetNbinsX() * dnn_pbsm_bins_histo->GetNbinsX() * (dnn_pbkg_bin-1) +
         dnn_pbsm_bins_histo->GetNbinsX() * (mela_bin-1) + 
         dnn_pbsm_bin;
  
  //  std::cout << "category = " << cat << std::endl;
  return cat;
}
