{
  gROOT->Reset();
  //gROOT->SetStyle("Plain");
  gStyle->SetOptStat(0);
  gStyle->SetOptFit(0);
/////////////////////////////////////////////////////////////////
  TFile *f00  = new TFile("W+Jets.root");
  TFile *f01  = new TFile("W+Jets.root");
  TFile *f02  = new TFile("DYJetsToLL.root");
  TFile *f03  = new TFile("tbar{t}.root");
  TFile *f04  = new TFile("SingleTop.root");
  TFile *f06  = new TFile("VV.root");

  TFile *f90 = new TFile("QCD.root");
  TFile *f91 = new TFile("W+Jets.root");
  TFile *f92 = new TFile("DYJetsToLL.root");
  TFile *f93 = new TFile("tbar{t}.root");
  TFile *f94 = new TFile("SingleTop.root");
  TFile *f95 = new TFile("VV.root");
///////////////////////////////////////////////////////////////////
  TH1F *h_wj;
  TH1F *h_dy;
  TH1F *h_tt;
  TH1F *h_st;
  TH1F *h_vv;
  TH1F *h_back;
  TH1F *h_qc;
  TH1F *h_wj9;
  TH1F *h_dy9;
  TH1F *h_tt9;
  TH1F *h_st9;
  TH1F *h_vv9;

  h_wj = (TH1F *)f01->Get("NRecoBJet/DiTauReconstructableMassDeltaPt");
    h_wj->SetLineColor(1);
    h_wj->SetFillColor(kAzure+6);
    h_wj->Scale(1);
  h_dy = (TH1F *)f02->Get("NRecoBJet/DiTauReconstructableMassDeltaPt");
    h_dy->SetLineColor(1);
    h_dy->SetFillColor(kOrange-4);
    h_dy->Scale(1);
  h_tt = (TH1F *)f03->Get("NRecoBJet/DiTauReconstructableMassDeltaPt");
    h_tt->SetLineColor(1);
    h_tt->SetFillColor(kBlue-8);
    h_tt->Scale(1);
  h_st = (TH1F *)f04->Get("NRecoBJet/DiTauReconstructableMassDeltaPt");
    h_st->SetLineColor(1);
    h_st->SetFillColor(8);
  h_vv = (TH1F *)f06->Get("NRecoBJet/DiTauReconstructableMassDeltaPt");
    h_vv->SetLineColor(1);
    h_vv->SetFillColor(46);

  h_qc = (TH1F *)f90->Get("NRecoBJet/DiTauReconstructableMassDeltaPt");
  h_qc->SetLineColor(1);
  h_qc->SetFillColor(kMagenta-10);
  //h_qc->Add(h_wj9,-1);
  //h_qc->Add(h_dy9,-1);
  //h_qc->Add(h_tt9,-1);
  //h_qc->Add(h_st9,-1);
  //h_qc->Add(h_vv9,-1);
  //h_qc->Scale(0.0244); //0.133, 0.183

  h_back = (TH1F *)f00->Get("NRecoBJet/DiTauReconstructableMassDeltaPt");
  h_back->Add(h_dy);
  h_back->Add(h_tt);
  h_back->Add(h_st);
  h_back->Add(h_vv);
  h_back->Add(h_qc);

//Double_t nbins = 10;
Double_t nbins = h_back->GetXaxis()->GetNbins();
double edy, ewj, ett, est, evv, eback, edata, enmc, eqc;
double intdy = h_dy->IntegralAndError(1,nbins,edy);
double intwj = h_wj->IntegralAndError(1,nbins,ewj);
double inttt = h_tt->IntegralAndError(1,nbins,ett);
double intst = h_st->IntegralAndError(1,nbins,est);
double intvv = h_vv->IntegralAndError(1,nbins,evv);
double intqc = h_qc->IntegralAndError(1,nbins,eqc);
double intback = h_back->IntegralAndError(1,nbins,eback);

cout<<"-----------  beamer table ------------------"<<endl;
cout<<"\\begin{tabular}{lcr}"<<endl;
cout<<"\\hline\\hline"<<endl;
cout<<"Process   & Events          & \\%   \\\\ "<<endl; 
cout<<"\\hline"<<endl;
cout<<"W+jets    & "<<Form("%.1f",intwj)<<"$\\pm$"<<Form("%.1f",ewj)<<" & "<<Form("%.0f",100*intwj/intback)<<"\\% \\\\ "<<endl;
cout<<"Drell-Yan & "<<Form("%.1f",intdy)<<"$\\pm$"<<Form("%.1f",edy)<<" & "<<Form("%.0f",100*intdy/intback)<<"\\% \\\\ "<<endl;
cout<<"ttbar     & "<<Form("%.1f",inttt)<<"$\\pm$"<<Form("%.1f",ett)<<" & "<<Form("%.0f",100*inttt/intback)<<"\\% \\\\ "<<endl;
cout<<"SingleTop & "<<Form("%.1f",intst)<<"$\\pm$"<<Form("%.1f",est)<<" & "<<Form("%.0f",100*intst/intback)<<"\\% \\\\ "<<endl;
cout<<"Diboson   & "<<Form("%.1f",intvv)<<"$\\pm$"<<Form("%.1f",evv)<<" & "<<Form("%.0f",100*intvv/intback)<<"\\% \\\\ "<<endl;    
cout<<"QCD       & "<<Form("%.1f",intqc)<<"$\\pm$"<<Form("%.1f",eqc)<<" & "<<Form("%.0f",100*intqc/intback)<<"\\% \\\\ "<<endl;
cout<<"\\hline"<<endl;
cout<<"TotalBack.& "<<Form("%.1f",intback)<<"$\\pm$"<<Form("%.1f",eback)<<" & "<<Form("%.0f",100*intback/intback)<<"\\% \\\\ "<<endl;
cout<<"\\hline\\hline"<<endl;
cout<<"\\end{tabular}"<<endl;
cout<<"--------------------------------------------"<<endl;

double xbins[14] = {0,100,200,300,400,500,600,700,800,900,1000,1100,1200,2000};
    h_wj->Rebin(13,"hb_wj",xbins);
    h_dy->Rebin(13,"hb_dy",xbins);
    h_tt->Rebin(13,"hb_tt",xbins);
    h_st->Rebin(13,"hb_st",xbins);
    h_vv->Rebin(13,"hb_vv",xbins);
    h_qc->Rebin(13,"hb_qc",xbins);
    h_back->Rebin(13,"hb_back",xbins);

/*
    h_wj->Rebin(10);
    h_dy->Rebin(10);
    h_tt->Rebin(10);
    h_st->Rebin(10);
    h_vv->Rebin(10);
    h_data->Rebin(10);
    h_back->Rebin(10);
    h_ratio->Rebin(10);
*/
  THStack *hs_n = new THStack("hs_n","");
    hs_n->Add(hb_st);
    hs_n->Add(hb_vv);
    hs_n->Add(hb_tt);
    hs_n->Add(hb_wj);
    hs_n->Add(hb_qc);
    hs_n->Add(hb_dy);
    

  TCanvas *c0=new TCanvas("c0","rm_nn",500,500);
  //c0->Divide(2,1);

   c0->cd();
        gPad->SetLogy();
        gPad->SetTickx();
        gPad->SetTicky();
    	gPad->SetLeftMargin(0.13);
    	gPad->SetRightMargin(0.05);
    	gPad->SetBottomMargin(0.12);
    	gPad->SetTopMargin(0.08);
   hs_n->Draw(); 
   hs_n->SetMaximum(750);
   hs_n->SetMinimum(0.25);
   hs_n->GetXaxis()->SetRangeUser(0,2000);
   hs_n->Draw("HIST");
   hb_back->Draw("E2same");
   hb_back->SetLineWidth(1);
   hb_back->SetFillColor(1);
   hb_back->SetFillStyle(3002);
  hs_n->GetXaxis()->SetLabelSize(0.035);
  hs_n->GetYaxis()->SetLabelSize(0.04);
  hs_n->GetXaxis()->SetTitleSize(0.05);
  hs_n->GetYaxis()->SetTitleSize(0.05);
  hs_n->GetYaxis()->SetTitleOffset(1.3);
  hs_n->GetXaxis()->SetTitleOffset(1.1);
   hs_n->SetTitle(" ; m_{rec}(#tau, #tau, #Delta p_{T}) [GeV]  ; Events");
          leg1=new TLegend(0.6,0.6,0.93,0.89);
          //leg1->SetHeader("CMS Preliminary","C");
	  leg1->SetBorderSize(0);
          leg1->AddEntry(h_dy,"Drell-Yan","fp");
          leg1->AddEntry(h_wj,"W+jets","fp");
          leg1->AddEntry(h_tt,"t#bar{t}","fp");
          leg1->AddEntry(h_vv,"Diboson","fp");
          leg1->AddEntry(hb_qc,"QCD","fp");
          leg1->AddEntry(h_st,"Single Top","fp");
          leg1->AddEntry(hb_back,"Stat. Uncert.","fp");
          leg1->Draw();

    TPaveText *pt = new TPaveText(0.11,0.92,0.99,0.99,"NBNDC");
      pt->AddText("Signal Region                 59.7 fb^{-1} (13 TeV, 2018)");
      pt->SetTextFont(42);
      pt->SetTextAlign(32);
      pt->SetFillStyle(0);
      pt->SetBorderSize(0);
      pt->Draw();                                                                                             
    TPaveText *pt2 = new TPaveText(0.39,0.83,0.51,0.9,"NBNDC");
      pt2->AddText("CMS ");
      pt2->SetTextAlign(12);
      pt2->SetFillStyle(0);
      pt2->SetBorderSize(0);
      pt2->Draw();
    TPaveText *pt3 = new TPaveText(0.38,0.8,0.51,0.86,"NBNDC");
      pt3->AddText("Work in Progress");
      pt3->SetTextAlign(12);
      pt3->SetTextFont(52);
      pt3->SetFillStyle(0);
      pt3->SetBorderSize(0);
      pt3->Draw();

c0->SaveAs("tautau_SR18_092320.pdf");
TFile *f = new TFile("SR18_092320.root","RECREATE");
h_qc->Write("SR18_092320");
f->Write();
f->Close();
}
