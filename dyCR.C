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
  TFile *f900 = new TFile("Data.root");
  TFile *f901 = new TFile("Data.root");
///////////////////////////////////////////////////////////////////
  TH1F *h_wj;
  TH1F *h_dy;
  TH1F *h_tt;
  TH1F *h_st;
  TH1F *h_vv;
  TH1F *h_data;
  TH1F *h_ratio;
  TH1F *h_back;

  h_wj = (TH1F *)f01->Get("NRecoBJet/DiTauReconstructableMassDeltaPt");
    h_wj->SetLineColor(1);
    h_wj->SetFillColor(kAzure+6);
    h_wj->Scale(1);
  h_dy = (TH1F *)f02->Get("NRecoBJet/DiTauReconstructableMassDeltaPt");
    h_dy->SetLineColor(1);
    h_dy->SetFillColor(kOrange-4);
    //h_dy->Scale(0.98);
  h_tt = (TH1F *)f03->Get("NRecoBJet/DiTauReconstructableMassDeltaPt");
    h_tt->SetLineColor(1);
    h_tt->SetFillColor(kBlue-8);
  h_st = (TH1F *)f04->Get("NRecoBJet/DiTauReconstructableMassDeltaPt");
    h_st->SetLineColor(1);
    h_st->SetFillColor(8);
  h_vv = (TH1F *)f06->Get("NRecoBJet/DiTauReconstructableMassDeltaPt");
    h_vv->SetLineColor(1);
    h_vv->SetFillColor(46);
  h_ratio = (TH1F *)f901->Get("NRecoBJet/DiTauReconstructableMassDeltaPt");
  h_data = (TH1F *)f900->Get("NRecoBJet/DiTauReconstructableMassDeltaPt");
    //h_data->Scale(8.9);
    h_data->SetLineColor(1);
    h_data->SetLineStyle(1);
    h_data->SetMarkerStyle(20);
    h_data->SetLineWidth(1);


    h_dy->Sumw2();
    h_tt->Sumw2();
    h_st->Sumw2();
    h_vv->Sumw2();
    h_wj->Sumw2();


  h_back = (TH1F *)f00->Get("NRecoBJet/DiTauReconstructableMassDeltaPt");
  h_back->Sumw2();
  h_back->Add(h_dy);
  h_back->Add(h_tt);
  h_back->Add(h_st);
  h_back->Add(h_vv);

//Double_t nbins = 8;
Double_t nbins = h_back->GetXaxis()->GetNbins();
TH1F *h_nmc=(TH1F*)h_data->Clone();
h_nmc->SetName("h_nmc");
h_nmc->Add(h_back,-1);
h_nmc->Add(h_dy,1);
double edy, ewj, ett, est, evv, eback, edata, enmc;
double intdy = h_dy->IntegralAndError(1,nbins,edy);
double intwj = h_wj->IntegralAndError(1,nbins,ewj);
double inttt = h_tt->IntegralAndError(1,nbins,ett);
double intst = h_st->IntegralAndError(1,nbins,est);
double intvv = h_vv->IntegralAndError(1,nbins,evv);
double intback = h_back->IntegralAndError(1,nbins,eback);
double intdat = h_data->IntegralAndError(1,nbins,edata);
double intnmc = h_nmc->IntegralAndError(1,nbins,enmc);
double sf = intnmc/intdy;
double esf = sqrt( pow((enmc/intdy),2)
	          +pow((intnmc*edy/(intdy*intdy)),2));

//cout<<" "<<Form("%.2f",intnmc)<<" "<<endl;
cout<<"-----------  beamer table ------------------"<<endl;
cout<<"\\begin{tabular}{lc}"<<endl;
cout<<"\\hline\\hline"<<endl;
cout<<"Process   & Events             \\\\ "<<endl; 
cout<<"\\hline"<<endl;
cout<<"W+jets    & "<<Form("%.1f",intwj)<<"$\\pm$"<<Form("%.1f",ewj)<<" \\\\ "<<endl;
cout<<"Drell-Yan & "<<Form("%.1f",intdy)<<"$\\pm$"<<Form("%.1f",edy)<<" \\\\ "<<endl;
cout<<"ttbar     & "<<Form("%.1f",inttt)<<"$\\pm$"<<Form("%.1f",ett)<<" \\\\ "<<endl;
cout<<"SingleTop & "<<Form("%.1f",intst)<<"$\\pm$"<<Form("%.1f",est)<<" \\\\ "<<endl;
cout<<"Diboson   & "<<Form("%.1f",intvv)<<"$\\pm$"<<Form("%.1f",evv)<<" \\\\ "<<endl;    
cout<<"\\hline"<<endl;
cout<<"TotalBack.& "<<Form("%.1f",intback)<<"$\\pm$"<<Form("%.1f",eback)<<" \\\\ "<<endl;
cout<<"Data      & "<<intdat<<"       \\\\ "<<endl;
cout<<"\\hline"<<endl;
cout<<"Purity & "<<Form("%.0f",100*intdy/intback)<<"\\%   \\\\ "<<endl; 
cout<<"\\hline"<<endl;
cout<<"Scale Factor & "<<Form("%.2f",sf)<<" $\\pm$"<<Form("%.2f",esf)<<"   \\\\ "<<endl; 
cout<<"\\hline\\hline"<<endl;
cout<<"\\end{tabular}"<<endl;
cout<<"--------------------------------------------"<<endl;

 double xbins_massdeltapt[8] = {200,300,400,500,600,700,800,1200};
 double xbins_mass[13] = {40,45,50,55,60,65,70,75,80,85,90,95,100};//200,300,400,500,600,700,800,1200};
 double xbins_taupt[11] = {70,90,110,130,150,170,190,210,230,250,270};
 double xbins_jet[8] = {0,1,2,3,4,5,6,7};
 double xbins_jetpt[18] = {30,50,70,90,110,130,150,170,190,210,230,250,270,290,310,330,350,370};
 double xbins_ht[13] = {0,100,200,300,400,500,600,700,800,900,1000,1200,2500};
 h_wj->Rebin(7,"hb_wj",xbins_massdeltapt);
    h_dy->Rebin(7,"hb_dy",xbins_massdeltapt);
    h_tt->Rebin(7,"hb_tt",xbins_massdeltapt);
    h_st->Rebin(7,"hb_st",xbins_massdeltapt);
    h_vv->Rebin(7,"hb_vv",xbins_massdeltapt);
    h_data->Rebin(7,"hb_data",xbins_massdeltapt);
    h_back->Rebin(7,"hb_back",xbins_massdeltapt);
    h_ratio->Rebin(7,"hb_ratio",xbins_massdeltapt);
 // h_wj->Rebin(1);
 // h_dy->Rebin(1);
 // h_tt->Rebin(1);
 // h_st->Rebin(1);
 //    h_vv->Rebin(1);
 //    h_data->Rebin(1);
 //    h_back->Rebin(1);
 //    h_ratio->Rebin(1);

  THStack *hs_n = new THStack("hs_n","");
    hs_n->Add(hb_st);    
    hs_n->Add(hb_tt);
    hs_n->Add(hb_wj);
    hs_n->Add(hb_vv);
    hs_n->Add(hb_dy);

  TCanvas *c0=new TCanvas("c0","rm_nn",500,500);
  //c0->Divide(2,1);

   c0->cd();
TPad *pad1 = new TPad("pad1", "pad1", 0, 0, 1, 1);
  pad1->SetBottomMargin(0.3); // Upper and lower plot are joined
  pad1->Draw();             // Draw the upper pad: pad1
  pad1->cd();               // pad1 becomes the current pad
  //gPad->SetLogy();
        gPad->SetTickx();
        gPad->SetTicky();
    	gPad->SetLeftMargin(0.15);
    	gPad->SetRightMargin(0.05);
    	//gPad->SetBottomMargin(0.15);
    	gPad->SetTopMargin(0.06);
   hs_n->Draw("same"); 
   hs_n->SetMaximum(5000000);
   hs_n->SetMinimum(0.001);
   //hs_n->GetXaxis()->SetRangeUser(70,110);
   hs_n->Draw("HIST");
   hb_data->Draw("same");
   hb_back->Draw("E2same");
   hb_back->SetLineWidth(1);
   hb_back->SetFillColor(1);
   hb_back->SetFillStyle(3002);
   hs_n->SetTitle(" ;  ; Events");
   hs_n->GetYaxis()->SetTitleSize(0.04);
   //hs_n->SetTitle("CMS Preliminary                 35.9 fb^{-1} (13 TeV, 2016) ;  ; Events");
          leg1=new TLegend(0.67,0.63,0.93,0.91);
          //leg1->SetHeader("CMS Preliminary","C");
	  leg1->SetBorderSize(0);
          leg1->AddEntry(hb_data,"Data","lp");
          leg1->AddEntry(hb_dy,"Drell-Yan","fp");
          leg1->AddEntry(hb_wj,"W+jets","fp");
          leg1->AddEntry(hb_tt,"t#bar{t}","fp");
          leg1->AddEntry(hb_vv,"Diboson","fp");
          //leg1->AddEntry(h_qc,"data-driven QCD","fp");
          leg1->AddEntry(hb_st,"Single Top","fp");
          leg1->AddEntry(hb_back,"Stat. Uncert.","fp");
          leg1->Draw();

    TPaveText *pt = new TPaveText(0.11,0.94,0.99,0.99,"NBNDC");
      pt->AddText("QCD CR                 35.9 fb^{-1} (13 TeV, 2016)");
      pt->SetTextFont(42);
      pt->SetTextAlign(32);
      pt->SetFillStyle(0);
      pt->SetBorderSize(0);
      pt->Draw();                                                                                             
    TPaveText *pt2 = new TPaveText(0.19,0.85,0.31,0.92,"NBNDC");
      pt2->AddText("CMS ");
      pt2->SetTextAlign(12);
      pt2->SetFillStyle(0);
      pt2->SetBorderSize(0);
      pt2->Draw();
    TPaveText *pt3 = new TPaveText(0.18,0.82,0.31,0.88,"NBNDC");
      pt3->AddText("Preliminary");
      pt3->SetTextAlign(12);
      pt3->SetTextFont(52);
      pt3->SetFillStyle(0);
      pt3->SetBorderSize(0);
      pt3->Draw();

      c0->cd(1);
TPad *pad2 = new TPad("pad2", "pad2", 0, 0, 1, 0.25);
  pad2->SetTopMargin(0);
  pad2->SetBottomMargin(0.35);
  pad2->Draw();
  pad2->cd();       // pad2 becomes the current pad
        gPad->SetTickx();
        gPad->SetTicky();
	gPad->SetLeftMargin(0.15);
    	gPad->SetRightMargin(0.05);
	TLine *line = new TLine(200,1,1200,1);
	line->SetLineColor(2);

hb_ratio->Divide(hb_back);
//hb_ratio->Draw();
int Nbins =  hb_ratio->GetXaxis()->GetNbins();               
Double_t* mcX = new Double_t[Nbins];                      
Double_t* mcY = new Double_t[Nbins];                      
Double_t* mcErrorX = new Double_t[Nbins];                 
Double_t* mcErrorY = new Double_t[Nbins];

for(int bin=0; bin < Nbins; bin++) {       
  mcY[bin] = 1.0;                                              
	mcErrorY[bin] =  hb_ratio->GetBinError(bin+1);       
	mcX[bin] = hb_ratio->GetBinCenter(bin+1);                    
	mcErrorX[bin] = hb_ratio->GetBinWidth(bin+1) * 0.5;   
}
TGraphErrors *mcError = new TGraphErrors(h_ratio->GetXaxis()->GetNbins(),mcX,mcY,mcErrorX,mcErrorY);   
mcError->SetLineWidth(1);                                 
mcError->SetFillColor(1);                                 
mcError->SetFillStyle(3002);

hb_ratio->Draw();
mcError->Draw("E2same");
line->Draw("same");
hb_ratio->SetMaximum(2.2);
hb_ratio->SetMinimum(-0.1);
hb_ratio->SetLineWidth(1);
hb_ratio->SetFillColor(1);
hb_ratio->SetLineColor(1);
hb_ratio->SetMarkerStyle(20);
//hb_ratio->GetXaxis()->SetRangeUser(0,2500);
hb_ratio->SetTitle("   ; m_{rec}(#tau, #tau, #Delta p_{T}) [GeV] ; #frac{Data}{Background} ");
//hb_ratio->SetTitle("   ; m_{rec}(#tau, #tau) [GeV] ; #frac{Data}{Background} ");
//hb_ratio->SetTitle("   ; p_{T}(#tau) [GeV] ; #frac{Data}{Background} ");
//hb_ratio->SetTitle("   ; p_{T}(j) [GeV] ; #frac{Data}{Background} ");
//hb_ratio->SetTitle("   ; H_{T} [GeV]; #frac{Data}{Background} ");
hb_ratio->SetMarkerStyle(33);
hb_ratio->GetXaxis()->SetLabelSize(0.14);
hb_ratio->GetYaxis()->SetLabelSize(0.14);
hb_ratio->GetXaxis()->SetTitleSize(0.16);
hb_ratio->GetYaxis()->SetTitleSize(0.14);
hb_ratio->GetYaxis()->SetTitleOffset(0.4);
hb_ratio->GetXaxis()->SetTitleOffset(0.9);

c0->SaveAs("diTau_dyCR2018_ditaurecomassdeltapt_100720.pdf");
}
