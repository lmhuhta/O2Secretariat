//#include "include/const_2pc.h"
#include <TFile.h>
#include <TH2D.h>
#include <TH1D.h>
#include <TString.h>

void Plot2DHist(){
	TFile* fin = new TFile("AnalysisResults.root", "read");
	TDirectoryFile* tdf = (TDirectoryFile*)fin->Get("lauras-code"); //twoparcorcombexample lauras-code
	TH2F* hDeltaPhiDeltaEta = (TH2F*) tdf->Get("hDeltaPhiDeltaEta")->Clone();
	TH1F* hTrigger = (TH1F*)tdf->Get("hTrigg")->Clone();

	TCanvas* c = new TCanvas("c","c",600,600);

	hDeltaPhiDeltaEta->GetZaxis()->SetTitle()

	hDeltaPhiDeltaEta->Draw("surf1");

	c-> Update();
	c->Show();
}