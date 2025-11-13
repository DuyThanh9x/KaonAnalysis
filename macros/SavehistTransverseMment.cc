void SavehistTransverseMment(TString name)
{
        TFile *file = TFile::Open(name);
        TTree *data;
        file->GetObject("BmnData",data);

	TFile *out = new TFile("25/histmasspt.root","recreate");
	
	data->Draw("sqrt(TrackPx*TrackPx+TrackPy*TrackPy):TrackTof400Mq2>>ms1GeV(250,-.2,1.5,14,0.0,.28)","TrackTof400Mq2<1.6 && TrackPq<=1. && TrackPq>=.9 && TrackTof400Mq2>-.3","goff");
	TH2D *mass1 = (TH2D*)gDirectory->Get("ms1GeV");

	data->Draw("sqrt(TrackPx*TrackPx+TrackPy*TrackPy):TrackTof400Mq2>>ms1GeV5(250,-.2,1.5,16,0.0,.32)","TrackTof400Mq2<1.6 && TrackPq>=1.4 && TrackPq<=1.5 && TrackTof400Mq2>-.3","goff");
        TH2D *mass1GeV5 = (TH2D*)gDirectory->Get("ms1GeV5");

	data->Draw("sqrt(TrackPx*TrackPx+TrackPy*TrackPy):TrackTof400Mq2>>ms2GeV(250,-.2,1.5,25,0.05,.55)","TrackTof400Mq2<1.6 && TrackPq>=1.9 && TrackPq<=2. && TrackTof400Mq2>-.3","goff");
        TH2D *mass2 = (TH2D*)gDirectory->Get("ms2GeV");
        
	for (int i = 1; i < 15; i++) {
		mass1->ProjectionX(Form("mass1GeV%d",i),i,i);
                TH1D *mass = (TH1D*)gDirectory->Get(Form("mass1GeV%d",i));
                mass->SetTitle(Form("Mq2 distribution at TOF400 with pt range [%.2f, %.2f] GeV in momentum range [0.9, 1.0] GeV;GeV^{2};#Candidates",.02*(i-1),.02*i));
		mass->Write();
	}

	for (int i = 1; i < 17; i++) {
                mass1GeV5->ProjectionX(Form("mass15GeV%d",i),i,i);
                TH1D *mass15 = (TH1D*)gDirectory->Get(Form("mass15GeV%d",i));
                mass15->SetTitle(Form("Mq2 distribution at TOF400 with pt range [%.2f, %.2f] GeV in momentum range [1.4, 1.5] GeV;GeV^{2};#Candidates",.02*(i-1),.02*i));
                mass15->Write();
        }

	for (int i = 1; i < 26 ; i++) {
                mass2->ProjectionX(Form("mass2GeV%d",i),i,i);
                TH1D *mass2s = (TH1D*)gDirectory->Get(Form("mass2GeV%d",i));
                mass2s->SetTitle(Form("Mq2 distribution at TOF400 with pt range [%.2f, %.2f] GeV in momentum range [1.9, 2.0] GeV;GeV^{2};#Candidates",.02*(i-1)+.05,.02*i+.05));
                mass2s->Write();
        }
	
	out->Close();
	//auto c = new TCanvas();
	//mass1->Draw();
	//auto c1 = new TCanvas();
        //mass1GeV5->Draw();
	//auto c2 = new TCanvas();
	//mass2->Draw();
}
