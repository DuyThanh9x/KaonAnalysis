double skewt(double *x, double *par) {
    double mu    = par[1];
    double sigma = par[2];
    double alpha = par[3];
    double nu    = par[4];
    double y = (x[0] - mu) / sigma;
    double pdf_t = TMath::Student(y, nu);
    double sk = alpha * y * TMath::Sqrt((nu+1)/(nu + y*y));
    double cdf_t = TMath::StudentI(sk, nu+1.);
    return par[0]*(2.0/sigma) * pdf_t * cdf_t;
}

void msptfit(TString filename)
{
        TFile *file = TFile::Open(filename);

        vector<double> mean1G, dmean1G, momentum1G, dmomentum1G;
        mean1G.clear(); dmean1G.clear();momentum1G.clear(); dmomentum1G.clear();

        vector<double> mean1G5, dmean1G5, momentum1G5, dmomentum1G5;
        mean1G5.clear(); dmean1G5.clear();momentum1G5.clear(); dmomentum1G5.clear();

        vector<double> mean2G, dmean2G, momentum2G, dmomentum2G;
        mean2G.clear(); dmean2G.clear();momentum2G.clear(); dmomentum2G.clear();

        gStyle->SetOptStat(11);
        gStyle->SetOptFit(111);
        TH1D *mass = new TH1D();

        TF1 *z=new TF1("z",skewt,-9,5,5);
        z->SetNpx(900);
	z->SetParName(1,"Mean1");
	z->SetParName(2,"Sigma1");
	z->SetParName(3,"Skewness1");
	z->SetParName(4,"NDF1");
	for (int i = 5; i < 14; i++) {
                //auto c = new TCanvas();
                //c->SetLogy();
                mass = (TH1D*)file->Get(Form("mass1GeV%d",i));
                mass->SetLineColor(14);z->SetParameters(15,0.9,0.05,1., 2);
                z->SetParLimits(1,0.8,1.);
                z->SetParLimits(2,0.018,0.15);
                z->SetParLimits(3,0,2.9);
		if ((i ==5) || (i ==6)) {
                       if (i == 5) 
                               z->FixParameter(4,5);
                       //mass->Fit(z,"SEM0","",0.8,1.16); 
			     else  z->ReleaseParameter(4);
                       mass->Fit(z,"SEMQ0","",0.8,1.16);
                }
                else if (i == 13) {
                               z->FixParameter(4,18);
                               mass->Fit(z,"SEMQ0","",0.78,1.1);
                       }
                       else {
			       z->ReleaseParameter(4);
			       mass->Fit(z,"SEMQ0","",0.76,1.1);}

                mean1G.push_back(z->GetParameter(1));
                dmean1G.push_back(z->GetParError(1));

                double mment=(0.02*i+0.02*(i-1))/2;
                momentum1G.push_back(mment);
                dmomentum1G.push_back(0);
	}

	auto c = new TCanvas();
        TGraphErrors *mv1G = new TGraphErrors (mean1G.size(),momentum1G.data(),mean1G.data(),dmomentum1G.data(),dmean1G.data());
        mv1G->SetMarkerStyle(20);
        mv1G->SetMarkerSize(.7);
        mv1G->SetTitle("(mass)^{2} vs transverse momentum in momentum range [0.9, 1] GeV;GeV/Q;(GeV/Q)^{2}");
        mv1G->Draw("AP");

        for (int i = 2; i < 15; i++) {
                //auto c = new TCanvas();
                //c->SetLogy();
                mass = (TH1D*)file->Get(Form("mass15GeV%d",i));z->SetParameters(15,0.9,0.05,1., 2);
                z->SetParLimits(1,0.8,1.);
                z->SetParLimits(2,0.018,0.15);
                z->SetParLimits(3,0,2.9);
		if (i==6) {
                        z->FixParameter(4,15);
                        mass->Fit(z,"SEMQ0","",0.7,1.16);
                }
                else {
                    z->ReleaseParameter(4);
                    mass->Fit(z,"SEMQ0","",0.7,1.16);
                }

                mean1G5.push_back(z->GetParameter(1));
                dmean1G5.push_back(z->GetParError(1));

                double mment=(0.02*i+0.02*(i-1))/2;
                momentum1G5.push_back(mment);
                dmomentum1G5.push_back(0);
	}

	auto cv = new TCanvas();
        TGraphErrors *mv1G5 = new TGraphErrors (mean1G5.size(),momentum1G5.data(),mean1G5.data(),dmomentum1G5.data(),dmean1G5.data());
        mv1G5->SetMarkerStyle(20);
        mv1G5->SetMarkerSize(.7);
        mv1G5->SetTitle("(mass)^{2} vs transverse momentum in momentum range [1.4, 1.5] GeV;GeV/Q;(GeV/Q)^{2}");
        mv1G5->Draw("AP");

	for (int i = 3; i < 25; i++) {
                //auto c = new TCanvas();
                //c->SetLogy();
                mass = (TH1D*)file->Get(Form("mass2GeV%d",i));
                mass->SetLineColor(14);
		z->SetParameters(15,0.9,0.05,1., 2);
                z->SetParLimits(1,0.8,1.);
                z->SetParLimits(2,0.018,0.15);
                z->SetParLimits(3,0,2.9);
		mass->Fit(z,"SEMQ0","",0.7,1.16);
                mean2G.push_back(z->GetParameter(1));
                dmean2G.push_back(z->GetParError(1));
                double mment=(0.02*i+0.02*(i-1))/2+0.05;
                momentum2G.push_back(mment);
                dmomentum2G.push_back(0);
	}

	auto cvv = new TCanvas();
        TGraphErrors *mv1zz = new TGraphErrors (mean2G.size(),momentum2G.data(),mean2G.data(),dmomentum2G.data(),dmean2G.data());
        mv1zz->SetMarkerStyle(20);
        mv1zz->SetMarkerSize(.7);
        mv1zz->SetTitle("(mass)^{2} vs transverse momentum in momentum range [1.9, 2] GeV;GeV/Q;(GeV/Q)^{2}");
        mv1zz->Draw("AP");
}
