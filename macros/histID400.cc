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

double back(double *x, double *par)
{
        return par[0]*TMath::Exp(x[0]/par[1])+par[2]+par[3]*x[0];
}

double total(double *x, double *par)
{
        return skewt(x, par)+skewt(x,&par[5])+skewt(x,&par[10])+back(x, &par[15]);
}

void tf1s(TF1 *s[], TF1 *&bg, TF1 *&sum)
{
        s[0] = new TF1("s1",skewt,-9,5,5);
        s[0]->SetNpx(900);
        s[0]->SetLineColor(kBlue);
        s[0]->SetLineStyle(9);

        s[1] = new TF1("s2",skewt,-9,5,5);
        s[1]->SetNpx(900);
        s[1]->SetLineColor(kRed);
        s[1]->SetLineStyle(9);

        s[2] = new TF1("s3",skewt,-9,5,5);
        s[2]->SetNpx(900);
        s[2]->SetLineColor(kOrange);
        s[2]->SetLineStyle(9);

        bg = new TF1("bg",back,-9,5,4);
        bg->SetNpx(900);
        bg->SetLineColor(kBlack);
        bg->SetLineStyle(9);

        sum = new TF1("sum",total,-9,5,19);
        sum->SetNpx(900);
        sum->SetLineColor(9);
        sum->SetParName(0,"Coeff1");
        sum->SetParName(1,"Mean1");
        sum->SetParName(2,"Sigma1");
        sum->SetParName(3,"Skewness1");
        sum->SetParName(4,"NDF1");
        sum->SetParName(5,"Coeff2");
        sum->SetParName(6,"Mean2");
        sum->SetParName(7,"Sigma2");
        sum->SetParName(8,"Skewness2");
        sum->SetParName(9,"NDF2");
        sum->SetParName(10,"Coeff3");
	sum->SetParName(11,"Mean3");
        sum->SetParName(12,"Sigma3");
        sum->SetParName(13,"Skewness3");
        sum->SetParName(14,"NDF3");
        sum->SetParName(15,"CoeffExp");
        sum->SetParName(16,"ParExp");
        sum->SetParName(17,"Constant");
        sum->SetParName(18,"Coeff1st");
}

void histID400(TString name)
{
        TFile *file = TFile::Open(name);
        FILE* outfile =fopen("particleID400.txt","w");
        char *ln;
        string r;
        ln = Form("||     Momentum range      ||         Kaon window          | Kaon backgrnd   | Kaon loss   | Kaon weight   ||          Pion window         | Pion backgrnd  |  Pion loss    | Pion weight   ||        Proton window        | Proton backgrnd | Proton loss | Proton weight\n");
        r = ln;
        fwrite(ln,1,r.size(),outfile);
	double par[18], range[2];
        vector<double> mean1, mean2, mean3, sigma1, sigma2, sigma3,skew, momentum;
        vector<double> ermean1, ermean2, ermean3, ersigma1, ersigma2, ersigma3,erskew, ermomentum;

        vector<double> mean11, mean21, mean31, sigma11, sigma21, sigma31;
        vector<double> dmean11, dmean21, dmean31, dsigma11, dsigma21, dsigma31;

        mean1.clear(); mean2.clear(); mean3.clear(); sigma1.clear(); sigma2.clear(); sigma3.clear();skew.clear();
        mean11.clear(); mean21.clear(); mean31.clear(); sigma11.clear(); sigma21.clear(); sigma31.clear();
        dmean11.clear(); dmean21.clear(); dmean31.clear(); dsigma11.clear(); dsigma21.clear(); dsigma31.clear();
        momentum.clear();
        ermean1.clear(); ermean2.clear(); ermean3.clear(); ersigma1.clear(); ersigma2.clear(); ersigma3.clear();erskew.clear();ermomentum.clear();

        TH1D *mass = new TH1D();
        double meanvl[3], smvl[3];
        TF1 *s[3], *bg, *totall;
        double kwind1,kwind2,piwind1,piwind2,prwind1,prwind2;
        gStyle->SetOptStat(11);
        gStyle->SetOptFit(111);
        ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(100000);
        for (int i = 1; i < 22; i++) {
		//zauto c = new TCanvas();
                //zc->SetLogy();
                mass = (TH1D*)file->Get(Form("mass%d",i));
                mass->SetLineColor(14);
                double mment=(.1*i+.1*(i-1))/2+.5;

		meanvl[0]=0.031-0.018*TMath::Power(mment,1);
                meanvl[1]=-1.93+2.168*TMath::Power(mment,-0.0092);
                meanvl[2]=36.05-35.16*TMath::Power(mment,0.000709);

		smvl[0]=0.005-0.005*mment+0.01*mment*mment;
		smvl[1]=0.0369-0.0366*mment+0.0195*mment*mment;
		smvl[2]=0.096-0.056*mment+0.019*mment*mment;

		mean11.push_back(meanvl[0]); mean21.push_back(meanvl[1]);mean31.push_back(meanvl[2]);
                sigma11.push_back(smvl[0]); sigma21.push_back(smvl[1]);sigma31.push_back(smvl[2]);

		double par1[]={29, meanvl[0], smvl[0], 0.5, 2, 15, meanvl[1], smvl[1], 0.5, 1, 39, meanvl[2], smvl[2], 0.9, 2, 9, -6, 5, -9};
                tf1s(s,bg,totall);
                totall->SetParameters(par1);

                range[0]=-0.034;range[1]=1.68;
                if (i==1) {
                totall->SetParLimits(2,0.0048,0.02);
                totall->SetParLimits(3,0,1);
                totall->SetParLimits(6,0.21,0.29);
                totall->FixParameter(7,0.015);
                totall->FixParameter(8,0.01);
		totall->FixParameter(9,12);
		totall->SetParLimits(14,1,99);
                totall->SetParLimits(16,-18,0);
		totall->SetParLimits(18,-28,0);
		}

		if (i==2) {
                totall->SetParLimits(2,0.004,0.009);
                totall->SetParLimits(4,1.,41);
                totall->FixParameter(3,0.01);
                totall->SetParLimits(6,0.21,0.26);
                totall->SetParLimits(7,0.01,0.05);
		totall->FixParameter(9,12);
		totall->SetParLimits(16,-18,0);
		totall->SetParLimits(18,-28,0);
		}

		if (i==3) {
                        totall->SetParLimits(2,0.006,0.01);
                        totall->SetParLimits(7,0.01,0.03);
                        totall->FixParameter(3,0.01);
                        totall->FixParameter(8,0.02);
			totall->SetParLimits(9,1,332);
                        totall->SetParLimits(16,-18,0);
                        totall->SetParLimits(18,-28,0);
		}

		if (i==4) {
                totall->SetParLimits(2,0.005,0.01);
		totall->FixParameter(3,0.01);
                totall->SetParLimits(8,0.,1.5);
                totall->SetParameter(9,4);
                totall->SetParLimits(9,2,9);
                totall->SetParLimits(16,-18,0);
		totall->SetParLimits(18,-28,0);
		}

		if (i==5) {
                        totall->SetParLimits(2,0.008,0.02);
                        totall->SetParLimits(3,0.05,0.1);
			totall->SetParLimits(8,0.,2.);
			totall->SetParameter(16,-0.5);
                        totall->SetParLimits(16,-1,0);
		}

		if ((i==6) || (i==7)) {
			range[0]=-0.04;
		}

		if (i==8) {
			totall->SetParLimits(8,0.,2.);
                        range[0]=-0.044;
		}

		if (i==9) {
			totall->SetParLimits(8,0.,2);
                        totall->SetParLimits(16,-1,0);
                        range[0]=-0.052;
		}

		if (i==10) {
			totall->SetParameter(3,0.15);
                        totall->SetParLimits(3,0.0,0.3);
			totall->SetParameter(9,8);
                        totall->SetParLimits(9,5,19);
                        totall->SetParameter(8,0.2);
                        totall->SetParLimits(8,0.11,0.3);
                        totall->SetParameter(15,85);
                        totall->SetParLimits(15,60,89);
                        totall->SetParameter(16,-0.5);
			range[0]=-0.068;
		}

		if (i==11) {
			totall->SetParameter(1,0.0001);
			totall->FixParameter(8,0.01);
			totall->SetParameter(16,-0.5);
                        totall->SetParLimits(16,-1,0);
                        range[0]=-0.08;
		}
		
		if (i==12) {
			totall->FixParameter(8,0.03);
			totall->SetParameter(9,4);
                        totall->SetParLimits(9,2,9);
                        totall->SetParameter(15,54.9);
                        totall->SetParLimits(15,22,59);
                        totall->SetParameter(16,-0.5);
                        totall->SetParLimits(16,-1,0);
			range[0]=-0.1;
		}

		if (i==13) {
			totall->FixParameter(8,0.02);
                        totall->FixParameter(9,3);
			totall->SetParameter(15,1.2);
                        totall->SetParLimits(15,0,2);
			totall->SetParLimits(16,-18,0);
			range[0]=-0.14;
		}

		if (i==14) {
                        totall->FixParameter(8,0.01);
			totall->SetParLimits(15,0,19);
			totall->SetParLimits(16,-228,0);
			range[0]=-0.14;
		}

		if (i==15) {
			totall->SetParLimits(3,0.01,1.2);
			totall->FixParameter(8,0.01);
			totall->SetParameter(15,0.5);
                        totall->SetParLimits(15,0,1);
			totall->SetParLimits(16,-318,0);
                        totall->SetParLimits(18,-28,0);
                        range[0]=-0.14;
		}

		if (i==16) {
                        totall->SetParLimits(2,0.03,2);
			totall->FixParameter(3,0.01);
			totall->FixParameter(8,0.08);
                        totall->SetParameter(9,3.9);
			totall->SetParLimits(15,5,19);
			totall->SetParLimits(16,-38,0);
			range[0]=-0.14;
		}

		if (i==17) {
                        totall->SetParLimits(2,0.03,2);
			totall->FixParameter(8,0.02);
			totall->SetParameter(15,2);
                        totall->SetParLimits(15,0,4);
			totall->SetParLimits(16,-918,0);
			range[0]=-0.18;
		}

		if (i==18) {
			totall->SetParLimits(2,0.01,0.2);
			totall->FixParameter(3,0.01);
			totall->SetParLimits(7,0.04,0.06);
			totall->FixParameter(8,0.02);
			totall->SetParameter(9,3.9);
                        totall->SetParLimits(9,3,8);
                        totall->SetParameter(15,2);
                        totall->SetParLimits(15,0,4);
			totall->SetParLimits(16,-68,0);
                        totall->SetParLimits(18,-28,0);
                        range[0]=-0.18;
		}

		if (i==19) {
			totall->FixParameter(3,0.0);
                        totall->FixParameter(9,9);
                        totall->FixParameter(8,0.02);
                        totall->SetParameter(15,0.5);
                        totall->SetParLimits(15,0,9);
                        totall->SetParLimits(16,-91,0);
			range[0]=-0.2;
		}

		if (i==20) {
			totall->SetParLimits(1,-0.02,-0.01);
			totall->SetParameter(3,0.1);
                        totall->SetParLimits(3,0.0,0.2);
			totall->SetParLimits(7,0.03,0.08);
			totall->FixParameter(8,0.0);
			totall->FixParameter(9,6.5);
			totall->FixParameter(15,0.02);
                        totall->SetParameter(16,-0.1);
			totall->SetParLimits(16,-1,0);
			range[0]=-0.19;
		}

		if (i==21) {
			totall->SetParLimits(2,0.054,0.068);
                        totall->FixParameter(3,0.02);
			totall->SetParLimits(7,0.02,0.08);
                        totall->FixParameter(8,0.01);
			totall->SetParameter(9,12);
                        totall->SetParLimits(9,2,24);
                        totall->SetParameter(15,0.5);
                        totall->SetParLimits(15,0,15);
                        totall->SetParLimits(16,-28,0);
                        totall->SetParLimits(18,-28,0);
                        range[0]=-0.19;
		}

		mass->Fit(totall,"SEMQ0","",range[0],range[1]);

                totall->GetParameters(par);
                s[0]->SetParameters(par);
                s[1]->SetParameters(&par[5]);
                s[2]->SetParameters(&par[10]);
                bg->SetParameters(&par[15]);

                mean1.push_back(totall->GetParameter(1));
                ermean1.push_back(totall->GetParError(1));
                sigma1.push_back(totall->GetParameter(2));
                ersigma1.push_back(totall->GetParError(2));
                mean2.push_back(totall->GetParameter(6));
		ermean2.push_back(totall->GetParError(6));
		sigma2.push_back(totall->GetParameter(7));
                        ersigma2.push_back(totall->GetParError(7));
		mean3.push_back(totall->GetParameter(11));
                        ermean3.push_back(totall->GetParError(11));
                        ersigma3.push_back(totall->GetParError(12));
			sigma3.push_back(totall->GetParameter(12));
		momentum.push_back(mment);
                ermomentum.push_back(0);

		if (i==1) {
                        kwind1=0.22; kwind2=0.3;
                        piwind1=-0.01;piwind2=0.06;
                        prwind1=0.82;prwind2=1.22;
                }

                if (i==2) {
                        kwind1=0.22; kwind2=0.3;
                        piwind1=-0.02;piwind2=0.08;
                        prwind1=0.8;prwind2=1.16;
                }

                if (i==3) {
                        kwind1=0.22; kwind2=0.29;
                        piwind1=-0.02;piwind2=0.06;
                        prwind1=0.8;prwind2=1.14;
                }

                if (i==4) {
                        kwind1=0.22; kwind2=0.3;
                        piwind1=-0.03;piwind2=0.06;
                        prwind1=0.78;prwind2=1.16;
                }

                if (i==5) {
                        kwind1=0.22; kwind2=0.3;
                        piwind1=-0.03;piwind2=0.06;
                        prwind1=0.75;prwind2=1.22;
		}

		if (i==6) {
                        kwind1=0.2; kwind2=0.32;
                        piwind1=-0.03;piwind2=0.06;
                        prwind1=0.72;prwind2=1.26;
		}

		if (i==7) {
                kwind1=0.2; kwind2=0.32;
                piwind1=-0.04;piwind2=0.08;
		prwind1=0.73;prwind2=1.2;}

		if (i==8) {
                        kwind1=0.2; kwind2=0.32;
                        piwind1=-0.05;piwind2=0.1;
                        prwind1=0.73;prwind2=1.2;
		}

		if (i==9) {
                kwind1=0.18; kwind2=0.32;
                piwind1=-0.06;piwind2=0.1;
		prwind1=0.73;prwind2=1.18;
		}

		if (i==10) {
                kwind1=0.16; kwind2=0.32;
                piwind1=-0.06;piwind2=0.1;
		prwind1=0.68;prwind2=1.22;
		}

		if (i==11) {
                kwind1=0.16; kwind2=0.34;
                piwind1=-0.08;piwind2=0.12;
		prwind1=0.68;prwind2=1.22;
		}
	
		if (i==12) {
                kwind1=0.16; kwind2=0.38;
                piwind1=-0.1;piwind2=0.12;
		prwind1=0.68;prwind2=1.22;
		}
		
		if (i==13) {
		kwind1=0.16; kwind2=0.38;
                piwind1=-0.12;piwind2=0.12;
		prwind1=0.65;prwind2=1.22;
		}
		
		if (i==14) {
                kwind1=0.16; kwind2=0.38;
                piwind1=-0.14;piwind2=0.12;
		prwind1=0.65;prwind2=1.25;
		}
		
		if (i==15) {
                kwind1=0.16; kwind2=0.4;
                piwind1=-0.14;piwind2=0.12;
		prwind1=0.62;prwind2=1.25;
		}
		
		if (i==16) {
                kwind1=0.16; kwind2=0.4;
                piwind1=-0.14;piwind2=0.1;
		prwind1=0.6;prwind2=1.26;
		}
		
		if (i==17) {
                        kwind1=0.16; kwind2=0.4;
                        piwind1=-0.18;piwind2=0.1;
                        prwind1=0.56;prwind2=1.29;
                }

                        if (i==18) {
                kwind1=0.16; kwind2=0.4;
                piwind1=-0.18;piwind2=0.1;
                        prwind1=0.56;prwind2=1.29;
                }

                        if (i==19) {
                        kwind1=0.16; kwind2=0.4;
                        piwind1=-0.2;piwind2=0.1;
                        prwind1=0.5;prwind2=1.35;
                }

                if (i==20) {
                        kwind1=0.16; kwind2=0.4;
                        piwind1=-0.2;piwind2=0.1;
                        prwind1=0.5;prwind2=1.38;
                }

                if (i==21) {
                        kwind1=0.16; kwind2=0.38;
                        piwind1=-0.2;piwind2=0.1;
                        prwind1=0.5;prwind2=1.38;
                }

		double kback=(s[0]->Integral(kwind1,kwind2) + s[2]->Integral(kwind1,kwind2) + bg->Integral(kwind1,kwind2))/totall->Integral(kwind1,kwind2);
                double kloss=1 - s[1]->Integral(kwind1,kwind2)/s[1]->Integral(-0.5,1.9);

                double kweight=(1 - kback)/(1-kloss);

                double piback=(s[1]->Integral(piwind1,piwind2) + s[2]->Integral(piwind1,piwind2) + bg->Integral(piwind1,piwind2))/totall->Integral(piwind1,piwind2);
                double piloss=1 - s[0]->Integral(piwind1,piwind2)/s[0]->Integral(-0.5,1.9);

                double piweight=(1 - piback)/(1-piloss);

                double prback=(s[0]->Integral(prwind1,prwind2) + s[1]->Integral(prwind1,prwind2) + bg->Integral(prwind1,prwind2))/totall->Integral(prwind1,prwind2);
                double prloss=1 - s[2]->Integral(prwind1,prwind2)/s[2]->Integral(-0.5,1.9);

                double prweight=(1 - prback)/(1-prloss);

		if (i!=21) {
                ln = Form("%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n",.1*(i-1)+.5,.1*i+.5,kwind1,kwind2,kback,kloss,kweight,piwind1,piwind2,piback,piloss,piweight,prwind1,prwind2,prback,prloss,prweight);
                } else
                        ln = Form("%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f",.1*(i-1)+.5,.1*i+.5,kwind1,kwind2,kback,kloss,kweight,piwind1,piwind2,piback,piloss,piweight,prwind1,prwind2,prback,prloss,prweight);
                r = ln;
                fwrite(ln,1,r.size(),outfile);
		//s[0]->DrawF1(-0.5,1.4,"same");
                //s[1]->DrawF1(-0.5,1.4,"same");
                //s[2]->DrawF1(-0.5,1.9,"same");
                //bg->DrawF1(range[0],range[1],"same");

		//gPad->Update();
                //TPaveStats *stat=(TPaveStats*)mass->FindObject("stats");
                //stat->SetX1NDC(0.8);
                //stat->SetY1NDC(0.5);

		TLine *z=new TLine(kwind1,0.4,kwind1,600);
                TLine *z1=new TLine(kwind2,0.4,kwind2,600);
                z->SetLineColor(kRed);
                z1->SetLineColor(kRed);
		z->SetLineWidth(2);
                z1->SetLineWidth(2);
                //z->Draw("same");
                //z1->Draw("same");

		TLine *zz=new TLine(piwind1,0.4,piwind1,2000);
                TLine *zz1=new TLine(piwind2,0.4,piwind2,2000);
                zz->SetLineColor(kBlue);
                zz1->SetLineColor(kBlue);
                zz->SetLineWidth(2);
                zz1->SetLineWidth(2);
                //zz->Draw("same");
                //zz1->Draw("same");

		TLine *zv=new TLine(prwind1,0.4,prwind1,600);
                TLine *zv1=new TLine(prwind2,0.4,prwind2,600);
                zv->SetLineColor(12);
                zv1->SetLineColor(12);
                zv->SetLineWidth(2);
                zv1->SetLineWidth(2);
                //zv->Draw("same");
                //zv1->Draw("same");
                TLegend *legend=new TLegend(0.65,0.75,0.78,0.85);
                legend->AddEntry(s[0],"#pi^{+} Signal","l");
                legend->AddEntry(s[1],"K^{+} Signal","l");
                legend->AddEntry(s[2],"p Signal","l");
                legend->AddEntry(bg,"Background","l");
                legend->AddEntry(totall,"Total","l");
                //legend->Draw();
	}

	file->Close();
        fclose(outfile);

	TGraphErrors *mv1 = new TGraphErrors (mean1.size(),momentum.data(),mean1.data(),ermomentum.data(),ermean1.data());
        TGraphErrors *mv2 = new TGraphErrors (mean2.size(),momentum.data(),mean2.data(),ermomentum.data(),ermean2.data());
        TGraphErrors *mv3 = new TGraphErrors (mean3.size(),momentum.data(),mean3.data(),ermomentum.data(),ermean3.data());
	mv1->SetMarkerStyle(20);
        mv1->SetMarkerSize(.7);
        mv1->SetMarkerColor(kBlue);
        mv2->SetMarkerStyle(20);
        mv2->SetMarkerSize(.7);
        mv2->SetMarkerColor(kRed);
        mv3->SetMarkerStyle(20);
        mv3->SetMarkerSize(.7);
        mv3->SetMarkerColor(kOrange);
        mv3->SetTitle("(mass)^{2} vs momentum;GeV/Q;(GeV/Q)^{2}");

	TGraphErrors  *sv1 = new TGraphErrors (sigma1.size(),momentum.data(),sigma1.data(),ermomentum.data(),ersigma1.data());
        TGraphErrors  *sv2 = new TGraphErrors (sigma2.size(),momentum.data(),sigma2.data(),ermomentum.data(),ersigma2.data());
        TGraphErrors  *sv3 = new TGraphErrors (sigma3.size(),momentum.data(),sigma3.data(),ermomentum.data(),ersigma3.data());
        sv1->SetMarkerStyle(20);
        sv1->SetMarkerSize(.7);
        sv1->SetMarkerColor(kBlue);
        sv2->SetMarkerStyle(20);
        sv2->SetMarkerSize(.7);
        sv2->SetMarkerColor(kRed);
        sv3->SetMarkerStyle(20);
        sv3->SetMarkerSize(.7);
        sv3->SetMarkerColor(kOrange);
        sv3->SetTitle("#sigma(mass^{2}) vs momentum;GeV/Q;(GeV/Q)^{2}");
}

