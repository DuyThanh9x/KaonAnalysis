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

void histID700(TString name)
{
        TFile *file = TFile::Open(name);
        FILE* outfile =fopen("particleID700.txt","w");
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
        for (int i = 2; i < 22; i++) {
                //if ((i <19) || (i >21)) continue;
                //auto c = new TCanvas();
                //c->SetLogy();
                mass = (TH1D*)file->Get(Form("mass%d",i));
                mass->SetLineColor(14);
                double mment=(.1*i+.1*(i-1))/2+.5;
		smvl[0]=0.0073-0.0129*mment+0.0168*mment*mment;
		smvl[1]=0.024-0.028*mment+0.019*mment*mment;
		smvl[2]=0.082-0.054*mment+0.024*mment*mment;

		meanvl[0]=0.031-0.018*TMath::Power(mment,1);
                meanvl[1]=-1.93+2.168*TMath::Power(mment,-0.0092);
                meanvl[2]=36.05-35.16*TMath::Power(mment,0.000709);
		

		mean11.push_back(meanvl[0]); mean21.push_back(meanvl[1]);mean31.push_back(meanvl[2]);
                sigma11.push_back(smvl[0]); sigma21.push_back(smvl[1]);sigma31.push_back(smvl[2]);

		cout<<endl<<"++++++++ "<<i<<" Mean "<<meanvl[0]<<'\t'<<meanvl[1]<<'\t'<<meanvl[2]<<" Sm "<<smvl[0]<<'\t'<<smvl[1]<<'\t'<<smvl[2]<<endl<<endl;

		double par1[]={29, meanvl[0], smvl[0], 0.5, 2, 15, meanvl[1], smvl[1], 0.5, 1, 39, meanvl[2], smvl[2], 0.9, 2, 9, -6, 5, -9};
                tf1s(s,bg,totall);
                totall->SetParameters(par1);

		range[0]=-0.034;range[1]=1.68;
		cout<<endl<<"++++++++ "<<i<<endl<<endl;

		if (i==2) {
		totall->FixParameter(3,0.02);
		totall->SetParLimits(6,0.21,0.29);
		totall->FixParameter(7,0.01);
		totall->FixParameter(8,0.02);
		totall->FixParameter(9,3);
		totall->SetParLimits(16,-18,0);
		}

		if (i== 3) {
			totall->FixParameter(3,0.01);
			totall->FixParameter(4,2.9);
			totall->SetParLimits(6,0.21,0.29);
			totall->SetParLimits(8,0.,1.5);
			totall->FixParameter(9,3);
			totall->SetParLimits(16,-18,0);
		}
		if (i==4) {
			totall->FixParameter(3,0.02);
			totall->FixParameter(9,3);
		}

		if (i==5) {
			totall->FixParameter(3,0.01);
			totall->FixParameter(8,0.01);
		}
	
		if (i==6) {
			totall->SetParLimits(4,2,39);
			totall->SetParLimits(8,0.,2.9);
			totall->SetParLimits(9,2,41);
			totall->FixParameter(9,2);
		}

		if (i==7) {
			totall->SetParLimits(3,0,1);
			totall->SetParLimits(4,1,39);
			totall->FixParameter(4,2.8);
			totall->SetParLimits(6,0.21,0.29);
			totall->SetParLimits(8,0.,1.9);
			totall->FixParameter(8,0.01);
			totall->FixParameter(9,3);
		}

		if (i==8) {
			totall->SetParLimits(3,0,1);
			totall->SetParLimits(4,1,39);
			totall->FixParameter(4,2.8);
			totall->SetParLimits(8,0.,2.9);
			totall->FixParameter(8,0.02);
			totall->SetParLimits(9,1,41);
		}

		if (i==9) {
			totall->SetParLimits(3,0,1);
			totall->SetParLimits(4,2,39);
			totall->FixParameter(4,2.8);
			totall->SetParLimits(8,0.,2.9);
			totall->FixParameter(8,0.0);
			totall->SetParLimits(9,2,41);
			totall->FixParameter(9,2);
			totall->SetParameter(16,-0.5);
		}
		
                if (i==10) {
			totall->SetParLimits(3,0,1);
			totall->FixParameter(4,2.8);
			totall->SetParLimits(8,0.,2.9);
			totall->FixParameter(8,0.02);
			totall->SetParLimits(9,2,41);
			totall->FixParameter(9,2);
			totall->SetParameter(16,-0.5);
		}

		if (i==11) {
			totall->SetParLimits(3,0,1);
			totall->FixParameter(4,2.8);
			totall->SetParLimits(8,0.,2.9);
			totall->FixParameter(8,0.02);
			totall->SetParLimits(9,2,41);
			totall->FixParameter(9,2);
			totall->SetParameter(16,-0.5);
		}

		if (i==12) {
                        totall->SetParLimits(3,0,1);
                        totall->FixParameter(4,2.8);
                        totall->SetParLimits(8,0.,2.9);
                        totall->FixParameter(8,0.02);
                        totall->SetParLimits(9,2,41);
                        totall->FixParameter(9,2);
                        totall->SetParameter(16,-0.5);
		}

		if (i==13) {
                        totall->SetParLimits(3,0,2.9);
                        totall->FixParameter(4,2.8);
			totall->SetParLimits(6,0.21,0.29);
                        totall->SetParLimits(8,0.,2.9);
                        totall->FixParameter(8,0.02);
                        totall->SetParLimits(9,2,41);
                        totall->FixParameter(9,2);
                        totall->SetParameter(16,-0.5);
		}

		if (i==14) {
                        totall->SetParLimits(3,0,2);
                        totall->SetParLimits(6,0.21,0.29);
			totall->SetParLimits(8,0.,2.9);
                        totall->FixParameter(9,2.5);
			totall->SetParLimits(15,2,39);
                        totall->SetParLimits(16,-28,0);
			totall->SetParLimits(18,-28,0);
			range[0]=-0.1;
		}

		if (i==15) {
                        totall->SetParLimits(3,0,2.9);
			totall->SetParLimits(6,0.21,0.29);
                        totall->SetParLimits(8,0.,2.9);
                        totall->FixParameter(9,2.5);
			totall->SetParLimits(15,25,43);
			totall->FixParameter(15,43);
			totall->SetParLimits(16,-28,0);
			totall->SetParLimits(18,-28,0);
			range[0]=-0.1;
		}

		if (i==16) {
			totall->SetParLimits(3,0,2.9);
			totall->SetParLimits(6,0.21,0.29);
			totall->SetParLimits(8,0.,2.9);
			totall->FixParameter(9,2.9);
			totall->SetParLimits(15,2,39);
			totall->SetParLimits(16,-28,0);
			range[0]=-0.1;
		}

		if (i==17) {
			totall->SetParLimits(3,0,2.9);
			totall->SetParLimits(6,0.21,0.29);
			totall->SetParLimits(8,0.,2.9);
			totall->SetParLimits(9,2,41);
			totall->FixParameter(9,3);
			totall->SetParLimits(15,2,39);
			totall->FixParameter(15,43);
			totall->SetParLimits(16,-28,0);
			range[0]=-0.1;
		}

		if (i==18) {
			totall->SetParLimits(3,0,2.9);
			totall->SetParLimits(6,0.21,0.29);
			totall->SetParLimits(8,0.,2.9);
			totall->SetParLimits(9,2,41);
			totall->FixParameter(9,3);
			totall->FixParameter(15,43);
			totall->SetParLimits(16,-28,0);
			range[0]=-0.1;
		}

		if (i==19) {
			totall->SetParLimits(3,0,2.9);
			totall->SetParLimits(6,0.21,0.29);
			totall->SetParLimits(8,0.,2.9);
			totall->SetParLimits(9,2,41);
			totall->FixParameter(9,3.2);
			totall->SetParLimits(15,2,39);
			totall->SetParLimits(16,-28,0);
		range[0]=-0.1;
		}

		if (i==20) {
                        totall->SetParLimits(3,0,2.9);			
			totall->SetParLimits(6,0.21,0.29);
                        totall->SetParLimits(8,0.,2.9);
                        totall->FixParameter(9,3.5);
			totall->SetParLimits(15,2,39);
                        totall->SetParLimits(16,-28,0);
                range[0]=-0.1;
                }

		if (i==21) {
                        totall->SetParLimits(3,0,2.9);
			totall->SetParLimits(6,0.21,0.29);
                        totall->SetParLimits(8,0.,2.9);
                        totall->FixParameter(9,3.2);
			totall->SetParLimits(15,2,39);
                        totall->SetParLimits(16,-28,0);
                range[0]=-0.1;
                }
		
		mass->Fit(totall,"SEMq0","",range[0],range[1]);

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
                                skew.push_back(totall->GetParameter(13));
                                erskew.push_back(totall->GetParError(13));

				momentum.push_back(mment);
                ermomentum.push_back(0);

		//s[0]->DrawF1(-0.5,1.4,"same");
                //s[1]->DrawF1(-0.5,1.4,"same");
                //s[2]->DrawF1(-0.5,range[1],"same");
                //bg->DrawF1(range[0],range[1],"same");

		//gPad->Update();
                //TPaveStats *stat=(TPaveStats*)mass->FindObject("stats");
                //stat->SetX1NDC(0.8);
                //stat->SetY1NDC(0.5);

		TLegend *legend=new TLegend(0.65,0.75,0.78,0.85);
                legend->AddEntry(s[0],"#pi^{+} Signal","l");
                legend->AddEntry(s[1],"K^{+} Signal","l");
                legend->AddEntry(s[2],"p Signal","l");
                legend->AddEntry(bg,"Background","l");
                legend->AddEntry(totall,"Total","l");
                //legend->Draw();

		if (i==2) {
                        kwind1=0.23; kwind2=0.3;
                        piwind1=-0.0;piwind2=0.06;
                        prwind1=0.8;prwind2=1.32;
                }

                if (i==3) {
                        kwind1=0.22; kwind2=0.32;
                        piwind1=-0.005;piwind2=0.074;
                        prwind1=0.8;prwind2=1.22;
                }

                if (i==4) {
                        kwind1=0.22; kwind2=0.32;
                        piwind1=-0.01;piwind2=0.06;
                        prwind1=0.78;prwind2=1.18;
                }

                if (i==5) {
                        kwind1=0.23; kwind2=0.3;
                        piwind1=-0.02;piwind2=0.068;
                        prwind1=0.78;prwind2=1.15;
		}

		if (i==6) {
                        kwind1=0.22; kwind2=0.3;
                        piwind1=-0.03;piwind2=0.084;
                        prwind1=0.78;prwind2=1.15;
			}

		if (i==7) {
                kwind1=0.22; kwind2=0.3;
                piwind1=-0.035;piwind2=0.1;
		prwind1=0.78;prwind2=1.15;
		}

		if (i==8) {
                        kwind1=0.22; kwind2=0.3;
                        piwind1=-0.02;piwind2=0.12;
                        prwind1=0.76;prwind2=1.15;
		}

		if (i==9) {
                kwind1=0.22; kwind2=0.32;
                piwind1=-0.03;piwind2=0.14;
		prwind1=0.76;prwind2=1.15;
		}

		if (i==10) {
                kwind1=0.21; kwind2=0.32;
                piwind1=-0.03;piwind2=0.14;
		prwind1=0.74;prwind2=1.2;
		}

		if (i==11) {
                kwind1=0.2; kwind2=0.32;
                piwind1=-0.03;piwind2=0.17;
                        prwind1=0.74;prwind2=1.22;
                }

                if (i==12) {
                kwind1=0.19; kwind2=0.34;
                piwind1=-0.05;piwind2=0.15;
                        prwind1=0.72;prwind2=1.24;
                }

                if (i==13) {
                kwind1=0.19; kwind2=0.34;
                piwind1=-0.05;piwind2=0.15;
                        prwind1=0.7;prwind2=1.3;
                }

                if (i==14) {
                kwind1=0.17; kwind2=0.43;
                piwind1=-0.09;piwind2=0.15;
                        prwind1=0.7;prwind2=1.3;
                }

                if (i==15) {
                kwind1=0.19; kwind2=0.4;
                piwind1=-0.09;piwind2=0.15;
                        prwind1=0.68;prwind2=1.35;
                }

                if (i==16) {
                kwind1=0.18; kwind2=0.4;
                piwind1=-0.09;piwind2=0.15;
                        prwind1=0.67;prwind2=1.39;
                }

                if (i==17) {
                        kwind1=0.19; kwind2=0.43;
                        piwind1=-0.09;piwind2=0.15;
			prwind1=0.63;prwind2=1.47;
		}

                        if (i==18) {
                kwind1=0.19; kwind2=0.45;
                piwind1=-0.1;piwind2=0.15;
                        prwind1=0.62;prwind2=1.55;
                }

                        if (i==19) {
                        kwind1=0.19; kwind2=0.45;
                        piwind1=-0.1;piwind2=0.15;
                        prwind1=0.6;prwind2=1.57;
                }

                if (i==20) {
                        kwind1=0.175; kwind2=0.43;
                        piwind1=-0.085;piwind2=0.15;
                        prwind1=0.6;prwind2=1.55;
                }

                if (i==21) {
                        kwind1=0.18; kwind2=0.42;
                        piwind1=-0.08;piwind2=0.15;
                        prwind1=0.58;prwind2=1.55;
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

		cout<<">>>> K window "<<.1*(i-1)+.5<<"\t"<<.1*i+.5<<"\t"<<kwind1<<"\t"<<kwind2<<"\t"<<kback<<"\t"<<kloss<<"\t"<<kweight<<"\t"<<piwind1<<"\t"<<piwind2<<"\t"<<piback<<"\t"<<piloss<<"\t"<<piweight<<"\t"<<prwind1<<"\t"<<prwind2<<"\t"<<prback<<"\t"<<prloss<<"\t"<<prweight<<endl<<" >>>> Kback = "<<kback<<" klss = "<<kloss<<" Kweight = "<<kweight<<endl
                        <<">>>> pi window "<<piwind1<<"\t"<<piwind2<<"\t"<<piback<<"\t"<<piloss<<"\t"<<piweight<<">>>>>> piback "<<piback<<" pilss = "<<piloss<<" piweight = "<<piweight<<endl
                        <<">>>> pr window "<<prwind1<<"\t"<<prwind2<<"\t"<<prback<<"\t"<<prloss<<"\t"<<prweight<<">>>>>> prback "<<prback<<" prlss = "<<prloss<<" prweight = "<<prweight<<endl<<endl;

		if (i!=21) {
                ln = Form("%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n",.1*(i-1)+.5,.1*i+.5,kwind1,kwind2,kback,kloss,kweight,piwind1,piwind2,piback,piloss,piweight,prwind1,prwind2,prback,prloss,prweight);
                } else
                        ln = Form("%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f",.1*(i-1)+.5,.1*i+.5,kwind1,kwind2,kback,kloss,kweight,piwind1,piwind2,piback,piloss,piweight,prwind1,prwind2,prback,prloss,prweight);
                r = ln;
                fwrite(ln,1,r.size(),outfile);
		TLine *z=new TLine(kwind1,0.4,kwind1,600);
                TLine *z1=new TLine(kwind2,0.4,kwind2,600);
                z->SetLineColor(kRed);
                z1->SetLineColor(kRed);
                //z->SetLineColor(12);
                //z1->SetLineColor(12);
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

		//c->SaveAs(Form("25/m%d700.pdf",i));
	}
		file->Close();
        fclose(outfile);
	//auto c = new TCanvas();
	TGraphErrors *mv1 = new TGraphErrors (mean1.size(),momentum.data(),mean1.data(),ermomentum.data(),ermean1.data());
        TGraphErrors *mv2 = new TGraphErrors (mean2.size(),momentum.data(),mean2.data(),ermomentum.data(),ermean2.data());
        TGraphErrors *mv3 = new TGraphErrors (mean3.size(),momentum.data(),mean3.data(),ermomentum.data(),ermean3.data());
	
	TF1 *mv=new TF1("mv","[0]+[1]*TMath::Power(x,[2])");
        mv->SetParameters(0.02,0.02,1.1);
        //mv1->Fit(mv,"EM","",0.9,2.2);
        mv->SetParameter(0,0.2);
        //mv2->Fit(mv,"EM","",0.9,2.2);
        mv->SetParameter(0,0.9);
        //mv3->Fit(mv,"EM","",0.9,2.2);

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
        mv3->GetYaxis()->SetRangeUser(-.2,1.2);
        mv3->GetXaxis()->SetRangeUser(.2,5.9);

	//mv3->Draw("AP");
	//mv1->Draw("SAMEP");
        //mv2->Draw("SAMEP");

	auto legend = new TLegend();
        legend->AddEntry(mv1,"#pi^{+}","p");
        legend->AddEntry(mv2,"K^{+}","p");
        legend->AddEntry(mv3,"p","p");
	//legend->Draw();

	//auto c1 = new TCanvas();
	TGraph *mr1 = new TGraph (mean11.size(),momentum.data(),mean11.data());
        mr1->SetMarkerStyle(20);
        mr1->SetMarkerSize(.7);
	mr1->SetTitle("#pi^{+} (mass)^{2} vs momentum;GeV/Q;(GeV/Q)^{2}");
	//mr1->Draw("AP");
        //mv1->Draw("SAMEP");
	
	auto legend2 = new TLegend();
        legend2->AddEntry(mv1,"Fit value","p");
        legend2->AddEntry(mr1,"Initial value","p");
	//legend2->Draw();

	//auto c2 = new TCanvas();
	TGraph *mr2 = new TGraph (mean21.size(),momentum.data(),mean21.data());
        mr2->SetMarkerStyle(20);
        mr2->SetMarkerSize(.7);
	mr2->SetTitle("K^{+} (mass)^{2} vs momentum;GeV/Q;(GeV/Q)^{2}");
	//mr2->Draw("AP");
        //mv2->Draw("SAMEP");

	auto legend3 = new TLegend();
        legend3->AddEntry(mv2,"Fit value","p");
        legend3->AddEntry(mr2,"Initial value","p");

	//auto c4 = new TCanvas();
	TGraph *mr3 = new TGraph (mean31.size(),momentum.data(),mean31.data());
        mr3->SetMarkerStyle(20);
        mr3->SetMarkerSize(.7);
	mr3->SetTitle("proton (mass)^{2} vs momentum;GeV/Q;(GeV/Q)^{2}");
	//mr3->Draw("AP");
        //mv3->Draw("SAMEP");
	
	auto legend4 = new TLegend();
        legend4->AddEntry(mv3,"Fit value","p");
        legend4->AddEntry(mr3,"Initial value","p");

	//auto c5 = new TCanvas();
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
	sv3->GetYaxis()->SetRangeUser(-.02,0.15);
	//sv3->Draw("AP");
        //sv1->Draw("SAMEP");
        //sv2->Draw("SAMEP");
	
	auto legend1 = new TLegend();
        legend1->AddEntry(sv1,"#pi^{+}","p");
        legend1->AddEntry(sv2,"K^{+}","p");
        legend1->AddEntry(sv3,"p","p");
	//legend1->Draw();

	TF1 *smv=new TF1("smv","[0]+[1]*x+[2]*x*x");
	smv->SetParameters(0.5,0.2,2);
	//sv1->Fit(smv,"EM","",0.6,1.9);
        //sv2->Fit(smv,"EM","",0.9,1.6);
	//sv3->Fit(smv,"EM","",0.9,2.5);
	//auto c11 = new TCanvas();
	TGraphErrors *sk = new TGraphErrors (skew.size(),momentum.data(),skew.data(),ermomentum.data(),erskew.data());
        sk->SetMarkerStyle(20);
        sk->SetMarkerSize(.7);
        sk->SetTitle("Proton distribution skewness vs momentum;GeV/Q");
	//sk->Draw("AP");

	//auto c12 = new TCanvas();	
	TGraph *sr1 = new TGraph (sigma11.size(),momentum.data(),sigma11.data());
        sr1->SetMarkerStyle(20);
        sr1->SetMarkerSize(.7);

	sr1->SetTitle("#pi^{+} #sigma(mass^{2}) vs momentum;GeV/Q;(GeV/Q)^{2}");

	//sr1->Draw("AP");
        //sv1->Draw("SAMEP");

	//auto c14 = new TCanvas();
	TGraph *sr2 = new TGraph (sigma21.size(),momentum.data(),sigma21.data());
        sr2->SetMarkerStyle(20);
        sr2->SetMarkerSize(.7);
	sr2->SetTitle("K^{+} #sigma(mass^{2}) vs momentum;GeV/Q;(GeV/Q)^{2}");

	//sr2->Draw("AP");
        //sv2->Draw("SAMEP");

	auto legend11 = new TLegend();
        legend11->AddEntry(sv2,"Fit value","p");
        legend11->AddEntry(sr2,"Initial value","p");

	//auto c15 = new TCanvas();
	TGraph *sr3 = new TGraph (sigma31.size(),momentum.data(),sigma31.data());
        sr3->SetMarkerStyle(20);
        sr3->SetMarkerSize(.7);
	sr3->SetTitle("proton #sigma(mass)^{2} vs momentum;GeV/Q;(GeV/Q)^{2}");

	//sr3->Draw("AP");
        //sv3->Draw("SAMEP");
	
	auto legend12 = new TLegend();
        legend12->AddEntry(sv3,"Fit value","p");
        legend12->AddEntry(sr3,"Initial value","p");

}
