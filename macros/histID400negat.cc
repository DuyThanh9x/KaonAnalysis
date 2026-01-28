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
        return skewt(x, par) + skewt(x,&par[5]) + back(x, &par[10]);
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

	bg = new TF1("bg",back,-9,5,4);
        bg->SetNpx(900);
        bg->SetLineColor(kBlack);
        bg->SetLineStyle(9);

        sum = new TF1("sum",total,-9,5,14);
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
	sum->SetParName(10,"CoeffExp");
	sum->SetParName(11,"ParExp");
	sum->SetParName(12,"Constant");
	sum->SetParName(13,"Coeff1st");
}

void histID400negat(TString name)
{
        TFile *file = TFile::Open(name);
        //FILE* outfile =fopen("particleID400negat.txt","w");
	char *ln;
 	string r;
	ln = Form("||     Momentum range      ||         Kaon window          | Kaon backgrnd   | Kaon loss   | Kaon weight   ||          Pion window         | Pion backgrnd  |  Pion loss    | Pion weight\n");
	r = ln;
	//fwrite(ln,1,r.size(),outfile);
	double par[14], range[2];
        vector<double> mean1, mean2,ermean1, ermean2,sigma1, sigma2, ersigma1, ersigma2, momentum, ermomentum;
	mean1.clear();mean2.clear();sigma1.clear();sigma2.clear();momentum.clear();
	ermean1.clear();ermean2.clear();ersigma1.clear();ersigma2.clear();ermomentum.clear();
	double kwind1,kwind2,piwind1,piwind2;
	TH1D *mass = new TH1D();
	TF1 *s[2], *bg, *totall;
	gStyle->SetOptStat(11);
        gStyle->SetOptFit(111);
        ROOT::Math::MinimizerOptions::SetDefaultMaxFunctionCalls(100000);
	for (int i = 1; i < 22; i++) {
		//if ((i <19) || (i >21)) continue;
		auto c = new TCanvas();
                c->SetLogy();
                mass = (TH1D*)file->Get(Form("mass%d",i));
                mass->SetLineColor(14);
                double mment=(.1*i+.1*(i-1))/2+.5;
		double par1[]={29,0.01,0.1,1,2,15,0.24,0.2,1.2,5,9, -6, 5, -9};
		tf1s(s,bg,totall);
                totall->SetParameters(par1);
		range[0]=-0.2;range[1]=1.5;
		if (i==21) {
			totall->SetParLimits(2,0.01,0.07);
			totall->SetParLimits(6,0.21,0.27);
			totall->SetParLimits(7,0.005,0.05);
			totall->FixParameter(7,0.05);
			totall->SetParLimits(8,0,1.5);totall->FixParameter(8,0.01);
			totall->SetParLimits(9,1,18);
			totall->FixParameter(9,18);
			totall->SetParLimits(10,2,38);
			totall->SetParLimits(11,-28,0);totall->FixParameter(12,3);range[0]=-0.22;
		}

		if (i==20) {
			totall->SetParLimits(2,0.02,0.08);
			totall->SetParLimits(6,0.21,0.27);
			totall->SetParLimits(7,0.005,0.06);
			totall->FixParameter(7,0.05);
			totall->SetParLimits(8,0,1.5);
			totall->FixParameter(8,0.01);
			totall->SetParLimits(9,1,89);
                        totall->FixParameter(9,89);
			totall->SetParLimits(10,2,38);
                        totall->SetParLimits(11,-28,0);totall->FixParameter(12,3);
		}

		if (i== 19) {
                        totall->SetParLimits(2,0.01,0.07);
                        totall->SetParLimits(6,0.21,0.27);
                        totall->SetParLimits(7,0.005,0.05);
			totall->FixParameter(7,0.05);
			totall->SetParLimits(8,0,1.5);totall->FixParameter(8,0.01);
                        totall->SetParLimits(9,1,58);
			totall->FixParameter(9,18);
			totall->SetParLimits(10,2,38);
                        totall->SetParLimits(11,-28,0);totall->FixParameter(12,3);
		}

		if (i== 18) {
                        totall->SetParLimits(2,0.01,0.07);
			totall->SetParLimits(3,0,1.5);
                        totall->FixParameter(3,0.01);
			totall->SetParLimits(6,0.21,0.27);
                        totall->SetParLimits(7,0.005,0.09);
                        totall->SetParLimits(8,0,1.5);
                        totall->FixParameter(8,0.01);
			totall->SetParLimits(9,1,58);
                        totall->FixParameter(9,18);
			totall->SetParLimits(10,2,58);
                        totall->SetParLimits(11,-28,0);totall->FixParameter(12,3);//totall->SetParLimits(12,0,138);
                }

		if (i== 17) {
			//totall->SetParLimits(3,0,1.5);
			totall->FixParameter(3,0.01);
			totall->SetParLimits(6,0.19,0.27);
                        totall->SetParLimits(7,0.02,0.07);
                        totall->SetParLimits(8,0,3.);
                        totall->FixParameter(8,0.03);
			totall->SetParLimits(9,1,58);
			totall->FixParameter(9,18);
                        totall->SetParLimits(10,1,18);
                        totall->SetParLimits(11,-28,0);
			totall->SetParLimits(12,0,138);
			totall->FixParameter(12,1);
		}

		if (i== 16) {
			totall->SetParLimits(3,0,1.5);
			totall->FixParameter(3,0.01);
			totall->SetParLimits(6,0.19,0.27);
			totall->SetParLimits(7,0.02,0.07);
			//totall->FixParameter(7,0.035);
			totall->SetParLimits(8,0,1.5);
			totall->FixParameter(8,0.01);
			totall->SetParLimits(9,1,18);
			totall->FixParameter(9,18);
			totall->SetParLimits(10,2,18);
			totall->FixParameter(10,3);
			totall->SetParLimits(11,-38,0);
			totall->SetParLimits(12,0,138);
			totall->FixParameter(12,1);
		}

		if (i== 15) {
			totall->SetParLimits(3,0,1.5);
			totall->FixParameter(3,0.01);
			totall->SetParLimits(6,0.19,0.27);
                        totall->SetParLimits(7,0.02,0.07);
			//totall->FixParameter(7,0.035);
                        totall->SetParLimits(8,0,1.9);
			totall->FixParameter(8,0.01);
                        totall->SetParLimits(9,1,38);
			//totall->FixParameter(9,18);
                        totall->SetParLimits(10,2,18);
			totall->FixParameter(10,3);
                        totall->SetParLimits(11,-38,0);
                        totall->SetParLimits(12,0,138);
			totall->FixParameter(12,1);
			//totall->SetParLimits(13,-28,0);
		}

		if (i== 14) {
			totall->SetParLimits(3,0,1.5);
			totall->FixParameter(3,0.01);
			totall->SetParLimits(6,0.19,0.27);
			totall->SetParLimits(8,0,1.9);
			totall->FixParameter(8,0.01);
			totall->SetParLimits(10,1,38);
			//totall->FixParameter(10,3);
			totall->SetParLimits(11,-28,0);
			range[0]=-0.15;
		}

		if (i== 13) {
                        //totall->SetParLimits(3,0,1.5);
                        totall->FixParameter(3,0.01);
                        totall->SetParLimits(6,0.19,0.27);
                        totall->SetParLimits(7,0.02,0.07);
			totall->SetParLimits(8,0,1.9);
                        totall->SetParLimits(9,1,58);
			totall->FixParameter(9,58);
			totall->FixParameter(8,0.01);
                        totall->SetParLimits(10,1,38);
                        totall->FixParameter(10,1); 
			totall->SetParLimits(11,-28,0);
			totall->SetParLimits(12,0,138);
			//totall->FixParameter(12,3);
                        range[0]=-0.15;
                }

		if (i== 12) {
			totall->SetParLimits(2,0.01,0.07);
			totall->SetParLimits(3,0,1.5);
			totall->SetParLimits(6,0.19,0.27);
			totall->SetParLimits(7,0.02,0.07);
			totall->SetParLimits(8,0,1.9);
			totall->FixParameter(8,0.01);
			totall->SetParLimits(9,3,139);
			totall->SetParLimits(10,1,58);
			totall->FixParameter(10,19);
			totall->FixParameter(9,139);
                        totall->SetParLimits(11,-28,0);
			totall->SetParLimits(12,0,138);
			totall->FixParameter(12,1);
			range[0]=-0.1;
		}

		if (i== 11) {
			totall->SetParLimits(3,0,1.5);
			totall->SetParLimits(6,0.19,0.27);
			totall->SetParLimits(7,0.02,0.07);
			totall->SetParLimits(8,0,1.9);
			totall->FixParameter(8,0.01);
			totall->SetParLimits(9,1,139);
			totall->FixParameter(9,139);
			totall->SetParLimits(10,1,118);
			//totall->FixParameter(10,19);
                        totall->SetParLimits(11,-28,0);
			totall->SetParLimits(12,0,138);
			//totall->FixParameter(12,1);
			range[0]=-0.09;
		}

		if (i== 10) {
			totall->SetParLimits(2,0.009,0.03);
			totall->SetParLimits(3,0,1.5);
                        totall->SetParLimits(6,0.19,0.27);
			totall->SetParLimits(7,0.005,0.05);
			totall->FixParameter(7,0.03);
			totall->SetParLimits(8,0,3.9);
			totall->SetParLimits(9,1,139);
			totall->FixParameter(9,139);
			totall->SetParLimits(10,1,138);
			totall->SetParLimits(11,-28,0);
			totall->FixParameter(10,118);
			range[0]=-0.05;
                }

		if (i== 9) {
			totall->SetParLimits(2,0.009,0.03);
			totall->SetParLimits(3,0,1.5);
                        totall->SetParLimits(6,0.19,0.27);
			totall->SetParLimits(7,0.01,0.05);
			totall->SetParLimits(8,0,1.9);
			totall->SetParLimits(9,1,139);
			totall->FixParameter(9,139);
			totall->SetParLimits(10,1,138);
			totall->FixParameter(10,138);
			totall->SetParLimits(11,-28,0);
			range[0]=-0.05;
                }

		if (i== 8) {
			totall->SetParLimits(2,0.009,0.03);
			totall->SetParLimits(3,0,1.5);
                        totall->SetParLimits(6,0.19,0.27);
			totall->SetParLimits(7,0.01,0.05);
			totall->SetParLimits(8,0,1.9);
			totall->SetParLimits(9,1,119);
			totall->FixParameter(9,119);
			totall->SetParLimits(10,1,118);
			totall->FixParameter(10,118);
			totall->SetParLimits(11,-28,0);
			range[0]=-0.05;
                }

		if (i== 7) {
			totall->SetParLimits(2,0.009,0.03);
			totall->SetParLimits(3,0,1.5);
                        totall->SetParLimits(6,0.19,0.27);
			totall->SetParLimits(7,0.01,0.05);
			totall->SetParLimits(8,0,3.9);totall->FixParameter(8,0.01);
			totall->SetParLimits(9,1,119);
			totall->SetParLimits(10,1,158);
			//totall->FixParameter(10,138);
			totall->SetParLimits(11,-28,0);
			range[0]=-0.05;
                }

		if (i== 6) {
			totall->SetParLimits(2,0.009,0.03);
			totall->SetParLimits(3,0,1.5);
                        totall->SetParLimits(6,0.19,0.27);
			totall->SetParLimits(7,0.01,0.05);
			totall->SetParLimits(8,0,3.9);
			totall->SetParLimits(9,1,119);
			totall->FixParameter(9,119);
                        totall->SetParLimits(10,1,138);
                        totall->FixParameter(10,138);
			totall->SetParLimits(11,-28,0);
			range[0]=-0.05;
                }

		if (i== 5) {
			totall->SetParLimits(2,0.005,0.05);
			totall->SetParLimits(3,0,1.5);
                        totall->SetParLimits(6,0.19,0.27);
			totall->SetParLimits(7,0.01,0.05);
			totall->SetParLimits(8,0,3.9);
			totall->SetParLimits(9,1,119);
                        totall->FixParameter(9,119);
			totall->SetParLimits(10,1,158);
                        totall->SetParLimits(11,-28,0);
			totall->SetParLimits(12,0,119);
			totall->FixParameter(12,1);
			range[0]=-0.03;
                }

		if (i== 4) {
			totall->SetParLimits(2,0.007,0.03);
			//totall->FixParameter(2,0.009);
			totall->SetParLimits(3,0,1.);
			totall->FixParameter(3,0.01);
			totall->SetParLimits(6,0.17,0.27);
			totall->SetParLimits(7,0.01,0.03);
			totall->SetParLimits(8,1,3.9);
			totall->SetParLimits(9,1,119);
			//totall->FixParameter(9,119);
                        totall->SetParLimits(10,1,138);
			totall->FixParameter(10,158);
                        totall->SetParLimits(11,-28,0);
			//totall->SetParLimits(12,0,119);
			range[0]=-0.03;
                }

		if (i== 3) {
			totall->SetParLimits(2,0.005,0.03);
			totall->SetParLimits(3,0,1.5);
                       totall->FixParameter(3,0.01); 
			totall->SetParLimits(6,0.17,0.27);
			totall->SetParLimits(7,0.01,0.05);
			totall->SetParLimits(8,0,3.);
			totall->SetParLimits(9,1,119);
			totall->FixParameter(9,119);
                        totall->SetParLimits(10,1,188);totall->FixParameter(10,158);
                        totall->SetParLimits(11,-28,0);
			range[0]=-0.03;
                }

		if (i== 2) {
			totall->SetParLimits(2,0.003,0.03);
			totall->SetParLimits(3,0,1.5);
                        totall->SetParLimits(6,0.19,0.25);
			totall->SetParLimits(7,0.01,0.05);
			totall->SetParLimits(8,0,3.);
			//totall->FixParameter(8,3.3);
			totall->SetParLimits(9,1,119);
			totall->FixParameter(9,119);
                        totall->SetParLimits(10,1,158);
                        totall->SetParLimits(11,-28,0);
			range[0]=-0.03;
                }

		if (i== 1) {
			totall->SetParLimits(2,0.003,0.03);
			totall->SetParLimits(3,0,1.5);
                        totall->SetParLimits(6,0.19,0.25);
			totall->FixParameter(6,0.25);
			totall->SetParLimits(7,0.01,0.05);
			totall->FixParameter(7,0.01);
			totall->SetParLimits(8,0,3.);
			totall->FixParameter(8,0.3);
			totall->SetParLimits(9,1,119);
			totall->FixParameter(9,119);
                        totall->SetParLimits(10,1,118);
                        totall->SetParLimits(11,-28,0);
			range[0]=-0.03;
                }
		mass->Fit(totall,"SEM","",range[0],range[1]);
		totall->GetParameters(par);
                s[0]->SetParameters(par);
                s[1]->SetParameters(&par[5]);
		bg->SetParameters(&par[10]);
		mean1.push_back(totall->GetParameter(1));
                ermean1.push_back(totall->GetParError(1));
                sigma1.push_back(totall->GetParameter(2));
                ersigma1.push_back(totall->GetParError(2));
                mean2.push_back(totall->GetParameter(6));
                ermean2.push_back(totall->GetParError(6));
                sigma2.push_back(totall->GetParameter(7));
                        ersigma2.push_back(totall->GetParError(7));
			momentum.push_back(mment);
                ermomentum.push_back(0);
		s[0]->DrawF1(-0.5,1.4,"same");
                s[1]->DrawF1(-0.5,1.4,"same");
		bg->DrawF1(range[0],range[1],"same");
		gPad->Update();
                TPaveStats *stat=(TPaveStats*)mass->FindObject("stats");
                stat->SetX1NDC(0.8);
                stat->SetY1NDC(0.5);
		TLegend *legend=new TLegend(0.65,0.75,0.78,0.85);
                legend->AddEntry(s[0],"#pi^{-} Signal","l");
                legend->AddEntry(s[1],"K^{-} Signal","l");
                //legend->AddEntry(s[2],"p Signal","l");
                legend->AddEntry(bg,"Background","l");
                legend->AddEntry(totall,"Total","l");
                legend->Draw();
		
		if (i==2) {
                        kwind1=0.23; kwind2=0.28;
                        piwind1=-0.01;piwind2=0.07;
                }
		
		if (i==3) {
                        kwind1=0.23; kwind2=0.28;
                        piwind1=-0.01;piwind2=0.07;
                }

                if (i==4) {
                        kwind1=0.23; kwind2=0.28;
                        piwind1=-0.017;piwind2=0.07;
                }

                if (i==5) {
                        kwind1=0.22; kwind2=0.28;
                        piwind1=-0.029;piwind2=0.09;
                }
                if (i==6) {
                        kwind1=0.23; kwind2=0.28;
                        piwind1=-0.04;piwind2=0.09;
                }
                if (i==7) {
                kwind1=0.22; kwind2=0.28;
                piwind1=-0.05;piwind2=0.09;
                }
                if (i==8) {
                        kwind1=0.22; kwind2=0.28;
                        piwind1=-0.05;piwind2=0.11;
                }
                if (i==9) {
                kwind1=0.21; kwind2=0.28;
                piwind1=-0.07;piwind2=0.12;
                }
                if (i==10) {
                kwind1=0.2; kwind2=0.3;
                piwind1=-0.07;piwind2=0.14;
                }
                if (i==11) {
                kwind1=0.18; kwind2=0.3;
                piwind1=-0.09;piwind2=0.17;
                }
                if (i==12) {
                kwind1=0.18; kwind2=0.34;
                piwind1=-0.12;piwind2=0.17;
                }

		if (i==13) {
                        kwind1=0.18; kwind2=0.38;
                piwind1=-0.15;piwind2=0.17;
                }
                if (i==14) {
                kwind1=0.2; kwind2=0.4;
                piwind1=-0.15;piwind2=0.17;
                }
                if (i==15) {
                kwind1=0.19; kwind2=0.38;
                piwind1=-0.18;piwind2=0.17;
                }
                if (i==16) {
                kwind1=0.18; kwind2=0.39;
                piwind1=-0.18;piwind2=0.15;
                }
                if (i==17) {
                        kwind1=0.19; kwind2=0.4;
                        piwind1=-0.17;piwind2=0.15;
                }
                if (i==18) {
                kwind1=0.19; kwind2=0.4;
                piwind1=-0.17;piwind2=0.15;
                }
                if (i==19) {
                        kwind1=0.18; kwind2=0.35;
                        piwind1=-0.19;piwind2=0.14;
                }
                if (i==20) {
                        kwind1=0.18; kwind2=0.35;
                        piwind1=-0.2;piwind2=0.15;
                }
                if (i==21) {
                        kwind1=0.18; kwind2=0.37;
                        piwind1=-0.22;piwind2=0.15;
                }

		if (i>1) {
                double kback=(s[0]->Integral(kwind1,kwind2) + bg->Integral(kwind1,kwind2))/totall->Integral(kwind1,kwind2);
                double kloss=1 - s[1]->Integral(kwind1,kwind2)/s[1]->Integral(-0.5,1.5);

                double kweight=(1 - kback)/(1-kloss);

                double piback=(s[1]->Integral(piwind1,piwind2)  + bg->Integral(piwind1,piwind2))/totall->Integral(piwind1,piwind2);
                double piloss=1 - s[0]->Integral(piwind1,piwind2)/s[0]->Integral(-0.5,1.5);

                double piweight=(1 - piback)/(1-piloss);

                cout<<">>>> "<<.1*(i-1)+.5<<"\t"<<.1*i+.5<<"\t"<<kwind1<<"\t"<<kwind2<<"\t"<<kback<<"\t"<<kloss<<"\t"<<kweight<<"\t"<<piwind1<<"\t"<<piwind2<<"\t"<<piback<<"\t"<<piloss<<"\t"<<piweight<<endl<<endl;
		if (i!=21) {
		ln = Form("%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n",.1*(i-1)+.5,.1*i+.5,kwind1,kwind2,kback,kloss,kweight,piwind1,piwind2,piback,piloss,piweight);
		} else
			ln = Form("%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f",.1*(i-1)+.5,.1*i+.5,kwind1,kwind2,kback,kloss,kweight,piwind1,piwind2,piback,piloss,piweight);
		r = ln;
       		//fwrite(ln,1,r.size(),outfile);
	
                TLine *z=new TLine(kwind1,0.,kwind1,600);
                TLine *z1=new TLine(kwind2,0.,kwind2,600);
                z->SetLineColor(kRed);
                z1->SetLineColor(kRed);
                z->SetLineWidth(2);
                z1->SetLineWidth(2);
                z->Draw("same");
                z1->Draw("same");

                TLine *zz=new TLine(piwind1,0.,piwind1,2000);
                TLine *zz1=new TLine(piwind2,0.,piwind2,2000);
                zz->SetLineColor(kBlue);
                zz1->SetLineColor(kBlue);
                zz->SetLineWidth(2);
                zz1->SetLineWidth(2);
                zz->Draw("same");
                zz1->Draw("same");
                }
		c->SaveAs(Form("25/mnegat%d.pdf",i));
	}
//	file->Close();
	//fclose(outfile);
	auto c = new TCanvas();
	TGraphErrors *mv1 = new TGraphErrors (mean1.size(),momentum.data(),mean1.data(),ermomentum.data(),ermean1.data());
        TGraphErrors *mv2 = new TGraphErrors (mean2.size(),momentum.data(),mean2.data(),ermomentum.data(),ermean2.data());
	mv1->SetMarkerStyle(20);
        mv1->SetMarkerSize(.7);
        mv1->SetMarkerColor(kBlue);
        mv2->SetMarkerStyle(20);
        mv2->SetMarkerSize(.7);
        mv2->SetMarkerColor(kRed);
	//mv2->SetTitle("(mass)^{2} vs momentum;GeV/Q;(GeV/Q)^{2}");
	mv2->GetYaxis()->SetRangeUser(-.2,0.3);
	//mv2->Draw("AP");
	mv1->SetTitle("#pi^{-} (mass)^{2} vs momentum;GeV/Q;(GeV/Q)^{2}");
	mv1->Draw("AP");//->Draw("SAMEP");
	auto legend = new TLegend();
        legend->AddEntry(mv1,"#pi^{+}","p");
        legend->AddEntry(mv2,"K^{+}","p");
	//legend->Draw();
	c->SaveAs("25/msnegat4.pdf");
	auto vc = new TCanvas();
	TGraphErrors  *sv1 = new TGraphErrors (sigma1.size(),momentum.data(),sigma1.data(),ermomentum.data(),ersigma1.data());
        TGraphErrors  *sv2 = new TGraphErrors (sigma2.size(),momentum.data(),sigma2.data(),ermomentum.data(),ersigma2.data());
	sv1->SetMarkerStyle(20);
        sv1->SetMarkerSize(.7);
        sv1->SetMarkerColor(kBlue);
        sv2->SetMarkerStyle(20);
        sv2->SetMarkerSize(.7);
        sv2->SetMarkerColor(kRed);
	sv1->SetTitle("#pi^{-} #sigma(mass^{2}) vs momentum;GeV/Q;(GeV/Q)^{2}");
	//sv2->Draw("AP");
        sv1->GetYaxis()->SetRangeUser(0,0.07);
	sv1->Draw("AP");vc->SaveAs("25/smsnegat4.pdf");
	//->Draw("SAMEP");

	auto legend1 = new TLegend();
        legend1->AddEntry(sv1,"#pi^{+}","p");
        legend1->AddEntry(sv2,"K^{+}","p");
	//legend1->Draw();
}
