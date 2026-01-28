void read(TString name,vector<double> &KaonLeftLimit, vector<double> &KaonRightLimit, vector<double> &KaonWeight,
                vector<double> &PionLeftLimit,vector<double> &PionRightLimit,vector<double> &PionWeight,vector<double> &ProtonLeftLimit,vector<double> &ProtonRightLimit,vector<double> &ProtonWeight) {
        char text[300];
	double id[18];
        FILE* file =fopen(name,"r");
	fgets(text,300,file);
	while ((file)) {
		fgets(text,300,file);
		if (feof(file))
        		break;
		sscanf(text,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf",&id[0],&id[1],&id[2],&id[3],&id[4],&id[5],&id[6],&id[7],&id[8],&id[9],&id[10],&id[11],&id[12],&id[13],&id[14],&id[15],&id[16]);
		KaonLeftLimit.push_back(id[2]);KaonRightLimit.push_back(id[3]);KaonWeight.push_back(id[6]);
                PionLeftLimit.push_back(id[7]);PionRightLimit.push_back(id[8]);PionWeight.push_back(id[11]);ProtonLeftLimit.push_back(id[12]);ProtonRightLimit.push_back(id[13]);ProtonWeight.push_back(id[16]);
	}
	fclose(file);
}

void ragidityVsTransverseMment(vector<TString> filename,vector<double> kw,vector<double> kwind1, vector<double> kwind2, vector<double> piw,vector<double> piwind1, vector<double> piwind2 ,
		vector<double> prw,vector<double> prwind1, vector<double> prwind2, TH2D *&rkaon, TH2D *&rpion,TH2D *&rproton)
{
			double kmss2=TMath::Sq(0.4937);
			 rkaon = new TH2D("rkaon","Phase space distribution of K^{+};y;pt [GeV/Q]",200,0.5,2.5,200,0,1.4);
		 	double pimss2=TMath::Sq(0.1396);
			 rpion= new TH2D("rpion","Phase space distribution of #pi^{+};y;pt [GeV/Q]",200,0.6,3.6,200,0,1.4);
		 double prmss2=TMath::Sq(0.9383);
		 rproton= new TH2D("rproton","Phase space distribution of p;y;pt [GeV/Q]",200,0.2,1.8,200,0,1.5);
	TChain *data = new TChain("BmnData");
	for (int id = 0 ; id <filename.size() ; id++ ) {
		data->Add(filename[id]);
	}
	vector<double> *px=nullptr, *py=nullptr,*pz=nullptr, *p=nullptr,*m2=nullptr,*hitz=nullptr;
	TBranch *brpx=nullptr, *brpy=nullptr,*brpz=nullptr,*brp=nullptr,*brm2=nullptr,*bhitz=nullptr;

	int ntrack = 0;
        vector<double> *TrackStartX=nullptr, *TrackStartY=nullptr;
        TVector3 *vertex=nullptr;
        double centerOfTargetX = 0.40;
        double centerOfTargetY = 0.15;
        double radiusOfTarget = 1.2;
	data->SetBranchAddress("PrimaryVertexNTracks",&ntrack);
        data->SetBranchAddress("PrimaryVertexPos",&vertex);
        data->SetBranchAddress("TrackStartPointX",&TrackStartX);
        data->SetBranchAddress("TrackStartPointY",&TrackStartY);
	data->SetBranchAddress("TrackPx",&px,&brpx);
        data->SetBranchAddress("TrackPy",&py,&brpy);
	data->SetBranchAddress("TrackPz",&pz,&brpz);
        data->SetBranchAddress("TrackPq",&p,&brp);
        data->SetBranchAddress("TrackTof400Mq2",&m2,&brm2);
	data->SetBranchAddress("Hit400PlaneNumber",&hitz,&bhitz);
	
	double y=0, E = 0,pt = 0;
	
	for (int i = 0; i <data->GetEntries();i++) {
		data->GetEntry(i);
		if (ntrack < 4) continue;
                if (vertex->Z() < -0.5 || vertex->Z() > +0.5) continue;
                if (TMath::Sqrt(TMath::Sq(vertex->X() - centerOfTargetX) + TMath::Sq(vertex->Y() - centerOfTargetY)) > radiusOfTarget) continue;
		
		for (int j=0;j<m2->size();j++) {
			if (TMath::Sqrt(TMath::Sq(TrackStartX->at(j) - vertex->X()) + TMath::Sq(TrackStartY->at(j) - vertex->Y())) > 1.) continue;
			for (int n = 1; n < 22; n++) {
			if ((p->at(j) >= .1*(n-1)+.5) && (p->at(j) <= .1*n+.5)) {
				if ((m2->at(j) >=kwind1[n-1]) && (m2->at(j) <= kwind2[n-1])) {
					E = sqrt(kmss2+p->at(j)*p->at(j));
					pt = sqrt(px->at(j)*px->at(j) + py->at(j)*py->at(j));
					y = (1/2.)*TMath::Log((E+pz->at(j))/(E-pz->at(j)));
					rkaon->Fill(y,pt,kw[n]);
				}
				if ((m2->at(j) >=piwind1[n-1]) && (m2->at(j) <= piwind2[n-1])) {
					E = sqrt(pimss2+p->at(j)*p->at(j));
                                        pt = sqrt(px->at(j)*px->at(j) + py->at(j)*py->at(j));
                                        y = (1/2.)*TMath::Log((E+pz->at(j))/(E-pz->at(j)));
                                        rpion->Fill(y,pt,piw[n]);
				}
				if ((m2->at(j) >=prwind1[n-1]) && (m2->at(j) <= prwind2[n-1])) {
                                        E = sqrt(prmss2+p->at(j)*p->at(j));
                                        pt = sqrt(px->at(j)*px->at(j) + py->at(j)*py->at(j));
                                        y = (1/2.)*TMath::Log((E+pz->at(j))/(E-pz->at(j)));
                                        rproton->Fill(y,pt,prw[n]);
                                }
			}
			}
		}
	}
	
}

void ragidityVsTransverseMment400(vector<TString> filename, TString particleIDfile) {
	vector<double> kw; vector<double> kwind1; vector<double> kwind2; vector<double> piw;vector<double> piwind1; vector<double> piwind2; vector<double> prw;vector<double> prwind1; vector<double> prwind2;
	read(particleIDfile,kwind1,kwind2,kw,piwind1,piwind2,piw,prwind1,prwind2,prw);
	cout<<">>>> "<<kw.size()<<endl;
	TFile *out = new TFile("25/histragidity400.root","recreate");
	TH2D *rk=nullptr,*rpi=nullptr,*rpr=nullptr;
	ragidityVsTransverseMment(filename,kw,kwind1, kwind2, piw,piwind1, piwind2,prw,prwind1,prwind2,rk,rpi,rpr);
	out->WriteObjectAny(rk,"TH2D","rkaon");
	out->WriteObjectAny(rpi,"TH2D","rpion");
	out->WriteObjectAny(rpr,"TH2D","rproton");
TString name[]= {"kaon","pion","proton"};
	double xcheck[100], ycheck[100], mss2[]={TMath::Sq(0.4937),TMath::Sq(0.1396),TMath::Sq(0.9383)};	
	TF1 *z= new TF1("z","(1/2.)*TMath::Log((sqrt([0]*[0]+[1])+sqrt([0]*[0]-x*x))/(sqrt([0]*[0]+[1])-sqrt([0]*[0]-x*x)))",0,1.5);

	for (int i = 0; i < 3 ; i++) {
	z->SetParameter(0,0.5);
	z->SetParameter(1,mss2[i]);
	for (int n = 1; n < 101; n++) {
		xcheck[n-1]=z->Eval(0.015*n);
		ycheck[n-1]=0.015*n;
	}

	TGraph *v= new TGraph(100,xcheck,ycheck);
	v->SetMarkerStyle(20);

	z->SetParameter(0,1.);
        for (int n = 1; n < 101; n++) {
                xcheck[n-1]=z->Eval(0.015*n);
                ycheck[n-1]=0.015*n;
        }

        TGraph *v1= new TGraph(100,xcheck,ycheck);
        v1->SetMarkerStyle(20);
	v1->SetMarkerColor(kOrange);

	z->SetParameter(0,1.5);
        for (int n = 1; n < 101; n++) {
                xcheck[n-1]=z->Eval(0.015*n);
                ycheck[n-1]=0.015*n;
        }

        TGraph *v2= new TGraph(100,xcheck,ycheck);
        v2->SetMarkerStyle(20);
	v2->SetMarkerColor(kGray);

	z->SetParameter(0,2.6);
        for (int n = 1; n < 101; n++) {
                xcheck[n-1]=z->Eval(0.015*n);
                ycheck[n-1]=0.015*n;
        }

        TGraph *v5= new TGraph(100,xcheck,ycheck);
        v5->SetMarkerStyle(20);
        v5->SetMarkerColor(kMagenta);
	out->WriteObjectAny(v,"TGraph",name[i]+"Momentum0GeV5");
	out->WriteObjectAny(v1,"TGraph",name[i]+"Momentum1GeV");
	out->WriteObjectAny(v2,"TGraph",name[i]+"Momentum1GeV5");
	out->WriteObjectAny(v5,"TGraph",name[i]+"Momentum2GeV6");
	}
	
			out->Close();
	
}

