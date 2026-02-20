void SaveHistMass2ID700 (vector<TString> filename) {
	TChain *data = new TChain("BmnData");
        for (int id = 0 ; id <filename.size() ; id++ ) {
                data->Add(filename[id]);
        }

	int ntrack = 0,VertexNDF = -1;
        vector<double> *TrackStartX=nullptr, *TrackStartY=nullptr,*m2=nullptr, *p=nullptr;
        TVector3 *vertex=nullptr;
        double centerOfTargetX = 0.40;
        double centerOfTargetY = 0.15;
        double radiusOfTarget = 1.2;double VertexChi2 = -1;
	TH2D *masv = new TH2D("masv","",250,-.5,2,21,0.5,2.6);
	data->SetBranchAddress("PrimaryVertexNTracks",&ntrack);
        data->SetBranchAddress("PrimaryVertexPos",&vertex);
	data->SetBranchAddress("PrimaryVertexChi2",&VertexChi2);
        data->SetBranchAddress("PrimaryVertexNDF",&VertexNDF);
        data->SetBranchAddress("TrackStartPointX",&TrackStartX);
        data->SetBranchAddress("TrackStartPointY",&TrackStartY);
	data->SetBranchAddress("TrackTof700Mq2",&m2);
	data->SetBranchAddress("TrackPq",&p);
	for (int i = 0; i <data->GetEntries();i++) {
                data->GetEntry(i);
                if (ntrack < 4) continue;
                if (vertex->Z() < -0.26 || vertex->Z() > +0.2) continue;
                if (TMath::Sqrt(TMath::Sq(vertex->X() - centerOfTargetX) + TMath::Sq(vertex->Y() - centerOfTargetY)) > radiusOfTarget) continue;
		double chi2ndf = VertexChi2/VertexNDF;
                if ((chi2ndf < 0.1) || (chi2ndf >10)) continue;
		for (int j=0;j<m2->size();j++) {
                        if (TMath::Sqrt(TMath::Sq(TrackStartX->at(j) - vertex->X()) + TMath::Sq(TrackStartY->at(j) - vertex->Y())) > 1.) continue;
			if ((m2->at(j) >= -0.5) && (m2->at(j) <= 2) && (p->at(j) >= 0.5) && (p->at(j) <= 2.6) ) {
				masv->Fill(m2->at(j),p->at(j));
			}
		}
	}
	TFile *out = new TFile("25/histmass700.root","recreate");masv->Write();
	for (int i = 1; i < 22; i++) {
                masv->ProjectionX(Form("mass%d",i),i,i);
                TH1D *mass1d = (TH1D*)gDirectory->Get(Form("mass%d",i));
		mass1d->SetTitle(Form("Mq2 distribution at TOF700 with momentum range [%.2f, %.2f] GeV;GeV^{2};#Candidates",.1*(i-1)+.5,.1*i+.5));
		mass1d->Write();
	}
        out->Close();
}

