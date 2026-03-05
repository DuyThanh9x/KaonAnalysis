auto GetCentrBin_8120_8170 = [](Double_t _refMult){
//RunId_corr_factor_h2_RunId_nTracks_8120_8170
        if(_refMult < 6) return -1;
        else if ( _refMult < 12 ) return 7; // 70-80%
        else if ( _refMult < 22 ) return 6; // 60-70%
        else if ( _refMult < 33 ) return 5; // 50-60%
        else if ( _refMult < 49 ) return 4; // 40-50%
        else if ( _refMult < 71 ) return 3; // 30-40%
        else if ( _refMult < 99 ) return 2; // 20-30%
        else if ( _refMult < 137) return 1; // 10-20%
        else if ( _refMult < 236) return 0; //  0-10%
        else if ( _refMult >=236) return -1;
        return -1;
};
void TrackCorrec(vector<TString> filename, vector<double> runlist){
        TChain *data = new TChain("BmnData");
        for (int id = 0 ; id <filename.size() ; id++ ) {
                data->Add(filename[id]);
        }
        auto Cr = TFile::Open("CorrRunId.root");
        TGraphErrors *corr = (TGraphErrors*)Cr->Get("RunId_corr_factor_h2_RunId_nTracks_8120_8170");
        double nTrackCorr,centrality;
	int ntrackVertex = 0,VertexNDF = -1,n = 0;
        vector<double> *TrackStartX=nullptr, *TrackStartY=nullptr;
        TVector3 *vertex=nullptr;
        double centerOfTargetX = 0.40;
        double centerOfTargetY = 0.15;
        double radiusOfTarget = 1.2;
        double VertexChi2 = -1;
        TH1D *nTrck = new TH1D("nTrck","",50,0,236);
        data->SetBranchAddress("PrimaryVertexNTracks",&ntrackVertex);
        data->SetBranchAddress("PrimaryVertexPos",&vertex);
        data->SetBranchAddress("PrimaryVertexChi2",&VertexChi2);
        data->SetBranchAddress("PrimaryVertexNDF",&VertexNDF);
        data->SetBranchAddress("TrackStartPointX",&TrackStartX);
        data->SetBranchAddress("TrackStartPointY",&TrackStartY);
        for (int i = 0; i <data->GetEntries();i++) {
                data->GetEntry(i);
                if (ntrackVertex < 4) continue;
                if (vertex->Z() < -0.26 || vertex->Z() > +0.2) continue;
                if (TMath::Sqrt(TMath::Sq(vertex->X() - centerOfTargetX) + TMath::Sq(vertex->Y() - centerOfTargetY)) > radiusOfTarget) continue;
                double chi2ndf = VertexChi2/VertexNDF;
                if ((chi2ndf < 0.1) || (chi2ndf >10)) continue;
                double count = 0;
                for (int j=0;j<TrackStartX->size();j++) {
                        if (TMath::Sqrt(TMath::Sq(TrackStartX->at(j) - vertex->X()) + TMath::Sq(TrackStartY->at(j) - vertex->Y())) > 1.) continue;
                        count++;
                }
		if(data->GetCurrentFile()->GetName() == filename[n])
                        nTrackCorr =count * corr->Eval(runlist[n]);
                else {
                        n++;nTrackCorr =count * corr->Eval(runlist[n]);
                }
                centrality = GetCentrBin_8120_8170(nTrackCorr);
		nTrck->Fill(nTrackCorr);
        }
	nTrck->Draw();
}

