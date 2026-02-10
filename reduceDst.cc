FairRunAna* fRunAna = new FairRunAna();
#define TIME_CUT 66.0           // ns - upper bound (negative)
#define TIME_MIN -1766.0        // ns - minimum time  
#define BIN_WIDTH 50.0          // ns - bin width
#define NUM_BINS 34             // (1766 - 66) / 50 = 34 bins
#define PILEUP_THRESHOLD 10     // hits per bin threshold for pileup

bool isPileup(TTreeReaderArray<BmnTofHit> &TofHit) {
        int nHits = TofHit.GetSize();
 	int maxHitsTOF = 0;
        int tof_bins[NUM_BINS] = {0};
        for (int i = 0; i < nHits; i++) {
                BmnTofHit hit = TofHit.At(i);
                double time = hit.GetTimeStamp();
                if (time >= TIME_MIN && time <= -TIME_CUT) {
                    int bin = (int)((time - TIME_MIN) / BIN_WIDTH);
                    if (bin >= 0 && bin < NUM_BINS) {
                        tof_bins[bin]++;
                    }
                }
            }
        for (int bin = 0; bin < NUM_BINS; bin++) {
                if (tof_bins[bin] > maxHitsTOF) {
                    maxHitsTOF = tof_bins[bin];
                }
            }
        if (maxHitsTOF >= PILEUP_THRESHOLD) {
                 return true;
        } else return false;
}

void reduceData (int runPeriod, int runID, int evID, int partID)
{
	gRandom->SetSeed(0);
        TString geoFileName = Form("current_geo_file_%d.root", UInt_t(gRandom->Integer(UINT32_MAX)));
        Int_t res_code = UniRun::ReadGeometryFile(runPeriod, runID, (char*)geoFileName.Data());
        TGeoManager::Import(geoFileName);
       
        BmnFieldMap* magField = new BmnNewFieldMap("FieldMap_1900_extrap_noPed.root");
       
        UniRun* pCurrentRun = UniRun::GetRun(runPeriod, runID);
        Double_t* field_voltage = pCurrentRun->GetFieldVoltage();
        Double_t fieldScale = (*field_voltage) / 112.0;
        magField->SetScale(fieldScale);
        fRunAna->SetField(magField);
        magField->Init();

	TString targ;
        if (pCurrentRun->GetTargetParticle() == NULL) {
            targ = "-";
        } else {
            targ = (pCurrentRun->GetTargetParticle())[0];
        }
        TString beam = pCurrentRun->GetBeamParticle();
	cout<<">>>>>>>>>>>>>>>>>>>>>>>>Analyse data<<<<<<<<<<<<<<<<<<<"<<endl;
	cout<<">>>>>>>>>> Target: "<<targ<<endl
	    <<">>>>>>>>>> Beam:   "<<beam<<endl;

	TFile *input = TFile::Open(Form("mpd_run_Top_%d_ev%d_p%d.root",runID,evID,partID));
	TTreeReader Intree("bmndata",input);
	TTreeReaderValue<CbmVertex> vertex(Intree,"MpdVertex.");
	TTreeReaderValue<BmnTrigInfoDst> trigger(Intree, "BmnTrigInfo.");
	TTreeReaderArray<BmnGlobalTrack> track(Intree,"BmnGlobalTrack");
	TTreeReaderArray<CbmStsTrack> stsvectr(Intree,"StsVector");
	TTreeReaderArray<BmnTofHit> tof400Hit(Intree,"BmnTof400Hit");
	TTreeReaderArray<BmnTofHit> tof700Hit(Intree,"BmnTof700Hit");

	const auto nentries = Intree.GetEntries();

	TFile output(Form("ExtractData_run%d_Top_%d_ev%d_p%d.root",runPeriod,runID,evID,partID),"recreate");
	auto NewTree = new TTree("BmnData","BmnData");
	
	int nEvents = 0;
	TVector3 PVertexPos(-100, -100, -100);

	std::vector<double> TrackStartPointX, TrackStartPointY, TrackStartPointZ, TrackEndPointX, TrackEndPointY, TrackEndPointZ;
	TrackStartPointX.clear(); TrackStartPointY.clear(); TrackStartPointZ.clear(); TrackEndPointX.clear(); TrackEndPointY.clear(); TrackEndPointZ.clear();

	double PVChi2 = -1;
	bool isPileup400 = false, isPileup700  = false;
	vector<double> TrackChi2, TrackLength, TrackPq, TrackPx, TrackPy, TrackPz, TrackTof400Mq2, TrackTof700Mq2, Beta400, Beta700;
	TrackChi2.clear(); TrackLength.clear(); TrackPq.clear(); TrackPx.clear(); TrackPy.clear(); TrackPz.clear();
	TrackTof400Mq2.clear(); TrackTof700Mq2.clear(); Beta400.clear(); Beta700.clear();
	vector<int> TrackID;
	TrackID.clear();

	vector<double> MLtrack;
	MLtrack.clear();
	vector<int> StsVectrID;
        StsVectrID.clear();

	int PVNTracks = -1, PVNDF = -1;
	
	vector<int> TrackNHits, TrackNDF, TrackTof400HitID, TrackTof700HitID;
	TrackNHits.clear(); TrackNDF.clear(); TrackTof400HitID.clear(); TrackTof700HitID.clear();
	
	vector<double> HitPosX400, HitPosY400, HitPosZ400;
	HitPosX400.clear(); HitPosY400.clear(); HitPosZ400.clear();
	
	vector<int> Hit400PlaneNumber, Hit400StripNumber;
	Hit400PlaneNumber.clear();Hit400StripNumber.clear();

	vector<double> HitTimeOfFlight400;
	HitTimeOfFlight400.clear();
	vector<double> HitLength400;
	HitLength400.clear();

	vector<double> HitPosX700, HitPosY700, HitPosZ700;
        HitPosX700.clear(); HitPosY700.clear(); HitPosZ700.clear();
        vector<int> Hit700PlaneNumber, Hit700StripNumber;
        Hit700PlaneNumber.clear();Hit700StripNumber.clear();
	vector<double> HitTimeOfFlight700;
        HitTimeOfFlight700.clear();
	vector<double> HitLength700;
	HitLength700.clear();
	
	NewTree->Branch("PrimaryVertexPos",&PVertexPos);
	NewTree->Branch("PrimaryVertexChi2",&PVChi2);
	NewTree->Branch("PrimaryVertexNDF",&PVNDF);
	NewTree->Branch("PrimaryVertexNTracks",&PVNTracks);

	NewTree->Branch("TrackStartPointX",&TrackStartPointX);
	NewTree->Branch("TrackStartPointY",&TrackStartPointY);
	NewTree->Branch("TrackStartPointZ",&TrackStartPointZ);

	NewTree->Branch("TrackEndPointX",&TrackEndPointX);
	NewTree->Branch("TrackEndPointY",&TrackEndPointY);
	NewTree->Branch("TrackEndPointZ",&TrackEndPointZ);

	NewTree->Branch("TrackID",&TrackID);
	NewTree->Branch("TrackTof400HitID",&TrackTof400HitID);
	NewTree->Branch("TrackTof700HitID",&TrackTof700HitID);
	NewTree->Branch("TrackPq",&TrackPq);
	NewTree->Branch("TrackPx",&TrackPx);
	NewTree->Branch("TrackPy",&TrackPy);
	NewTree->Branch("TrackPz",&TrackPz);
	NewTree->Branch("TrackTof400Mq2",&TrackTof400Mq2);
	NewTree->Branch("TrackTof700Mq2",&TrackTof700Mq2);
	NewTree->Branch("Beta400",&Beta400);
	NewTree->Branch("Beta700",&Beta700);
	NewTree->Branch("TrackNHits",&TrackNHits);
	NewTree->Branch("TrackNDF",&TrackNDF);
	NewTree->Branch("TrackChi2",&TrackChi2);
	NewTree->Branch("TrackLength",&TrackLength);

	NewTree->Branch("TMVATrack",&MLtrack);

	NewTree->Branch("HitPosX400",&HitPosX400);
	NewTree->Branch("HitPosY400",&HitPosY400);
	NewTree->Branch("HitPosZ400",&HitPosZ400);
	NewTree->Branch("HitTimeOfFlight400",&HitTimeOfFlight400);
	NewTree->Branch("HitLength400",&HitLength400);
	NewTree->Branch("Hit400PlaneNumber",&Hit400PlaneNumber);
	NewTree->Branch("Hit400StripNumber",&Hit400StripNumber);
	NewTree->Branch("isPileup400",&isPileup400);

	NewTree->Branch("HitPosX700",&HitPosX700);
	NewTree->Branch("HitPosY700",&HitPosY700);
	NewTree->Branch("HitPosZ700",&HitPosZ700);
        NewTree->Branch("HitTimeOfFlight700",&HitTimeOfFlight700);
	NewTree->Branch("HitLength700",&HitLength700);
	NewTree->Branch("Hit700PlaneNumber",&Hit700PlaneNumber);
        NewTree->Branch("Hit700StripNumber",&Hit700StripNumber);
	NewTree->Branch("isPileup700",&isPileup700);
	
	int EventID = 0;
	while (Intree.Next()) {
		if (EventID % 1000 == 0)
			cout<<">>>>>>>>>>>>>>>>Event ID: "<<EventID<<endl;
		EventID++;
		/*if(vertex->GetNTracks() < 2)
			continue;

		if (vertex->GetZ() < -0.5 || vertex->GetZ() > +0.5)
           	        continue;

        	Float_t centerOfTargetX = 0.40;
        	Float_t centerOfTargetY = 0.15;
        	Float_t radiusOfTarget = 1.2;

        	if (TMath::Sqrt(TMath::Sq(vertex->GetX() - centerOfTargetX) + TMath::Sq(vertex->GetY() - centerOfTargetY)) > radiusOfTarget)
            		continue;*/
		if (!(trigger->IsTriggerBitTrueAR(7)))   // bit #7 is for CCT2
            		continue;
		if (isPileup(tof400Hit)) isPileup400 = true; else isPileup400 = false;
		if (isPileup(tof700Hit)) isPileup700 = true; else isPileup700  = false; //continue;
		PVertexPos.SetX(vertex->GetX());
		PVertexPos.SetY(vertex->GetY());
		PVertexPos.SetZ(vertex->GetZ());
		PVChi2 = vertex->GetChi2();
		PVNDF = vertex->GetNDF();
		PVNTracks = vertex->GetNTracks();
		for (BmnGlobalTrack& gl : track) {
			//int n = track.GetEntriesFast();
			
			/*Int_t pdg = (gl.GetParamFirst()->GetQp() > 0) ? 2212 : -211;
        		if (fKalman->TGeoTrackPropagate(gl.GetParamFirst(), vertex->GetZ(), pdg, nullptr, nullptr, kTRUE) == kBMNERROR)
               			 continue;
            		Float_t dtrack = TMath::Sqrt(TMath::Sq(gl.GetParamFirst()->GetX() - vertex->GetX()) + TMath::Sq(gl.GetParamFirst()->GetY() - vertex->GetY()));
            		if (dtrack > 1.0)
                		continue;*/

			TrackStartPointX.push_back(gl.GetParamFirst()->GetX());
			TrackStartPointY.push_back(gl.GetParamFirst()->GetY());
			TrackStartPointZ.push_back(gl.GetParamFirst()->GetZ());
			TrackEndPointX.push_back(gl.GetParamLast()->GetX());
			TrackEndPointY.push_back(gl.GetParamLast()->GetY());
			TrackEndPointZ.push_back(gl.GetParamLast()->GetZ());

			TrackID.push_back(gl.GetGemTrackIndex());
			CbmStsTrack &ststrack = stsvectr.At(gl.GetGemTrackIndex());
			MLtrack.push_back(ststrack.GetB());
			TrackTof400HitID.push_back(gl.GetTof1HitIndex());
			TrackTof700HitID.push_back(gl.GetTof2HitIndex());
			TrackPq.push_back(gl.GetP());
			double pzz = abs(gl.GetP())/sqrt(gl.GetParamFirst()->GetTx() *gl.GetParamFirst()->GetTx() + gl.GetParamFirst()->GetTy() * gl.GetParamFirst()->GetTy() + 1);
			TrackPz.push_back(pzz);
			TrackPx.push_back(gl.GetParamFirst()->GetTx() * pzz);
			TrackPy.push_back(gl.GetParamFirst()->GetTy() * pzz);
			TrackTof400Mq2.push_back(gl.GetMass2(1));
			TrackTof700Mq2.push_back(gl.GetMass2(2));
			Beta400.push_back(gl.GetBeta(1));
			Beta700.push_back(gl.GetBeta(2));
			TrackNHits.push_back(gl.GetNHits());
			TrackChi2.push_back(gl.GetChi2());
			TrackLength.push_back(gl.GetLength());
			TrackNDF.push_back(gl.GetNDF());

			if (gl.GetTof1HitIndex() != -1) {
				BmnTofHit &hit = tof400Hit.At(gl.GetTof1HitIndex());
				HitPosX400.push_back(hit.GetX());
				HitPosY400.push_back(hit.GetY());
				HitPosZ400.push_back(hit.GetZ());
				HitTimeOfFlight400.push_back(hit.GetTimeStamp());
				HitLength400.push_back(hit.GetLength());
				Hit400PlaneNumber.push_back(hit.GetModule());
				Hit400StripNumber.push_back(hit.GetStation());
				}
				else {
					HitPosX400.push_back(-10000);
					HitPosY400.push_back(-10000);
					HitPosZ400.push_back(-10000);
					HitLength400.push_back(-10000);
					HitTimeOfFlight400.push_back(-10000);
					Hit400PlaneNumber.push_back(-1);
					Hit400StripNumber.push_back(-1);
			
				}
			if (gl.GetTof2HitIndex() != -1) {
				BmnTofHit &hit2 = tof700Hit.At(gl.GetTof2HitIndex());
					HitPosX700.push_back(hit2.GetX());
					HitPosY700.push_back(hit2.GetY());
					HitPosZ700.push_back(hit2.GetZ());
                        		HitTimeOfFlight700.push_back(hit2.GetTimeStamp());
					HitLength700.push_back(hit2.GetLength());
               
					Hit700PlaneNumber.push_back(hit2.GetModule());
                        		Hit700StripNumber.push_back(hit2.GetStation());
				}
				else {
                                        HitPosX700.push_back(-10000);
                                        HitPosY700.push_back(-10000);
                                        HitPosZ700.push_back(-10000);
					HitLength700.push_back(-10000);
                                        HitTimeOfFlight700.push_back(-10000);
                                        Hit700PlaneNumber.push_back(-1);
                                        Hit700StripNumber.push_back(-1);

                        	}

		}

		NewTree->Fill();
		nEvents++;

		TrackStartPointX.clear(); TrackStartPointY.clear(); TrackStartPointZ.clear(); TrackEndPointX.clear(); TrackEndPointY.clear(); TrackEndPointZ.clear();
		TrackChi2.clear(); TrackLength.clear(); TrackPq.clear(); TrackTof400Mq2.clear(); TrackTof700Mq2.clear(); Beta400.clear(); Beta700.clear();
        	TrackPx.clear(); TrackPy.clear(); TrackPz.clear();
		TrackNHits.clear(); TrackNDF.clear(); TrackTof400HitID.clear(); TrackTof700HitID.clear();

		TrackID.clear();
		
		StsVectrID.clear();
		MLtrack.clear();

		HitPosX400.clear(); HitPosY400.clear(); HitPosZ400.clear();
	        HitTimeOfFlight400.clear();
		HitLength400.clear();
	        
		HitPosX700.clear(); HitPosY700.clear(); HitPosZ700.clear();
	        HitTimeOfFlight700.clear();
		HitLength700.clear();
		
		Hit400PlaneNumber.clear();Hit400StripNumber.clear();
		Hit700PlaneNumber.clear();Hit700StripNumber.clear();
		PVertexPos.SetXYZ(-100, -100, -100);
	}
	

	NewTree->Print();
	output.Write();
	output.Close();
}
