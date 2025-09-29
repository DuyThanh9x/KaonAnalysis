BmnKalmanFilter* fKalman = new BmnKalmanFilter();
FairRunAna* fRunAna = new FairRunAna();

void reduceData (int runPeriod, int runID, int evID, int partID)
{
	
	//gInterpreter->GenerateDictionary("vector<TVector3>","vector;TVector3.h");
	//auto clss=TClass::GetClass("vector<TVector3>");
	//if (!clss || !clss->GetCollectionProxy()) Error("check","no proxy");
	gRandom->SetSeed(0);
       // put proper run number
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
	TTreeReaderValue<BmnTrigInfoDst> trigger(Intree,"BmnTrigInfo.");
	TTreeReaderValue<CbmVertex> vertex(Intree,"MpdVertex.");
	TTreeReaderArray<BmnGlobalTrack> track(Intree,"BmnGlobalTrack");
	TTreeReaderArray<BmnTofHit> tof400Hit(Intree,"BmnTof400Hit");
	TTreeReaderArray<BmnTofHit> tof700Hit(Intree,"BmnTof700Hit");

	const auto nentries = Intree.GetEntries();

	TFile output(Form("ExtractData_run%d_Top_%d_ev%d_p%d.root",runPeriod,runID,evID,partID),"recreate");
	auto NewTree = new TTree("BmnData","BmnData");
	
	int nEvents = 0;
	TVector3 PVertexPos(-100, -100, -100);

	//gSystem->Load("AutoDict_vector_TVector3__cxx.so");
	std::vector<double> TrackStartPointX, TrackStartPointY, TrackStartPointZ, TrackEndPointX, TrackEndPointY, TrackEndPointZ;
	TrackStartPointX.clear(); TrackStartPointY.clear(); TrackStartPointZ.clear(); TrackEndPointX.clear(); TrackEndPointY.clear(); TrackEndPointZ.clear();

	double PVChi2 = -1;

	vector<double> TrackChi2, TrackLength, TrackPq, TrackPx, TrackPy, TrackPz, TrackTof400Mq2, TrackTof700Mq2, Beta400, Beta700;
	TrackChi2.clear(); TrackLength.clear(); TrackPq.clear(); TrackPx.clear(); TrackPy.clear(); TrackPz.clear();
	TrackTof400Mq2.clear(); TrackTof700Mq2.clear(); Beta400.clear(); Beta700.clear();
	
	int PVNTracks = -1, PVNDF = -1;
	
	vector<int> TrackNHits, TrackNDF, TrackTof400HitID, TrackTof700HitID;
	TrackNHits.clear(); TrackNDF.clear(); TrackTof400HitID.clear(); TrackTof700HitID.clear();
	
	vector<double> HitPosX400, HitPosY400, HitPosZ400;
	HitPosX400.clear(); HitPosY400.clear(); HitPosZ400.clear();
	vector<double> HitTimeOfFlight400;
	HitTimeOfFlight400.clear();
	vector<double> HitLength400;
	HitLength400.clear();

	vector<double> HitPosX700, HitPosY700, HitPosZ700;
        HitPosX700.clear(); HitPosY700.clear(); HitPosZ700.clear();
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

	NewTree->Branch("HitPosX400",&HitPosX400);
	NewTree->Branch("HitPosY400",&HitPosY400);
	NewTree->Branch("HitPosZ400",&HitPosZ400);
	NewTree->Branch("HitTimeOfFlight400",&HitTimeOfFlight400);
	NewTree->Branch("HitLength400",&HitLength400);

	NewTree->Branch("HitPosX700",&HitPosX700);
	NewTree->Branch("HitPosY700",&HitPosY700);
	NewTree->Branch("HitPosZ700",&HitPosZ700);
        NewTree->Branch("HitTimeOfFlight700",&HitTimeOfFlight700);
	NewTree->Branch("HitLength700",&HitLength700);
	
	int EventID = 0;
	while (Intree.Next()) {
		if (EventID % 1000 == 0)
			cout<<">>>>>>>>>>>>>>>>Event ID: "<<EventID<<endl;
		//Intree.GetEntry(EventID);
		EventID++;
		//if (!(trigger->IsTriggerBitTrueAR(7)))
			//continue;
	
		if(vertex->GetNTracks() < 4)
			continue;

		if (vertex->GetZ() < -0.5 || vertex->GetZ() > +0.5)
           	        continue;

        	Float_t centerOfTargetX = 0.40;
        	Float_t centerOfTargetY = 0.15;
        	Float_t radiusOfTarget = 1.2;

        	if (TMath::Sqrt(TMath::Sq(vertex->GetX() - centerOfTargetX) + TMath::Sq(vertex->GetY() - centerOfTargetY)) > radiusOfTarget)
            		continue;

		PVertexPos.SetX(vertex->GetX());
		PVertexPos.SetY(vertex->GetY());
		PVertexPos.SetZ(vertex->GetZ());
		PVChi2 = vertex->GetChi2();
		PVNDF = vertex->GetNDF();
		PVNTracks = vertex->GetNTracks();
		for (BmnGlobalTrack& gl : track) {
			//int n = track.GetEntriesFast();
			
			Int_t pdg = (gl.GetParamFirst()->GetQp() > 0) ? 2212 : -211;
        		if (fKalman->TGeoTrackPropagate(gl.GetParamFirst(), vertex->GetZ(), pdg, nullptr, nullptr, kTRUE) == kBMNERROR)
               			 continue;
            		Float_t dtrack = TMath::Sqrt(TMath::Sq(gl.GetParamFirst()->GetX() - vertex->GetX()) + TMath::Sq(gl.GetParamFirst()->GetY() - vertex->GetY()));
            		if (dtrack > 1.0)
                		continue;

			TrackStartPointX.push_back(gl.GetParamFirst()->GetX());
			TrackStartPointY.push_back(gl.GetParamFirst()->GetY());
			TrackStartPointZ.push_back(gl.GetParamFirst()->GetZ());
			TrackEndPointX.push_back(gl.GetParamLast()->GetX());
			TrackEndPointY.push_back(gl.GetParamLast()->GetY());
			TrackEndPointZ.push_back(gl.GetParamLast()->GetZ());

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
		}

		for (BmnTofHit& hit : tof400Hit) {
			HitPosX400.push_back(hit.GetX());
			HitPosY400.push_back(hit.GetY());
			HitPosZ400.push_back(hit.GetZ());
			HitTimeOfFlight400.push_back(hit.GetTimeStamp());
			HitLength400.push_back(hit.GetLength());
		}

		for (BmnTofHit& hit2 : tof700Hit) {
                        HitPosX700.push_back(hit2.GetX());
			HitPosY700.push_back(hit2.GetY());
			HitPosZ700.push_back(hit2.GetZ());
                        HitTimeOfFlight700.push_back(hit2.GetTimeStamp());
			HitLength700.push_back(hit2.GetLength());
                }

		NewTree->Fill();
		nEvents++;

		TrackStartPointX.clear(); TrackStartPointY.clear(); TrackStartPointZ.clear(); TrackEndPointX.clear(); TrackEndPointY.clear(); TrackEndPointZ.clear();
		TrackChi2.clear(); TrackLength.clear(); TrackPq.clear(); TrackTof400Mq2.clear(); TrackTof700Mq2.clear(); Beta400.clear(); Beta700.clear();
        	TrackPx.clear(); TrackPy.clear(); TrackPz.clear();
		TrackNHits.clear(); TrackNDF.clear(); TrackTof400HitID.clear(); TrackTof700HitID.clear();
		
		HitPosX400.clear(); HitPosY400.clear(); HitPosZ400.clear();
	        HitTimeOfFlight400.clear();
		HitLength400.clear();
	        
		HitPosX700.clear(); HitPosY700.clear(); HitPosZ700.clear();
	        HitTimeOfFlight700.clear();
		HitLength700.clear();
		
		PVertexPos.SetXYZ(-100, -100, -100);
	}
	

	NewTree->Print();
	output.Write();
	output.Close();
}
