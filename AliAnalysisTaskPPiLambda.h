#ifndef ALIANALYSISTASKPPILAMBDA_H
#define ALIANALYSISTASKPPILAMBDA_H

// ROOT includes
#include <TList.h>
#include <TH1.h>

// AliRoot includes
#include <AliAnalysisTaskSE.h>
#include <AliAODEvent.h>
#include <AliAnalysisUtils.h>
#include <AliPIDResponse.h>
#include <AliAODv0.h>
#include <AliAODcascade.h>
#include <AliEventPoolManager.h>

class AliAnalysisTaskPPiLambda : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskPPiLambda();
  AliAnalysisTaskPPiLambda(const char *name);
  virtual ~AliAnalysisTaskPPiLambda();
  
  virtual void  SetTrigger(UInt_t ktriggerInt=AliVEvent::kINT7){ ftrigBit=ktriggerInt; }
  virtual void  UserCreateOutputObjects();
  virtual void  UserExec(Option_t *option);

  Bool_t   EventSelection(AliAODEvent *data);
  Bool_t   ProtonSelection(AliAODTrack *trk);
  Double_t InvMassLambda(AliAODcascade *casc);
  Double_t InvMassAntiLambda(AliAODcascade *casc);
  Double_t InvMassXi(AliAODcascade *casc);
  Double_t InvMassOmega(AliAODcascade *casc);
  Double_t LambdaCosPointingAngle(AliAODcascade *casc,const Double_t *DecayVtx,const Float_t *point) const;
  Double_t DecayLengthXY(const Double_t *DecayVtx,const Float_t *point) const;
  Double_t xiDecayLengthXY(const Double_t *xiDecayVtx,const Float_t *point) const;
  Double_t InvMasslambda(AliAODv0 *v0);
  Double_t InvMassK0(AliAODv0 *v0);
  Double_t InvMassAntilambda(AliAODv0 *v0);
  Double_t CalculateInvMassAntilambda(AliAODv0 *v0,AliAODTrack *antiprotontrk,AliAODTrack *piontrk);
  Double_t CosPointingAngle(AliAODv0 *v0,const Double_t *DecayVtx,const Float_t *point) const;
  Double_t OpenAngle(Double_t px1,Double_t py1,Double_t pz1,
		     Double_t px2,Double_t py2,Double_t pz2);
  Double_t InvariantMass(Double_t px1,Double_t py1,Double_t pz1,
			 Double_t px2,Double_t py2,Double_t pz2,Double_t energysum);

 private:
  AliAODEvent       *fAOD;    //! AOD object
  AliAODHeader      *fHeader; //! AOD header
  AliAnalysisUtils  *fUtils;
  AliPIDResponse    *fPIDResponse;

  UInt_t            ftrigBit;

  Int_t             nevt;
  Int_t             lambdacounter;
  Int_t             eventdepth[5][10];
  Int_t             eventdepth2[5][10];
  Int_t             mnLambda[5][10][5];
  Double_t          mEnergy[5][10][5][10];
  Double_t          mPx[5][10][5][10];
  Double_t          mPy[5][10][5][10];
  Double_t          mPz[5][10][5][10];  

  AliEventPoolManager *fPoolManagerlambda;
  AliEventPoolManager *fPoolManagerantilambda;
  AliEventPoolManager *fPoolManagerxi;
  AliEventPoolManager *fPoolManagerxip;
  AliEventPoolManager *fPoolManagerproton;
  AliEventPoolManager *fPoolManagerantiproton;
  AliEventPoolManager *fPoolManagerlam;
  AliEventPoolManager *fPoolManagerantilam;
  AliEventPoolManager *fPoolManagerx;
  AliEventPoolManager *fPoolManagerxp;

  TList             *fOutputList; //! Output list
  
  TH1F              *fEvtCounter;
  TH1F              *fEvtPassCut;
  TH1F              *fEvtVtxX;
  TH1F              *fEvtVtxY;
  TH1F              *fEvtVtxZ;
  TH1F              *fEvtVtxTrk;
  TH1F              *fSPDVtxZ;
  TH1F              *fSPDVtxTrk;
  TH2F              *fSPDVtxCor;
  TH1F              *fSPDVtxDisp;
  //multiplicity
  TH2F              *fmultiplicity;
  TH1F              *fmultiplicity_lambda;
  TH1F              *fmultiplicity_antilambda;
  TH1F              *fmultiplicity_xim;
  TH1F              *fmultiplicity_xip;
  TH1F              *fmultiplicity_proton;
  TH1F              *fmultiplicity_antiproton;
  TH1F              *fmultiplicity_lambdalambda;
  TH1F              *fmultiplicity_antilambdaantilambda;
  TH1F              *fmultiplicity_protonxi;
  TH1F              *fmultiplicity_antiprotonxi;
  TH1F              *fmultiplicity_lambdaxi;
  TH1F              *fmultiplicity_antilambdaxi;
  //proton
  TH1F              *fTrkEta;
  TH1F              *fTrkPt;
  TH1F              *fTPCclus;
  TH1F              *fTPCcrsR;
  TH1F              *fTPCclusF;
  TH1F              *fTrkDCAxy;
  TH1F              *fTrkDCAz;  
  TH2F              *fPvsTPCsignal;
  TH2F              *fPvsTOFsignal;
  //proton select
  TH1F              *fPt_proton;
  TH1F              *fPt_antiproton;
  TH2F              *PID_cut;
  //v0 selection
  TH2F              *fPtv0;
  TH2F              *fPtv0_lambda;
  TH1F              *fTransRadius;
  TH1F              *fDCAdaugTov0Vtx;
  TH1F              *fCPA;
  TH1F              *hInvMassLambda;  
  TH1F              *fPtv0_antilambda;
  TH2F              *fPtv0_antilambda_mass_masscut;
  TH1F              *hInvMassAntiLambda;
  //invariant mass lambda
  TH1F              *hInvMassXidecayLambda;
  TH1F              *hInvMassXidecayantiLambda;
  //invariant mass xi
  TH1F              *hInvMassXi;
  TH1F              *hInvMassposXi;
  TH1F              *XiCPA;
  TH1F              *XiTranverseradius;
  TH1F              *dcaXiDaughters;
  TH1F              *fPtcascade;
  TH1F              *fPtv0casc;
  TH1F              *fPtBachelor;
  TH1F              *fPhiXidecayLambda;
  TH1F              *fPtcascade_xip;
  TH1F              *fPtv0casc_xip;
  TH1F              *fPtBachelor_xip;
  TH1F              *fPhiXidecayLambda_xip;
  //invariant mass dibaryon
  TH2F              *hInvMassProtonXi;
  TH2F              *hInvMassantiProtonXi;  //anti proton+xi+
  TH2F              *hInvMassLambdaXi;
  TH2F              *hInvMassAntiLambdaXip;
  //decay rejection single
  TH1F              *fPt_allLambda;
  TH1F              *fPt_xidecaylambda;
  TH1F              *fPt_allAntiLambda;
  TH1F              *fPt_xidecayAntiLambda;
  //decay rejection pair
  TH2F              *hInvMassLambdaLambda_all;
  TH2F              *hInvMassLambdaLambda_onlyprompt;
  TH2F              *hInvMassAntiLambdaAntiLambda_all;
  TH2F              *hInvMassAntiLambdaAntiLambda_onlyprompt;
  //eventmixing
  //TH2F              *hInvMassLambdaLambda_combinatorial;
  TH2F              *hInvMassLambdaLambda_evtpool;
  TH2F              *hInvMassAntiLambdaAntiLambda_evtpool;
  TH2F              *hInvMassProtonXi_evtpool;
  TH2F              *hInvMassAntiProtonXi_evtpool;
  TH2F              *hInvMassLambdaXi_evtpool;
  TH2F              *hInvMassAntiLambdaXi_evtpool;
  //vs relative momentum
  TH2F              *hInvMassLambdaLambda_relmomentum;
  TH2F              *hInvMassLambdaLambda_openangle;
  TH2F              *hInvMassAntiLambdaAntiLambda_relmomentum;
  TH2F              *hInvMassAntiLambdaAntiLambda_openangle;
  TH2F              *hInvMassProtonXi_relmomentum;
  TH2F              *hInvMassProtonXi_openangle;
  TH2F              *hInvMassantiProtonXi_relmomentum;
  TH2F              *hInvMassantiProtonXi_openangle;
  TH2F              *hInvMassLambdaXi_relmomentum;
  TH2F              *hInvMassLambdaXi_openangle;
  TH2F              *hInvMassAntiLambdaXi_relmomentum;
  TH2F              *hInvMassAntiLambdaXi_openangle;
  TH2F              *hAngle_relmom;
  //pt correlation
  TH2F              *hPtcorrelation_lamlam;
  TH2F              *hPtcorrelation_pxi;
  TH2F              *hPtcorrelation_lamxi;
  //lowmass
  TH2F              *hInvMassProtonXi_lowmassx;
  TH2F              *hInvMassProtonXi_lowmassy;
  TH2F              *hInvMassProtonXi_lowmassz;
  TH2F              *hInvMassantiProtonXi_lowmassx;
  TH2F              *hInvMassantiProtonXi_lowmassy;
  TH2F              *hInvMassantiProtonXi_lowmassz;
  //FG LambdaLambda
  TH1F              *hInvMassLambdaLambda_1;
  TH1F              *hInvMassLambdaLambda_2;
  TH1F              *hInvMassLambdaLambda_3;
  TH1F              *hInvMassLambdaLambda_4;
  TH1F              *hInvMassLambdaLambda_5;
  TH1F              *hInvMassLambdaLambda_6;
  TH1F              *hInvMassLambdaLambda_7;
  TH1F              *hInvMassLambdaLambda_8;
  TH1F              *hInvMassLambdaLambda_9;
  TH1F              *hInvMassLambdaLambda_10;
  TH1F              *hInvMassLambdaLambda_11;
  TH1F              *hInvMassLambdaLambda_12;
  TH1F              *hInvMassLambdaLambda_13;
  TH1F              *hInvMassLambdaLambda_14;
  TH1F              *hInvMassLambdaLambda_15;
  TH1F              *hInvMassLambdaLambda_16;
  TH1F              *hInvMassLambdaLambda_17;
  TH1F              *hInvMassAntiLambdaAntiLambda_1;
  TH1F              *hInvMassAntiLambdaAntiLambda_2;
  TH1F              *hInvMassAntiLambdaAntiLambda_3;
  TH1F              *hInvMassAntiLambdaAntiLambda_4;
  TH1F              *hInvMassAntiLambdaAntiLambda_5;
  TH1F              *hInvMassAntiLambdaAntiLambda_6;
  TH1F              *hInvMassAntiLambdaAntiLambda_7;
  TH1F              *hInvMassAntiLambdaAntiLambda_8;
  TH1F              *hInvMassAntiLambdaAntiLambda_9;
  TH1F              *hInvMassAntiLambdaAntiLambda_10;
  TH1F              *hInvMassAntiLambdaAntiLambda_11;
  TH1F              *hInvMassAntiLambdaAntiLambda_12;
  TH1F              *hInvMassAntiLambdaAntiLambda_13;
  TH1F              *hInvMassAntiLambdaAntiLambda_14;
  TH1F              *hInvMassAntiLambdaAntiLambda_15;
  TH1F              *hInvMassAntiLambdaAntiLambda_16;
  TH1F              *hInvMassAntiLambdaAntiLambda_17;
  //FG ProtonXi
  TH1F              *hInvMassProtonXi_1;
  TH1F              *hInvMassProtonXi_2;
  TH1F              *hInvMassProtonXi_3;
  TH1F              *hInvMassProtonXi_4;
  TH1F              *hInvMassProtonXi_5;
  TH1F              *hInvMassProtonXi_6;
  TH1F              *hInvMassProtonXi_7;
  TH1F              *hInvMassProtonXi_8;
  TH1F              *hInvMassProtonXi_9;
  TH1F              *hInvMassProtonXi_10;
  TH1F              *hInvMassProtonXi_11;
  TH1F              *hInvMassProtonXi_12;
  TH1F              *hInvMassProtonXi_13;
  TH1F              *hInvMassProtonXi_14;
  TH1F              *hInvMassProtonXi_15;
  TH1F              *hInvMassProtonXi_16;
  TH1F              *hInvMassProtonXi_17;
  TH1F              *hInvMassAntiProtonXi_1;
  TH1F              *hInvMassAntiProtonXi_2;
  TH1F              *hInvMassAntiProtonXi_3;
  TH1F              *hInvMassAntiProtonXi_4;
  TH1F              *hInvMassAntiProtonXi_5;
  TH1F              *hInvMassAntiProtonXi_6;
  TH1F              *hInvMassAntiProtonXi_7;
  TH1F              *hInvMassAntiProtonXi_8;
  TH1F              *hInvMassAntiProtonXi_9;
  TH1F              *hInvMassAntiProtonXi_10;
  TH1F              *hInvMassAntiProtonXi_11;
  TH1F              *hInvMassAntiProtonXi_12;
  TH1F              *hInvMassAntiProtonXi_13;
  TH1F              *hInvMassAntiProtonXi_14;
  TH1F              *hInvMassAntiProtonXi_15;
  TH1F              *hInvMassAntiProtonXi_16;
  TH1F              *hInvMassAntiProtonXi_17;
  //FG LambdaXi
  TH1F              *hInvMassLambdaXi_1;
  TH1F              *hInvMassLambdaXi_2;
  TH1F              *hInvMassLambdaXi_3;
  TH1F              *hInvMassLambdaXi_4;
  TH1F              *hInvMassLambdaXi_5;
  TH1F              *hInvMassLambdaXi_6;
  TH1F              *hInvMassLambdaXi_7;
  TH1F              *hInvMassLambdaXi_8;
  TH1F              *hInvMassLambdaXi_9;
  TH1F              *hInvMassLambdaXi_10;
  TH1F              *hInvMassLambdaXi_11;
  TH1F              *hInvMassLambdaXi_12;
  TH1F              *hInvMassLambdaXi_13;
  TH1F              *hInvMassLambdaXi_14;
  TH1F              *hInvMassLambdaXi_15;
  TH1F              *hInvMassAntiLambdaXi_1;
  TH1F              *hInvMassAntiLambdaXi_2;
  TH1F              *hInvMassAntiLambdaXi_3;
  TH1F              *hInvMassAntiLambdaXi_4;
  TH1F              *hInvMassAntiLambdaXi_5;
  TH1F              *hInvMassAntiLambdaXi_6;
  TH1F              *hInvMassAntiLambdaXi_7;
  TH1F              *hInvMassAntiLambdaXi_8;
  TH1F              *hInvMassAntiLambdaXi_9;
  TH1F              *hInvMassAntiLambdaXi_10;
  TH1F              *hInvMassAntiLambdaXi_11;
  TH1F              *hInvMassAntiLambdaXi_12;
  TH1F              *hInvMassAntiLambdaXi_13;
  TH1F              *hInvMassAntiLambdaXi_14;
  TH1F              *hInvMassAntiLambdaXi_15;
  //BG LambdaLambda                                                                                                                  
  TH1F              *hInvMassLambdaLambda_evtpool1;
  TH1F              *hInvMassLambdaLambda_evtpool2;
  TH1F              *hInvMassLambdaLambda_evtpool3;
  TH1F              *hInvMassLambdaLambda_evtpool4;
  TH1F              *hInvMassLambdaLambda_evtpool5;
  TH1F              *hInvMassLambdaLambda_evtpool6;
  TH1F              *hInvMassLambdaLambda_evtpool7;
  TH1F              *hInvMassLambdaLambda_evtpool8;
  TH1F              *hInvMassLambdaLambda_evtpool9;
  TH1F              *hInvMassLambdaLambda_evtpool10;
  TH1F              *hInvMassLambdaLambda_evtpool11;
  TH1F              *hInvMassLambdaLambda_evtpool12;
  TH1F              *hInvMassLambdaLambda_evtpool13;
  TH1F              *hInvMassLambdaLambda_evtpool14;
  TH1F              *hInvMassLambdaLambda_evtpool15;
  TH1F              *hInvMassLambdaLambda_evtpool16;
  TH1F              *hInvMassLambdaLambda_evtpool17;
  TH1F              *hInvMassAntiLambdaAntiLambda_evtpool1;
  TH1F              *hInvMassAntiLambdaAntiLambda_evtpool2;
  TH1F              *hInvMassAntiLambdaAntiLambda_evtpool3;
  TH1F              *hInvMassAntiLambdaAntiLambda_evtpool4;
  TH1F              *hInvMassAntiLambdaAntiLambda_evtpool5;
  TH1F              *hInvMassAntiLambdaAntiLambda_evtpool6;
  TH1F              *hInvMassAntiLambdaAntiLambda_evtpool7;
  TH1F              *hInvMassAntiLambdaAntiLambda_evtpool8;
  TH1F              *hInvMassAntiLambdaAntiLambda_evtpool9;
  TH1F              *hInvMassAntiLambdaAntiLambda_evtpool10;
  TH1F              *hInvMassAntiLambdaAntiLambda_evtpool11;
  TH1F              *hInvMassAntiLambdaAntiLambda_evtpool12;
  TH1F              *hInvMassAntiLambdaAntiLambda_evtpool13;
  TH1F              *hInvMassAntiLambdaAntiLambda_evtpool14;
  TH1F              *hInvMassAntiLambdaAntiLambda_evtpool15;
  TH1F              *hInvMassAntiLambdaAntiLambda_evtpool16;
  TH1F              *hInvMassAntiLambdaAntiLambda_evtpool17;
  //BG ProtonXi                                                                                                                    
  TH1F              *hInvMassProtonXi_evtpool1;
  TH1F              *hInvMassProtonXi_evtpool2;
  TH1F              *hInvMassProtonXi_evtpool3;
  TH1F              *hInvMassProtonXi_evtpool4;
  TH1F              *hInvMassProtonXi_evtpool5;
  TH1F              *hInvMassProtonXi_evtpool6;
  TH1F              *hInvMassProtonXi_evtpool7;
  TH1F              *hInvMassProtonXi_evtpool8;
  TH1F              *hInvMassProtonXi_evtpool9;
  TH1F              *hInvMassProtonXi_evtpool10;
  TH1F              *hInvMassProtonXi_evtpool11;
  TH1F              *hInvMassProtonXi_evtpool12;
  TH1F              *hInvMassProtonXi_evtpool13;
  TH1F              *hInvMassProtonXi_evtpool14;
  TH1F              *hInvMassProtonXi_evtpool15;
  TH1F              *hInvMassProtonXi_evtpool16;
  TH1F              *hInvMassProtonXi_evtpool17;
  TH1F              *hInvMassAntiProtonXi_evtpool1;
  TH1F              *hInvMassAntiProtonXi_evtpool2;
  TH1F              *hInvMassAntiProtonXi_evtpool3;
  TH1F              *hInvMassAntiProtonXi_evtpool4;
  TH1F              *hInvMassAntiProtonXi_evtpool5;
  TH1F              *hInvMassAntiProtonXi_evtpool6;
  TH1F              *hInvMassAntiProtonXi_evtpool7;
  TH1F              *hInvMassAntiProtonXi_evtpool8;
  TH1F              *hInvMassAntiProtonXi_evtpool9;
  TH1F              *hInvMassAntiProtonXi_evtpool10;
  TH1F              *hInvMassAntiProtonXi_evtpool11;
  TH1F              *hInvMassAntiProtonXi_evtpool12;
  TH1F              *hInvMassAntiProtonXi_evtpool13;
  TH1F              *hInvMassAntiProtonXi_evtpool14;
  TH1F              *hInvMassAntiProtonXi_evtpool15;
  TH1F              *hInvMassAntiProtonXi_evtpool16;
  TH1F              *hInvMassAntiProtonXi_evtpool17;
  //BG LambdaXi                                                                                                                     
  TH1F              *hInvMassLambdaXi_evtpool1;
  TH1F              *hInvMassLambdaXi_evtpool2;
  TH1F              *hInvMassLambdaXi_evtpool3;
  TH1F              *hInvMassLambdaXi_evtpool4;
  TH1F              *hInvMassLambdaXi_evtpool5;
  TH1F              *hInvMassLambdaXi_evtpool6;
  TH1F              *hInvMassLambdaXi_evtpool7;
  TH1F              *hInvMassLambdaXi_evtpool8;
  TH1F              *hInvMassLambdaXi_evtpool9;
  TH1F              *hInvMassLambdaXi_evtpool10;
  TH1F              *hInvMassLambdaXi_evtpool11;
  TH1F              *hInvMassLambdaXi_evtpool12;
  TH1F              *hInvMassLambdaXi_evtpool13;
  TH1F              *hInvMassLambdaXi_evtpool14;
  TH1F              *hInvMassLambdaXi_evtpool15;
  TH1F              *hInvMassAntiLambdaXi_evtpool1;
  TH1F              *hInvMassAntiLambdaXi_evtpool2;
  TH1F              *hInvMassAntiLambdaXi_evtpool3;
  TH1F              *hInvMassAntiLambdaXi_evtpool4;
  TH1F              *hInvMassAntiLambdaXi_evtpool5;
  TH1F              *hInvMassAntiLambdaXi_evtpool6;
  TH1F              *hInvMassAntiLambdaXi_evtpool7;
  TH1F              *hInvMassAntiLambdaXi_evtpool8;
  TH1F              *hInvMassAntiLambdaXi_evtpool9;
  TH1F              *hInvMassAntiLambdaXi_evtpool10;
  TH1F              *hInvMassAntiLambdaXi_evtpool11;
  TH1F              *hInvMassAntiLambdaXi_evtpool12;
  TH1F              *hInvMassAntiLambdaXi_evtpool13;
  TH1F              *hInvMassAntiLambdaXi_evtpool14;
  TH1F              *hInvMassAntiLambdaXi_evtpool15;


  //======
  AliAnalysisTaskPPiLambda(const AliAnalysisTaskPPiLambda&); // not implemented
  AliAnalysisTaskPPiLambda& operator=(const AliAnalysisTaskPPiLambda&); // not implemented

  ClassDef(AliAnalysisTaskPPiLambda,1);

};

#endif
