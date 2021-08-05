//Github test
#include "AliAnalysisTaskPPiLambda.h"

// std
#include <iostream>
#include <stdio.h>

// ROOT includes
#include <TList.h>
#include <TChain.h>
#include <TMath.h>
#include <TH1.h>
#include <TFile.h>

// AliRoot includes
#include <AliAnalysisTask.h>
#include <AliAnalysisManager.h>
#include <AliAnalysisUtils.h>
#include <AliAODEvent.h>
#include <AliESDEvent.h>
#include <AliVEvent.h>
#include <AliAODVZERO.h>
#include <AliAODVertex.h>
#include <AliAODInputHandler.h> 
#include <AliAODHandler.h> 
#include <AliAODHeader.h>
#include <AliAODTrack.h> 
#include <AliAODPid.h> 
#include <AliPIDResponse.h>
#include <AliAODv0.h>
#include <AliAODcascade.h>
#include <AliEventPoolManager.h>
#include <AliAODTracklets.h>

//#include <AliCentrality>

class AliAnalysisTaskPPiLambda;

using namespace std;

ClassImp(AliAnalysisTaskPPiLambda)

//_____________________________________________________________________________
AliAnalysisTaskPPiLambda::AliAnalysisTaskPPiLambda() : 
AliAnalysisTaskSE(),
  fAOD(0),
  fHeader(0),
  fUtils(0),
  fPIDResponse(0),
  ftrigBit(0),
  nevt(0),
  lambdacounter(0),
//
  fPoolManagerlambda(0),
  fPoolManagerantilambda(0),
  fPoolManagerxi(0),
  fPoolManagerxip(0),
  fPoolManagerproton(0),
  fPoolManagerantiproton(0),
  fPoolManagerlam(0),
  fPoolManagerantilam(0),
  fPoolManagerx(0),
  fPoolManagerxp(0),
  fOutputList(0),
  fEvtCounter(0),
  fEvtPassCut(0),
  fEvtVtxX(0),
  fEvtVtxY(0),
  fEvtVtxZ(0),
  fEvtVtxTrk(0),
  fSPDVtxZ(0),
  fSPDVtxTrk(0),
  fSPDVtxCor(0),
  fSPDVtxDisp(0),
//multiplicity                                                                                                         
  fmultiplicity(0),
  fmultiplicity_lambda(0),
  fmultiplicity_antilambda(0),
  fmultiplicity_xim(0),
  fmultiplicity_xip(0),
  fmultiplicity_proton(0),
  fmultiplicity_antiproton(0),
  fmultiplicity_lambdalambda(0),
  fmultiplicity_antilambdaantilambda(0),
  fmultiplicity_protonxi(0),
  fmultiplicity_antiprotonxi(0),
  fmultiplicity_lambdaxi(0),
  fmultiplicity_antilambdaxi(0),
//proton
  fTrkEta(0),
  fTrkPt(0),
  fTPCclus(0),
  fTPCcrsR(0),
  fTPCclusF(0),
  fTrkDCAxy(0),
  fTrkDCAz(0),
  fPvsTPCsignal(0),
  fPvsTOFsignal(0),
//proton select 
  fPt_proton(0),
  fPt_antiproton(0),
  PID_cut(0),
//v0 selection                                                                                    
  fPtv0(0),
  fPtv0_lambda(0),
  fTransRadius(0),
  fDCAdaugTov0Vtx(0),
  fCPA(0),
  hInvMassLambda(0),
  fPtv0_antilambda(0),
  fPtv0_antilambda_mass_masscut(0),
  hInvMassAntiLambda(0),
//invariant mass lambda    
  hInvMassXidecayLambda(0),
  hInvMassXidecayantiLambda(0),
//invariant mass xi                                                                                
  hInvMassXi(0),
  hInvMassposXi(0),
  XiCPA(0),
  XiTranverseradius(0),
  dcaXiDaughters(0),
  fPtcascade(0),
  fPtv0casc(0),
  fPtBachelor(0),
  fPhiXidecayLambda(0),
  fPtcascade_xip(0),
  fPtv0casc_xip(0),
  fPtBachelor_xip(0),
  fPhiXidecayLambda_xip(0),
//invariant mass dibaryon
  hInvMassProtonXi(0),
  hInvMassantiProtonXi(0),
  hInvMassLambdaXi(0),
  hInvMassAntiLambdaXip(0),
//decay rejection
  fPt_allLambda(0),                                                                            
  fPt_xidecaylambda(0),
  fPt_allAntiLambda(0),
  fPt_xidecayAntiLambda(0),
//decay rejection pair                                                                            
  hInvMassLambdaLambda_all(0),
  hInvMassLambdaLambda_onlyprompt(0),
  hInvMassAntiLambdaAntiLambda_all(0),
  hInvMassAntiLambdaAntiLambda_onlyprompt(0),
//eventmixing                                                                                                      
//hInvMassLambdaLambda_combinatorial(0),
  hInvMassLambdaLambda_evtpool(0),
  hInvMassAntiLambdaAntiLambda_evtpool(0),
  hInvMassProtonXi_evtpool(0),
  hInvMassAntiProtonXi_evtpool(0),
  hInvMassLambdaXi_evtpool(0),
  hInvMassAntiLambdaXi_evtpool(0),
//vs relative momentum                                                                                   
  hInvMassLambdaLambda_relmomentum(0),
  hInvMassLambdaLambda_openangle(0),
  hInvMassAntiLambdaAntiLambda_relmomentum(0),
  hInvMassAntiLambdaAntiLambda_openangle(0),
  hAngle_relmom(0),
  hInvMassProtonXi_relmomentum(0),
  hInvMassProtonXi_openangle(0),
  hInvMassantiProtonXi_relmomentum(0),
  hInvMassantiProtonXi_openangle(0),
  hInvMassLambdaXi_relmomentum(0),
  hInvMassLambdaXi_openangle(0),
  hInvMassAntiLambdaXi_relmomentum(0),
  hInvMassAntiLambdaXi_openangle(0),
//pt correlation
  hPtcorrelation_lamlam(0),
  hPtcorrelation_pxi(0),
  hPtcorrelation_lamxi(0),
//lowmass                                                                                             
  hInvMassProtonXi_lowmassx(0),
  hInvMassProtonXi_lowmassy(0),
  hInvMassProtonXi_lowmassz(0),
  hInvMassantiProtonXi_lowmassx(0),
  hInvMassantiProtonXi_lowmassy(0),
  hInvMassantiProtonXi_lowmassz(0),
  //FG LambdaLambda                                                                                                                     
  hInvMassLambdaLambda_1(0),
  hInvMassLambdaLambda_2(0),
  hInvMassLambdaLambda_3(0),
  hInvMassLambdaLambda_4(0),
  hInvMassLambdaLambda_5(0),
  hInvMassLambdaLambda_6(0),
  hInvMassLambdaLambda_7(0),
  hInvMassLambdaLambda_8(0),
  hInvMassLambdaLambda_9(0),
  hInvMassLambdaLambda_10(0),
  hInvMassLambdaLambda_11(0),
  hInvMassLambdaLambda_12(0),
  hInvMassLambdaLambda_13(0),
  hInvMassLambdaLambda_14(0),
  hInvMassLambdaLambda_15(0),
  hInvMassLambdaLambda_16(0),
  hInvMassLambdaLambda_17(0),
  hInvMassAntiLambdaAntiLambda_1(0),
  hInvMassAntiLambdaAntiLambda_2(0),
  hInvMassAntiLambdaAntiLambda_3(0),
  hInvMassAntiLambdaAntiLambda_4(0),
  hInvMassAntiLambdaAntiLambda_5(0),
  hInvMassAntiLambdaAntiLambda_6(0),
  hInvMassAntiLambdaAntiLambda_7(0),
  hInvMassAntiLambdaAntiLambda_8(0),
  hInvMassAntiLambdaAntiLambda_9(0),
  hInvMassAntiLambdaAntiLambda_10(0),
  hInvMassAntiLambdaAntiLambda_11(0),
  hInvMassAntiLambdaAntiLambda_12(0),
  hInvMassAntiLambdaAntiLambda_13(0),
  hInvMassAntiLambdaAntiLambda_14(0),
  hInvMassAntiLambdaAntiLambda_15(0),
  hInvMassAntiLambdaAntiLambda_16(0),
  hInvMassAntiLambdaAntiLambda_17(0),
  //FG ProtonXi                                                                                                                     
  hInvMassProtonXi_1(0),
  hInvMassProtonXi_2(0),
  hInvMassProtonXi_3(0),
  hInvMassProtonXi_4(0),
  hInvMassProtonXi_5(0),
  hInvMassProtonXi_6(0),
  hInvMassProtonXi_7(0),
  hInvMassProtonXi_8(0),
  hInvMassProtonXi_9(0),
  hInvMassProtonXi_10(0),
  hInvMassProtonXi_11(0),
  hInvMassProtonXi_12(0),
  hInvMassProtonXi_13(0),
  hInvMassProtonXi_14(0),
  hInvMassProtonXi_15(0),
  hInvMassProtonXi_16(0),
  hInvMassProtonXi_17(0),
  hInvMassAntiProtonXi_1(0),
  hInvMassAntiProtonXi_2(0),
  hInvMassAntiProtonXi_3(0),
  hInvMassAntiProtonXi_4(0),
  hInvMassAntiProtonXi_5(0),
  hInvMassAntiProtonXi_6(0),
  hInvMassAntiProtonXi_7(0),
  hInvMassAntiProtonXi_8(0),
  hInvMassAntiProtonXi_9(0),
  hInvMassAntiProtonXi_10(0),
  hInvMassAntiProtonXi_11(0),
  hInvMassAntiProtonXi_12(0),
  hInvMassAntiProtonXi_13(0),
  hInvMassAntiProtonXi_14(0),
  hInvMassAntiProtonXi_15(0),
  hInvMassAntiProtonXi_16(0),
  hInvMassAntiProtonXi_17(0),
  //FG LambdaXi                                                                                                                     
  hInvMassLambdaXi_1(0),
  hInvMassLambdaXi_2(0),
  hInvMassLambdaXi_3(0),
  hInvMassLambdaXi_4(0),
  hInvMassLambdaXi_5(0),
  hInvMassLambdaXi_6(0),
  hInvMassLambdaXi_7(0),
  hInvMassLambdaXi_8(0),
  hInvMassLambdaXi_9(0),
  hInvMassLambdaXi_10(0),
  hInvMassLambdaXi_11(0),
  hInvMassLambdaXi_12(0),
  hInvMassLambdaXi_13(0),
  hInvMassLambdaXi_14(0),
  hInvMassLambdaXi_15(0),
  hInvMassAntiLambdaXi_1(0),
  hInvMassAntiLambdaXi_2(0),
  hInvMassAntiLambdaXi_3(0),
  hInvMassAntiLambdaXi_4(0),
  hInvMassAntiLambdaXi_5(0),
  hInvMassAntiLambdaXi_6(0),
  hInvMassAntiLambdaXi_7(0),
  hInvMassAntiLambdaXi_8(0),
  hInvMassAntiLambdaXi_9(0),
  hInvMassAntiLambdaXi_10(0),
  hInvMassAntiLambdaXi_11(0),
  hInvMassAntiLambdaXi_12(0),
  hInvMassAntiLambdaXi_13(0),
  hInvMassAntiLambdaXi_14(0),
  hInvMassAntiLambdaXi_15(0),
  //BG LambdaLambda
  hInvMassLambdaLambda_evtpool1(0),
  hInvMassLambdaLambda_evtpool2(0),
  hInvMassLambdaLambda_evtpool3(0),
  hInvMassLambdaLambda_evtpool4(0),
  hInvMassLambdaLambda_evtpool5(0),
  hInvMassLambdaLambda_evtpool6(0),
  hInvMassLambdaLambda_evtpool7(0),
  hInvMassLambdaLambda_evtpool8(0),
  hInvMassLambdaLambda_evtpool9(0),
  hInvMassLambdaLambda_evtpool10(0),
  hInvMassLambdaLambda_evtpool11(0),
  hInvMassLambdaLambda_evtpool12(0),
  hInvMassLambdaLambda_evtpool13(0),
  hInvMassLambdaLambda_evtpool14(0),
  hInvMassLambdaLambda_evtpool15(0),
  hInvMassLambdaLambda_evtpool16(0),
  hInvMassLambdaLambda_evtpool17(0),
  hInvMassAntiLambdaAntiLambda_evtpool1(0),
  hInvMassAntiLambdaAntiLambda_evtpool2(0),
  hInvMassAntiLambdaAntiLambda_evtpool3(0),
  hInvMassAntiLambdaAntiLambda_evtpool4(0),
  hInvMassAntiLambdaAntiLambda_evtpool5(0),
  hInvMassAntiLambdaAntiLambda_evtpool6(0),
  hInvMassAntiLambdaAntiLambda_evtpool7(0),
  hInvMassAntiLambdaAntiLambda_evtpool8(0),
  hInvMassAntiLambdaAntiLambda_evtpool9(0),
  hInvMassAntiLambdaAntiLambda_evtpool10(0),
  hInvMassAntiLambdaAntiLambda_evtpool11(0),
  hInvMassAntiLambdaAntiLambda_evtpool12(0),
  hInvMassAntiLambdaAntiLambda_evtpool13(0),
  hInvMassAntiLambdaAntiLambda_evtpool14(0),
  hInvMassAntiLambdaAntiLambda_evtpool15(0),
  hInvMassAntiLambdaAntiLambda_evtpool16(0),
  hInvMassAntiLambdaAntiLambda_evtpool17(0),
//BG ProtonXi
  hInvMassProtonXi_evtpool1(0),
  hInvMassProtonXi_evtpool2(0),
  hInvMassProtonXi_evtpool3(0),
  hInvMassProtonXi_evtpool4(0),
  hInvMassProtonXi_evtpool5(0),
  hInvMassProtonXi_evtpool6(0),
  hInvMassProtonXi_evtpool7(0),
  hInvMassProtonXi_evtpool8(0),
  hInvMassProtonXi_evtpool9(0),
  hInvMassProtonXi_evtpool10(0),
  hInvMassProtonXi_evtpool11(0),
  hInvMassProtonXi_evtpool12(0),
  hInvMassProtonXi_evtpool13(0),
  hInvMassProtonXi_evtpool14(0),
  hInvMassProtonXi_evtpool15(0),
  hInvMassProtonXi_evtpool16(0),
  hInvMassProtonXi_evtpool17(0),
  hInvMassAntiProtonXi_evtpool1(0),
  hInvMassAntiProtonXi_evtpool2(0),
  hInvMassAntiProtonXi_evtpool3(0),
  hInvMassAntiProtonXi_evtpool4(0),
  hInvMassAntiProtonXi_evtpool5(0),
  hInvMassAntiProtonXi_evtpool6(0),
  hInvMassAntiProtonXi_evtpool7(0),
  hInvMassAntiProtonXi_evtpool8(0),
  hInvMassAntiProtonXi_evtpool9(0),
  hInvMassAntiProtonXi_evtpool10(0),
  hInvMassAntiProtonXi_evtpool11(0),
  hInvMassAntiProtonXi_evtpool12(0),
  hInvMassAntiProtonXi_evtpool13(0),
  hInvMassAntiProtonXi_evtpool14(0),
  hInvMassAntiProtonXi_evtpool15(0),
  hInvMassAntiProtonXi_evtpool16(0),
  hInvMassAntiProtonXi_evtpool17(0),
//BG LambdaXi
  hInvMassLambdaXi_evtpool1(0),
  hInvMassLambdaXi_evtpool2(0),
  hInvMassLambdaXi_evtpool3(0),
  hInvMassLambdaXi_evtpool4(0),
  hInvMassLambdaXi_evtpool5(0),
  hInvMassLambdaXi_evtpool6(0),
  hInvMassLambdaXi_evtpool7(0),
  hInvMassLambdaXi_evtpool8(0),
  hInvMassLambdaXi_evtpool9(0),
  hInvMassLambdaXi_evtpool10(0),
  hInvMassLambdaXi_evtpool11(0),
  hInvMassLambdaXi_evtpool12(0),
  hInvMassLambdaXi_evtpool13(0),
  hInvMassLambdaXi_evtpool14(0),
  hInvMassLambdaXi_evtpool15(0),
  hInvMassAntiLambdaXi_evtpool1(0),
  hInvMassAntiLambdaXi_evtpool2(0),
  hInvMassAntiLambdaXi_evtpool3(0),
  hInvMassAntiLambdaXi_evtpool4(0),
  hInvMassAntiLambdaXi_evtpool5(0),
  hInvMassAntiLambdaXi_evtpool6(0),
  hInvMassAntiLambdaXi_evtpool7(0),
  hInvMassAntiLambdaXi_evtpool8(0),
  hInvMassAntiLambdaXi_evtpool9(0),
  hInvMassAntiLambdaXi_evtpool10(0),
  hInvMassAntiLambdaXi_evtpool11(0),
  hInvMassAntiLambdaXi_evtpool12(0),
  hInvMassAntiLambdaXi_evtpool13(0),
  hInvMassAntiLambdaXi_evtpool14(0),
  hInvMassAntiLambdaXi_evtpool15(0)

{
 
}
 
//_____________________________________________________________________________
AliAnalysisTaskPPiLambda::AliAnalysisTaskPPiLambda(const char *name) : 
  AliAnalysisTaskSE(name),
  fAOD(0),
  fHeader(0),
  fUtils(0),
  fPIDResponse(0),
  ftrigBit(0),
  nevt(0),
  lambdacounter(0),
  //
  fPoolManagerlambda(0),
  fPoolManagerantilambda(0),
  fPoolManagerxi(0),
  fPoolManagerxip(0),
  fPoolManagerproton(0),
  fPoolManagerantiproton(0),
  fPoolManagerlam(0),
  fPoolManagerantilam(0),
  fPoolManagerx(0),
  fPoolManagerxp(0),
  fOutputList(0),
  fEvtCounter(0),
  fEvtPassCut(0),
  fEvtVtxX(0),
  fEvtVtxY(0),
  fEvtVtxZ(0),
  fEvtVtxTrk(0),
  fSPDVtxZ(0),
  fSPDVtxTrk(0),
  fSPDVtxCor(0),
  fSPDVtxDisp(0),
  //multiplicity                                                                                                            
  fmultiplicity(0),
  fmultiplicity_lambda(0),
  fmultiplicity_antilambda(0),
  fmultiplicity_xim(0),
  fmultiplicity_xip(0),
  fmultiplicity_proton(0),
  fmultiplicity_antiproton(0),
  fmultiplicity_lambdalambda(0),
  fmultiplicity_antilambdaantilambda(0),
  fmultiplicity_protonxi(0),
  fmultiplicity_antiprotonxi(0),
  fmultiplicity_lambdaxi(0),
  fmultiplicity_antilambdaxi(0),
  //proton
  fTrkEta(0),
  fTrkPt(0),
  fTPCclus(0),
  fTPCcrsR(0),
  fTPCclusF(0),
  fTrkDCAxy(0),
  fTrkDCAz(0),
  fPvsTPCsignal(0),
  fPvsTOFsignal(0),
  //proton select 
  fPt_proton(0),
  fPt_antiproton(0),
  PID_cut(0),
  //v0 selection                                                                                   
  fPtv0(0),
  fPtv0_lambda(0),
  fTransRadius(0),
  fDCAdaugTov0Vtx(0),
  fCPA(0),
  hInvMassLambda(0),
  fPtv0_antilambda(0),
  fPtv0_antilambda_mass_masscut(0),
  hInvMassAntiLambda(0),
  //invariant mass lambda                                                                          
  hInvMassXidecayLambda(0),
  hInvMassXidecayantiLambda(0),
  //invariant mass xi                                                                              
  hInvMassXi(0),
  hInvMassposXi(0),
  XiCPA(0),
  XiTranverseradius(0),
  dcaXiDaughters(0),
  fPtcascade(0),
  fPtv0casc(0),
  fPtBachelor(0),
  fPhiXidecayLambda(0),
  fPtcascade_xip(0),
  fPtv0casc_xip(0),
  fPtBachelor_xip(0),
  fPhiXidecayLambda_xip(0),
  //invariant mass dibaryon
  hInvMassProtonXi(0),
  hInvMassantiProtonXi(0),
  hInvMassLambdaXi(0),
  hInvMassAntiLambdaXip(0),
  //decay rejection
  fPt_allLambda(0),
  fPt_xidecaylambda(0),
  fPt_allAntiLambda(0),
  fPt_xidecayAntiLambda(0),
  //decay rejection pair                                                                        
  hInvMassLambdaLambda_all(0),
  hInvMassLambdaLambda_onlyprompt(0),
  hInvMassAntiLambdaAntiLambda_all(0),
  hInvMassAntiLambdaAntiLambda_onlyprompt(0),
  //eventmixing                                                                                                      
  //hInvMassLambdaLambda_combinatorial(0),
  hInvMassLambdaLambda_evtpool(0),
  hInvMassAntiLambdaAntiLambda_evtpool(0),
  hInvMassProtonXi_evtpool(0),
  hInvMassAntiProtonXi_evtpool(0),
  hInvMassLambdaXi_evtpool(0),
  hInvMassAntiLambdaXi_evtpool(0),
  //vs relative momentum                                                                                    
  hInvMassLambdaLambda_relmomentum(0),
  hInvMassLambdaLambda_openangle(0),
  hInvMassAntiLambdaAntiLambda_relmomentum(0),
  hInvMassAntiLambdaAntiLambda_openangle(0),
  hAngle_relmom(0),
  hInvMassProtonXi_relmomentum(0),
  hInvMassProtonXi_openangle(0),
  hInvMassantiProtonXi_relmomentum(0),
  hInvMassantiProtonXi_openangle(0),
  hInvMassLambdaXi_relmomentum(0),
  hInvMassLambdaXi_openangle(0),
  hInvMassAntiLambdaXi_relmomentum(0),
  hInvMassAntiLambdaXi_openangle(0),
//pt correlation
  hPtcorrelation_lamlam(0),
  hPtcorrelation_pxi(0),
  hPtcorrelation_lamxi(0),
  //lowmass                                                                                             
  hInvMassProtonXi_lowmassx(0),
  hInvMassProtonXi_lowmassy(0),
  hInvMassProtonXi_lowmassz(0),
  hInvMassantiProtonXi_lowmassx(0),
  hInvMassantiProtonXi_lowmassy(0),
  hInvMassantiProtonXi_lowmassz(0),
  //FG LambdaLambda                                                                                                                     
  hInvMassLambdaLambda_1(0),
  hInvMassLambdaLambda_2(0),
  hInvMassLambdaLambda_3(0),
  hInvMassLambdaLambda_4(0),
  hInvMassLambdaLambda_5(0),
  hInvMassLambdaLambda_6(0),
  hInvMassLambdaLambda_7(0),
  hInvMassLambdaLambda_8(0),
  hInvMassLambdaLambda_9(0),
  hInvMassLambdaLambda_10(0),
  hInvMassLambdaLambda_11(0),
  hInvMassLambdaLambda_12(0),
  hInvMassLambdaLambda_13(0),
  hInvMassLambdaLambda_14(0),
  hInvMassLambdaLambda_15(0),
  hInvMassLambdaLambda_16(0),
  hInvMassLambdaLambda_17(0),
  hInvMassAntiLambdaAntiLambda_1(0),
  hInvMassAntiLambdaAntiLambda_2(0),
  hInvMassAntiLambdaAntiLambda_3(0),
  hInvMassAntiLambdaAntiLambda_4(0),
  hInvMassAntiLambdaAntiLambda_5(0),
  hInvMassAntiLambdaAntiLambda_6(0),
  hInvMassAntiLambdaAntiLambda_7(0),
  hInvMassAntiLambdaAntiLambda_8(0),
  hInvMassAntiLambdaAntiLambda_9(0),
  hInvMassAntiLambdaAntiLambda_10(0),
  hInvMassAntiLambdaAntiLambda_11(0),
  hInvMassAntiLambdaAntiLambda_12(0),
  hInvMassAntiLambdaAntiLambda_13(0),
  hInvMassAntiLambdaAntiLambda_14(0),
  hInvMassAntiLambdaAntiLambda_15(0),
  hInvMassAntiLambdaAntiLambda_16(0),
  hInvMassAntiLambdaAntiLambda_17(0),
  //FG ProtonXi                                                                                                                     
  hInvMassProtonXi_1(0),
  hInvMassProtonXi_2(0),
  hInvMassProtonXi_3(0),
  hInvMassProtonXi_4(0),
  hInvMassProtonXi_5(0),
  hInvMassProtonXi_6(0),
  hInvMassProtonXi_7(0),
  hInvMassProtonXi_8(0),
  hInvMassProtonXi_9(0),
  hInvMassProtonXi_10(0),
  hInvMassProtonXi_11(0),
  hInvMassProtonXi_12(0),
  hInvMassProtonXi_13(0),
  hInvMassProtonXi_14(0),
  hInvMassProtonXi_15(0),
  hInvMassProtonXi_16(0),
  hInvMassProtonXi_17(0),
  hInvMassAntiProtonXi_1(0),
  hInvMassAntiProtonXi_2(0),
  hInvMassAntiProtonXi_3(0),
  hInvMassAntiProtonXi_4(0),
  hInvMassAntiProtonXi_5(0),
  hInvMassAntiProtonXi_6(0),
  hInvMassAntiProtonXi_7(0),
  hInvMassAntiProtonXi_8(0),
  hInvMassAntiProtonXi_9(0),
  hInvMassAntiProtonXi_10(0),
  hInvMassAntiProtonXi_11(0),
  hInvMassAntiProtonXi_12(0),
  hInvMassAntiProtonXi_13(0),
  hInvMassAntiProtonXi_14(0),
  hInvMassAntiProtonXi_15(0),
  hInvMassAntiProtonXi_16(0),
  hInvMassAntiProtonXi_17(0),
  //FG LambdaXi                                                                                                                     
  hInvMassLambdaXi_1(0),
  hInvMassLambdaXi_2(0),
  hInvMassLambdaXi_3(0),
  hInvMassLambdaXi_4(0),
  hInvMassLambdaXi_5(0),
  hInvMassLambdaXi_6(0),
  hInvMassLambdaXi_7(0),
  hInvMassLambdaXi_8(0),
  hInvMassLambdaXi_9(0),
  hInvMassLambdaXi_10(0),
  hInvMassLambdaXi_11(0),
  hInvMassLambdaXi_12(0),
  hInvMassLambdaXi_13(0),
  hInvMassLambdaXi_14(0),
  hInvMassLambdaXi_15(0),
  hInvMassAntiLambdaXi_1(0),
  hInvMassAntiLambdaXi_2(0),
  hInvMassAntiLambdaXi_3(0),
  hInvMassAntiLambdaXi_4(0),
  hInvMassAntiLambdaXi_5(0),
  hInvMassAntiLambdaXi_6(0),
  hInvMassAntiLambdaXi_7(0),
  hInvMassAntiLambdaXi_8(0),
  hInvMassAntiLambdaXi_9(0),
  hInvMassAntiLambdaXi_10(0),
  hInvMassAntiLambdaXi_11(0),
  hInvMassAntiLambdaXi_12(0),
  hInvMassAntiLambdaXi_13(0),
  hInvMassAntiLambdaXi_14(0),
  hInvMassAntiLambdaXi_15(0),
  //BG LambdaLambda
  hInvMassLambdaLambda_evtpool1(0),
  hInvMassLambdaLambda_evtpool2(0),
  hInvMassLambdaLambda_evtpool3(0),
  hInvMassLambdaLambda_evtpool4(0),
  hInvMassLambdaLambda_evtpool5(0),
  hInvMassLambdaLambda_evtpool6(0),
  hInvMassLambdaLambda_evtpool7(0),
  hInvMassLambdaLambda_evtpool8(0),
  hInvMassLambdaLambda_evtpool9(0),
  hInvMassLambdaLambda_evtpool10(0),
  hInvMassLambdaLambda_evtpool11(0),
  hInvMassLambdaLambda_evtpool12(0),
  hInvMassLambdaLambda_evtpool13(0),
  hInvMassLambdaLambda_evtpool14(0),
  hInvMassLambdaLambda_evtpool15(0),
  hInvMassLambdaLambda_evtpool16(0),
  hInvMassLambdaLambda_evtpool17(0),
  hInvMassAntiLambdaAntiLambda_evtpool1(0),
  hInvMassAntiLambdaAntiLambda_evtpool2(0),
  hInvMassAntiLambdaAntiLambda_evtpool3(0),
  hInvMassAntiLambdaAntiLambda_evtpool4(0),
  hInvMassAntiLambdaAntiLambda_evtpool5(0),
  hInvMassAntiLambdaAntiLambda_evtpool6(0),
  hInvMassAntiLambdaAntiLambda_evtpool7(0),
  hInvMassAntiLambdaAntiLambda_evtpool8(0),
  hInvMassAntiLambdaAntiLambda_evtpool9(0),
  hInvMassAntiLambdaAntiLambda_evtpool10(0),
  hInvMassAntiLambdaAntiLambda_evtpool11(0),
  hInvMassAntiLambdaAntiLambda_evtpool12(0),
  hInvMassAntiLambdaAntiLambda_evtpool13(0),
  hInvMassAntiLambdaAntiLambda_evtpool14(0),
  hInvMassAntiLambdaAntiLambda_evtpool15(0),
  hInvMassAntiLambdaAntiLambda_evtpool16(0),
  hInvMassAntiLambdaAntiLambda_evtpool17(0),
//BG ProtonXi
  hInvMassProtonXi_evtpool1(0),
  hInvMassProtonXi_evtpool2(0),
  hInvMassProtonXi_evtpool3(0),
  hInvMassProtonXi_evtpool4(0),
  hInvMassProtonXi_evtpool5(0),
  hInvMassProtonXi_evtpool6(0),
  hInvMassProtonXi_evtpool7(0),
  hInvMassProtonXi_evtpool8(0),
  hInvMassProtonXi_evtpool9(0),
  hInvMassProtonXi_evtpool10(0),
  hInvMassProtonXi_evtpool11(0),
  hInvMassProtonXi_evtpool12(0),
  hInvMassProtonXi_evtpool13(0),
  hInvMassProtonXi_evtpool14(0),
  hInvMassProtonXi_evtpool15(0),
  hInvMassProtonXi_evtpool16(0),
  hInvMassProtonXi_evtpool17(0),
  hInvMassAntiProtonXi_evtpool1(0),
  hInvMassAntiProtonXi_evtpool2(0),
  hInvMassAntiProtonXi_evtpool3(0),
  hInvMassAntiProtonXi_evtpool4(0),
  hInvMassAntiProtonXi_evtpool5(0),
  hInvMassAntiProtonXi_evtpool6(0),
  hInvMassAntiProtonXi_evtpool7(0),
  hInvMassAntiProtonXi_evtpool8(0),
  hInvMassAntiProtonXi_evtpool9(0),
  hInvMassAntiProtonXi_evtpool10(0),
  hInvMassAntiProtonXi_evtpool11(0),
  hInvMassAntiProtonXi_evtpool12(0),
  hInvMassAntiProtonXi_evtpool13(0),
  hInvMassAntiProtonXi_evtpool14(0),
  hInvMassAntiProtonXi_evtpool15(0),
  hInvMassAntiProtonXi_evtpool16(0),
  hInvMassAntiProtonXi_evtpool17(0),
//BG LambdaXi
  hInvMassLambdaXi_evtpool1(0),
  hInvMassLambdaXi_evtpool2(0),
  hInvMassLambdaXi_evtpool3(0),
  hInvMassLambdaXi_evtpool4(0),
  hInvMassLambdaXi_evtpool5(0),
  hInvMassLambdaXi_evtpool6(0),
  hInvMassLambdaXi_evtpool7(0),
  hInvMassLambdaXi_evtpool8(0),
  hInvMassLambdaXi_evtpool9(0),
  hInvMassLambdaXi_evtpool10(0),
  hInvMassLambdaXi_evtpool11(0),
  hInvMassLambdaXi_evtpool12(0),
  hInvMassLambdaXi_evtpool13(0),
  hInvMassLambdaXi_evtpool14(0),
  hInvMassLambdaXi_evtpool15(0),
  hInvMassAntiLambdaXi_evtpool1(0),
  hInvMassAntiLambdaXi_evtpool2(0),
  hInvMassAntiLambdaXi_evtpool3(0),
  hInvMassAntiLambdaXi_evtpool4(0),
  hInvMassAntiLambdaXi_evtpool5(0),
  hInvMassAntiLambdaXi_evtpool6(0),
  hInvMassAntiLambdaXi_evtpool7(0),
  hInvMassAntiLambdaXi_evtpool8(0),
  hInvMassAntiLambdaXi_evtpool9(0),
  hInvMassAntiLambdaXi_evtpool10(0),
  hInvMassAntiLambdaXi_evtpool11(0),
  hInvMassAntiLambdaXi_evtpool12(0),
  hInvMassAntiLambdaXi_evtpool13(0),
  hInvMassAntiLambdaXi_evtpool14(0),
  hInvMassAntiLambdaXi_evtpool15(0)

{
  // Constructor
  DefineInput(0,TChain::Class());
  DefineOutput(1,TList::Class());

  nevt=0;
  lambdacounter=0;

  for(Int_t i=0; i<5; i++){
    for(Int_t j=0; j<10; j++){
      eventdepth[i][j] =0;
      eventdepth2[i][j]=0;
      for(Int_t k=0; k<5; k++){
	mnLambda[i][j][k]=0;
	for(Int_t l=0; l<10; l++){

	  mEnergy[i][j][k][l]=0;
	  mPx[i][j][k][l]=0;
	  mPy[i][j][k][l]=0;
	  mPz[i][j][k][l]=0;	  
	}
      }
    }
  }

  
}

//________________________________________________________________________
AliAnalysisTaskPPiLambda::~AliAnalysisTaskPPiLambda()
{
  // destructor
  if(fOutputList){
    delete fOutputList;
  }

}

//________________________________________________________________________
void AliAnalysisTaskPPiLambda::UserCreateOutputObjects()
{

  //==========define AliEventPoolMgr
  Int_t MaxNEvents = 5;
  Int_t MaxNLambda = 10;
  //Int_t nMultiBins = 12;
  //Double_t multBins[]={0,10,20,30,40,50,60,70,80,90,100,110,120};
  Int_t nMultiBins = 17;
  Double_t multBins[]={1,4,8,12,16,20,24,28,32,36,40,50,60,70,80,90,100,200};
  Int_t nZvtxBins = 4;
  Double_t vertexBins[]={-10,-5,0,5,10};

  Int_t nMultiBinslamxi = 15;
  Double_t multBinslamxi[]={1,5,10,15,20,25,30,35,40,50,60,70,80,90,100,200};

  fPoolManagerlambda     =new AliEventPoolManager(MaxNEvents,MaxNLambda,nMultiBins,(Double_t*)multBins,nZvtxBins,(Double_t *)vertexBins);
  fPoolManagerantilambda =new AliEventPoolManager(MaxNEvents,MaxNLambda,nMultiBins,(Double_t*)multBins,nZvtxBins,(Double_t *)vertexBins);

  fPoolManagerxi         =new AliEventPoolManager(MaxNEvents,MaxNLambda,nMultiBins,(Double_t*)multBins,nZvtxBins,(Double_t *)vertexBins);
  fPoolManagerxip        =new AliEventPoolManager(MaxNEvents,MaxNLambda,nMultiBins,(Double_t*)multBins,nZvtxBins,(Double_t *)vertexBins);
  fPoolManagerproton     =new AliEventPoolManager(MaxNEvents,MaxNLambda,nMultiBins,(Double_t*)multBins,nZvtxBins,(Double_t *)vertexBins);
  fPoolManagerantiproton =new AliEventPoolManager(MaxNEvents,MaxNLambda,nMultiBins,(Double_t*)multBins,nZvtxBins,(Double_t *)vertexBins);

  fPoolManagerlam        =new AliEventPoolManager(MaxNEvents,MaxNLambda,nMultiBinslamxi,(Double_t*)multBinslamxi,nZvtxBins,(Double_t *)vertexBins);
  fPoolManagerantilam    =new AliEventPoolManager(MaxNEvents,MaxNLambda,nMultiBinslamxi,(Double_t*)multBinslamxi,nZvtxBins,(Double_t *)vertexBins);
  fPoolManagerx          =new AliEventPoolManager(MaxNEvents,MaxNLambda,nMultiBinslamxi,(Double_t*)multBinslamxi,nZvtxBins,(Double_t *)vertexBins);
  fPoolManagerxp         =new AliEventPoolManager(MaxNEvents,MaxNLambda,nMultiBinslamxi,(Double_t*)multBinslamxi,nZvtxBins,(Double_t *)vertexBins);
  //==========

  fOutputList=new TList();
  fOutputList->SetOwner(kTRUE);

  // Histograms
  fEvtCounter                  =new TH1F("fEvtCounter","",1,0,1);
  fEvtPassCut                  =new TH1F("fEvtPassCut","",1,0,1);
  fEvtVtxX                     =new TH1F("fEvtVtxX","",400,-20.,20.);
  fEvtVtxY                     =new TH1F("fEvtVtxY","",400,-20.,20.);
  fEvtVtxZ                     =new TH1F("fEvtVtxZ","",400,-20.,20.);
  fEvtVtxTrk                   =new TH1F("fEvtVtxTrk","",1000,0.,1000.);
  fSPDVtxZ                     =new TH1F("fSPDVtxZ","",400,-20.,20.);
  fSPDVtxTrk                   =new TH1F("fSPDVtxTrk","",1000,0.,1000.);
  fSPDVtxCor                   =new TH2F("fSPDVtxCor","",400,-20.,20.,400,-20.,20.);
  fSPDVtxDisp                  =new TH1F("fSPDVtxDisp","",300,0.,1.5);
  //multiplicity
  fmultiplicity                      =new TH2F("fmultiplicity","",2000,0,200,40,-20,20);
  fmultiplicity_lambda               =new TH1F("fmultiplicity_lambda","",2000,0,200);
  fmultiplicity_antilambda           =new TH1F("fmultiplicity_antilambda","",2000,0,200);
  fmultiplicity_xim                  =new TH1F("fmultiplicity_xim","",2000,0,200);
  fmultiplicity_xip                  =new TH1F("fmultiplicity_xip","",2000,0,200);
  fmultiplicity_proton               =new TH1F("fmultiplicity_proton","",2000,0,200);
  fmultiplicity_antiproton           =new TH1F("fmultiplicity_antiproton","",2000,0,200);
  fmultiplicity_lambdalambda         =new TH1F("fmultiplicity_lambdalambda","",2000,0,200);
  fmultiplicity_antilambdaantilambda =new TH1F("fmultiplicity_antilambdaantilambda","",2000,0,200);
  fmultiplicity_protonxi             =new TH1F("fmultiplicity_protonxi","",2000,0,200);
  fmultiplicity_antiprotonxi         =new TH1F("fmultiplicity_antiprotonxi","",2000,0,200);
  fmultiplicity_lambdaxi             =new TH1F("fmultiplicity_lambdaxi","",2000,0,200);
  fmultiplicity_antilambdaxi         =new TH1F("fmultiplicity_antilambdaxi","",2000,0,200);
  //proton
  fTrkEta                      =new TH1F("fTrkEta","",100,-2,2);
  fTrkPt                       =new TH1F("fTrkPt","",1000,0,5);
  fTPCclus                     =new TH1F("fTrkTPCclus","",1000,0,200);
  fTPCcrsR                     =new TH1F("fTrkTPCcrsR","",1000,0,200);
  fTPCclusF                    =new TH1F("fTrkTPCclusF","",1000,0,200);
  fTrkDCAxy                    =new TH1F("fTrkDCAxy","",100,-3,3);
  fTrkDCAz                     =new TH1F("fTrkDCAz","",100,-3,3);
  fPvsTPCsignal                =new TH2F("fPvsTPCsignal","",2000,0,10,20000,0,1000);
  fPvsTOFsignal                =new TH2F("fPvsTOFsignal","",500,0,5,200,0,2);
  //proton select 
  fPt_proton                   =new TH1F("fPt_proton","",1000,0,10); 
  fPt_antiproton               =new TH1F("fPt_antiproton","",1000,0,10); 
  PID_cut                      =new TH2F("PID_cut","",2000,0,10,20000,0,1000);
  //v0 selection
  fPtv0                        =new TH2F("fPtv0","",1000,0,10,1000,0,10);
  fPtv0_lambda                 =new TH2F("fPtv0_lambda","",1000,0,10,1000,0,10);
  fTransRadius                 =new TH1F("fTransRadius","",1000,0,300);
  fDCAdaugTov0Vtx              =new TH1F("fDCAdaugTov0Vtx","",100,0,10);
  fCPA                         =new TH1F("fCPA","",100,-1,1);
  hInvMassLambda               =new TH1F("hInvMassLambda","",10000,0,10);
  fPtv0_antilambda             =new TH1F("fPtv0_antilambda","",10000,0,10);
  fPtv0_antilambda_mass_masscut=new TH2F("fPtv0_antilambda_mass_masscut","",1000,0,10,1000,0,10);
  hInvMassAntiLambda           =new TH1F("hInvMassAntiLambda","",10000,0,10);
  //invariant mass lambda
  hInvMassXidecayLambda        =new TH1F("hInvMassXidecayLambda","",10000,0,10);
  hInvMassXidecayantiLambda    =new TH1F("hInvMassXidecayantiLambda","",10000,0,10);
  //invariant mass xi                                                                           
  hInvMassXi                   =new TH1F("hInvMassXi","",10000,0,10); //xiのmassが0.001GeV/c オーダーで決定されてるから10000 bin
  hInvMassposXi                =new TH1F("hInvMassposXi","",10000,0,10);
  XiCPA                        =new TH1F("XiCPA","",100,-1,1);
  XiTranverseradius            =new TH1F("XiTranverseradius","",1000,0,300);
  dcaXiDaughters               =new TH1F("dcaXiDaughters","",100,0,10);
  fPtcascade                   =new TH1F("fPtcascade","",1000,0,10);
  fPtv0casc                    =new TH1F("fPtv0casc","",1000,0,10);
  fPtBachelor                  =new TH1F("fPtBachelor","",1000,0,10);
  fPhiXidecayLambda            =new TH1F("fPhiXidecayLambda","",1000,0,10);
  fPtcascade_xip               =new TH1F("fPtcascade_xip","",1000,0,10);
  fPtv0casc_xip                =new TH1F("fPtv0casc_xip","",1000,0,10);
  fPtBachelor_xip              =new TH1F("fPtBachelor_xip","",1000,0,10);
  fPhiXidecayLambda_xip        =new TH1F("fPhiXidecayLambda_xip","",1000,0,10);
  //invariant mass dibaryon
  hInvMassProtonXi             =new TH2F("hInvMassProtonXi","",1000,0,10,1000,0,10);
  hInvMassantiProtonXi         =new TH2F("hInvMassantiProtonXi","",1000,0,10,1000,0,10);
  hInvMassLambdaXi             =new TH2F("hInvMassLambdaXi","",1000,0,10,1000,0,10);
  hInvMassAntiLambdaXip        =new TH2F("hInvMassAntiLambdaXip","",1000,0,10,1000,0,10);
  //decay rejection                                                                     
  fPt_allLambda                =new TH1F("fPt_allLambda","",1000,0,10);
  fPt_xidecaylambda            =new TH1F("fPt_xidecaylambda","",1000,0,10);
  fPt_allAntiLambda            =new TH1F("fPt_allAntiLambda","",1000,0,10);
  fPt_xidecayAntiLambda        =new TH1F("fPt_xidecayAntiLambda","",1000,0,10);
  //decay rejection pair                                                                        
  hInvMassLambdaLambda_all                =new TH2F("hInvMassLambdaLambda_all","",1000,0,10,1000,0,10);
  hInvMassLambdaLambda_onlyprompt         =new TH2F("hInvMassLambdaLambda_onlyprompt","",1000,0,10,1000,0,10);
  hInvMassAntiLambdaAntiLambda_all        =new TH2F("hInvMassAntiLambdaAntiLambda_all","",1000,0,10,1000,0,10);
  hInvMassAntiLambdaAntiLambda_onlyprompt =new TH2F("hInvMassAntiLambdaAntiLambda_onlyprompt","",1000,0,10,1000,0,10);
  //eventmixing                                                                                              
  //hInvMassLambdaLambda_combinatorial      =new TH2F("hInvMassLambdaLambda_eventmixing","",1000,0,10,1000,0,10);
  hInvMassLambdaLambda_evtpool            =new TH2F("hInvMassLambdaLambda_evtpool","",1000,0,10,1000,0,10);
  hInvMassAntiLambdaAntiLambda_evtpool    =new TH2F("hInvMassAntiLambdaAntiLambda_evtpool","",1000,0,10,1000,0,10);
  hInvMassProtonXi_evtpool                =new TH2F("hInvMassProtonXi_evtpool","",1000,0,10,1000,0,10);
  hInvMassAntiProtonXi_evtpool            =new TH2F("hInvMassAntiProtonXi_evtpool","",1000,0,10,1000,0,10);
  hInvMassLambdaXi_evtpool                =new TH2F("hInvMassLambdaXi_evtpool","",1000,0,10,1000,0,10);
  hInvMassAntiLambdaXi_evtpool            =new TH2F("hInvMassAntiLambdaXi_evtpool","",1000,0,10,1000,0,10);
  //vs relative momentum                                                                                    
  hInvMassLambdaLambda_relmomentum        =new TH2F("hInvMassLambdaLambda_relmomentum","",1000,0,10,1000,0,10);
  hInvMassLambdaLambda_openangle          =new TH2F("hInvMassLambdaLambda_openangle","",1800,0,180,1000,0,10);
  hInvMassAntiLambdaAntiLambda_relmomentum=new TH2F("hInvMassAntiLambdaAntiLambda_relmomentum","",1000,0,10,1000,0,10);
  hInvMassAntiLambdaAntiLambda_openangle  =new TH2F("hInvMassAntiLambdaAntiLambda_openangle","",1800,0,180,1000,0,10);
  hAngle_relmom                           =new TH2F("hAngle_relmom","",1000,0,10,3600,0,360);
  hInvMassProtonXi_relmomentum            =new TH2F("hInvMassProtonXi_relmomentum","",1000,0,10,1000,0,10);
  hInvMassProtonXi_openangle              =new TH2F("hInvMassProtonXi_openangle","",1800,0,180,1000,0,10);
  hInvMassantiProtonXi_relmomentum        =new TH2F("hInvMassantiProtonXi_relmomentum","",1000,0,10,1000,0,10);
  hInvMassantiProtonXi_openangle          =new TH2F("hInvMassantiProtonXi_openangle","",1800,0,180,1000,0,10);
  hInvMassLambdaXi_relmomentum            =new TH2F("hInvMassLambdaXi_relmomentum","",1000,0,10,1000,0,10);
  hInvMassLambdaXi_openangle              =new TH2F("hInvMassLambdaXi_openangle","",1800,0,180,1000,0,10);
  hInvMassAntiLambdaXi_relmomentum        =new TH2F("hInvMassAntiLambdaXi_relmomentum","",1000,0,10,1000,0,10);
  hInvMassAntiLambdaXi_openangle          =new TH2F("hInvMassAntiLambdaXi_openangle","",1800,0,180,1000,0,10);
  //pt correlation
  hPtcorrelation_lamlam         =new TH2F("hPtcorrelation_lamlam","",1000,0,10,1000,0,10);
  hPtcorrelation_pxi            =new TH2F("hPtcorrelation_pxi","",1000,0,10,1000,0,10);
  hPtcorrelation_lamxi          =new TH2F("hPtcorrelation_lamxi","",1000,0,10,1000,0,10);
  //lowmass                                                                                             
  hInvMassProtonXi_lowmassx     =new TH2F("hInvMassProtonXi_lowmassx","",1000,0,10,1000,0,10);
  hInvMassProtonXi_lowmassy     =new TH2F("hInvMassProtonXi_lowmassy","",1000,0,10,1000,0,10);
  hInvMassProtonXi_lowmassz     =new TH2F("hInvMassProtonXi_lowmassz","",1000,0,10,1000,0,10);
  hInvMassantiProtonXi_lowmassx =new TH2F("hInvMassantiProtonXi_lowmassx","",1000,0,10,1000,0,10);
  hInvMassantiProtonXi_lowmassy =new TH2F("hInvMassantiProtonXi_lowmassy","",1000,0,10,1000,0,10);
  hInvMassantiProtonXi_lowmassz =new TH2F("hInvMassantiProtonXi_lowmassz","",1000,0,10,1000,0,10);
  //FG LambdaLambda 1MeV Bin 30 MeV                                                                                               
  hInvMassLambdaLambda_1  =new TH1F("hInvMassLambdaLambda_1","FG LambdaLambda(1,4)",10000,0,10);
  hInvMassLambdaLambda_2  =new TH1F("hInvMassLambdaLambda_2","FG LambdaLambda(5,8)",10000,0,10);
  hInvMassLambdaLambda_3  =new TH1F("hInvMassLambdaLambda_3","FG LambdaLambda(9,12)",10000,0,10);
  hInvMassLambdaLambda_4  =new TH1F("hInvMassLambdaLambda_4","FG LambdaLambda(13,16)",10000,0,10);
  hInvMassLambdaLambda_5  =new TH1F("hInvMassLambdaLambda_5","FG LambdaLambda(17,20)",10000,0,10);
  hInvMassLambdaLambda_6  =new TH1F("hInvMassLambdaLambda_6","FG LambdaLambda(21,24)",10000,0,10);
  hInvMassLambdaLambda_7  =new TH1F("hInvMassLambdaLambda_7","FG LambdaLambda(25,28)",10000,0,10);
  hInvMassLambdaLambda_8  =new TH1F("hInvMassLambdaLambda_8","FG LambdaLambda(29,32)",10000,0,10);
  hInvMassLambdaLambda_9  =new TH1F("hInvMassLambdaLambda_9","FG LambdaLambda(33,36)",10000,0,10);
  hInvMassLambdaLambda_10 =new TH1F("hInvMassLambdaLambda_10","FG LambdaLambda(37,40)",10000,0,10);
  hInvMassLambdaLambda_11 =new TH1F("hInvMassLambdaLambda_11","FG LambdaLambda(41,50)",10000,0,10);
  hInvMassLambdaLambda_12 =new TH1F("hInvMassLambdaLambda_12","FG LambdaLambda(51,60)",10000,0,10);
  hInvMassLambdaLambda_13 =new TH1F("hInvMassLambdaLambda_13","FG LambdaLambda(61,70)",10000,0,10);
  hInvMassLambdaLambda_14 =new TH1F("hInvMassLambdaLambda_14","FG LambdaLambda(71,80)",10000,0,10);
  hInvMassLambdaLambda_15 =new TH1F("hInvMassLambdaLambda_15","FG LambdaLambda(81,90)",10000,0,10);
  hInvMassLambdaLambda_16 =new TH1F("hInvMassLambdaLambda_16","FG LambdaLambda(91,100)",10000,0,10);
  hInvMassLambdaLambda_17 =new TH1F("hInvMassLambdaLambda_17","FG LambdaLambda(100>)",10000,0,10);
  
  Int_t n_hInvMassLL = 17;
  TH1F hInvMassLambdaLambda[n_hInvMassLL];
  for (Int_t i = 0; i < 10; i++) hInvMassLambdaLambda[i] = new TH1F(Form("hInvMassLambdaLambda_%d", i + 1), Form("FG LambdaLambda(%d,%d)", i * 4 + 1, (i + 1) * 4), 10000, 0, 10);
  for (Int_t i = 10; i < 16; i++) hInvMassLambdaLambda[i] = new TH1F(Form("hInvMassLambdaLambda_%d", i + 1), Form("FG LambdaLambda(%d,%d)", (i - 10) * 10 + 41, (i - 9) * 10 + 40), 10000, 0, 10);
  hInvMassLambdaLambda[16] = new TH1F("hInvMassLambdaLambda_16", "FG LambdaLambda(100>)", 10000, 0, 10);
  
  Int_t n_hInvMassLL = 17;
  Char_t multiplicity[n_hInvMassLL][10];
  for (Int_t i = 0; i < 10; i++) sprintf(multiplicity[i], "(%d,%d)", i * 4 + 1, (i + 1) * 4);
  for (Int_t i = 10; i < 16; i++) sprintf(multiplicity[i], "(%d,%d)", (i - 10) * 10 + 41, (i - 9) * 10 + 40);
  sprintf(multiplicity[i], "(100>)");
  
  TH1F hInvMassLambdaLambda[n_hInvMassLL];
  for (Int_t i = 0; i < n_hInvMassLL; i++) hInvMassLambdaLambda[i] = new TH1F(Form("hInvMassLambdaLambda_%d", i + 1), Form("FG LambdaLambda%s", multiplivity[i]), 10000, 0, 10);
  
  hInvMassAntiLambdaAntiLambda_1  =new TH1F("hInvMassAntiLambdaAntiLambda_1","FG AntiLambdaAntiLambda(1,4)",10000,0,10);
  hInvMassAntiLambdaAntiLambda_2  =new TH1F("hInvMassAntiLambdaAntiLambda_2","FG AntiLambdaAntiLambda(5,8)",10000,0,10);
  hInvMassAntiLambdaAntiLambda_3  =new TH1F("hInvMassAntiLambdaAntiLambda_3","FG AntiLambdaAntiLambda(9,12)",10000,0,10);
  hInvMassAntiLambdaAntiLambda_4  =new TH1F("hInvMassAntiLambdaAntiLambda_4","FG AntiLambdaAntiLambda(13,16)",10000,0,10);
  hInvMassAntiLambdaAntiLambda_5  =new TH1F("hInvMassAntiLambdaAntiLambda_5","FG AntiLambdaAntiLambda(17,20)",10000,0,10);
  hInvMassAntiLambdaAntiLambda_6  =new TH1F("hInvMassAntiLambdaAntiLambda_6","FG AntiLambdaAntiLambda(21,24)",10000,0,10);
  hInvMassAntiLambdaAntiLambda_7  =new TH1F("hInvMassAntiLambdaAntiLambda_7","FG AntiLambdaAntiLambda(25,28)",10000,0,10);
  hInvMassAntiLambdaAntiLambda_8  =new TH1F("hInvMassAntiLambdaAntiLambda_8","FG AntiLambdaAntiLambda(29,32)",10000,0,10);
  hInvMassAntiLambdaAntiLambda_9  =new TH1F("hInvMassAntiLambdaAntiLambda_9","FG AntiLambdaAntiLambda(33,36)",10000,0,10);
  hInvMassAntiLambdaAntiLambda_10 =new TH1F("hInvMassAntiLambdaAntiLambda_10","FG AntiLambdaAntiLambda(37,40)",10000,0,10);
  hInvMassAntiLambdaAntiLambda_11 =new TH1F("hInvMassAntiLambdaAntiLambda_11","FG AntiLambdaAntiLambda(41,50)",10000,0,10);
  hInvMassAntiLambdaAntiLambda_12 =new TH1F("hInvMassAntiLambdaAntiLambda_12","FG AntiLambdaAntiLambda(51,60)",10000,0,10);
  hInvMassAntiLambdaAntiLambda_13 =new TH1F("hInvMassAntiLambdaAntiLambda_13","FG AntiLambdaAntiLambda(61,70)",10000,0,10);
  hInvMassAntiLambdaAntiLambda_14 =new TH1F("hInvMassAntiLambdaAntiLambda_14","FG AntiLambdaAntiLambda(71,80)",10000,0,10);
  hInvMassAntiLambdaAntiLambda_15 =new TH1F("hInvMassAntiLambdaAntiLambda_15","FG AntiLambdaAntiLambda(81,90)",10000,0,10);
  hInvMassAntiLambdaAntiLambda_16 =new TH1F("hInvMassAntiLambdaAntiLambda_16","FG AntiLambdaAntiLambda(91,100)",10000,0,10);
  hInvMassAntiLambdaAntiLambda_17 =new TH1F("hInvMassAntiLambdaAntiLambda_17","FG AntiLambdaAntiLambda(100>)",10000,0,10);
  //FG ProtonXi
  hInvMassProtonXi_1  =new TH1F("hInvMassProtonXi_1","FG ProtonXi(1,4)",10000,0,10);
  hInvMassProtonXi_2  =new TH1F("hInvMassProtonXi_2","FG ProtonXi(5,8)",10000,0,10);
  hInvMassProtonXi_3  =new TH1F("hInvMassProtonXi_3","FG ProtonXi(9,12)",10000,0,10);
  hInvMassProtonXi_4  =new TH1F("hInvMassProtonXi_4","FG ProtonXi(13,16)",10000,0,10);
  hInvMassProtonXi_5  =new TH1F("hInvMassProtonXi_5","FG ProtonXi(17,20)",10000,0,10);
  hInvMassProtonXi_6  =new TH1F("hInvMassProtonXi_6","FG ProtonXi(21,24)",10000,0,10);
  hInvMassProtonXi_7  =new TH1F("hInvMassProtonXi_7","FG ProtonXi(25,28)",10000,0,10);
  hInvMassProtonXi_8  =new TH1F("hInvMassProtonXi_8","FG ProtonXi(29,32)",10000,0,10);
  hInvMassProtonXi_9  =new TH1F("hInvMassProtonXi_9","FG ProtonXi(33,36)",10000,0,10);
  hInvMassProtonXi_10 =new TH1F("hInvMassProtonXi_10","FG ProtonXi(37,40)",10000,0,10);
  hInvMassProtonXi_11 =new TH1F("hInvMassProtonXi_11","FG ProtonXi(41,50)",10000,0,10);
  hInvMassProtonXi_12 =new TH1F("hInvMassProtonXi_12","FG ProtonXi(51,60)",10000,0,10);
  hInvMassProtonXi_13 =new TH1F("hInvMassProtonXi_13","FG ProtonXi(61,70)",10000,0,10);
  hInvMassProtonXi_14 =new TH1F("hInvMassProtonXi_14","FG ProtonXi(71,80)",10000,0,10);
  hInvMassProtonXi_15 =new TH1F("hInvMassProtonXi_15","FG ProtonXi(81,90)",10000,0,10);
  hInvMassProtonXi_16 =new TH1F("hInvMassProtonXi_16","FG ProtonXi(91,100)",10000,0,10);
  hInvMassProtonXi_17 =new TH1F("hInvMassProtonXi_17","FG ProtonXi(100>)",10000,0,10);
  hInvMassAntiProtonXi_1  =new TH1F("hInvMassAntiProtonXi_1","FG AntiProtonXi(1,4)",10000,0,10);
  hInvMassAntiProtonXi_2  =new TH1F("hInvMassAntiProtonXi_2","FG AntiProtonXi(5,8)",10000,0,10);
  hInvMassAntiProtonXi_3  =new TH1F("hInvMassAntiProtonXi_3","FG AntiProtonXi(9,12)",10000,0,10);
  hInvMassAntiProtonXi_4  =new TH1F("hInvMassAntiProtonXi_4","FG AntiProtonXi(13,16)",10000,0,10);
  hInvMassAntiProtonXi_5  =new TH1F("hInvMassAntiProtonXi_5","FG AntiProtonXi(17,20)",10000,0,10);
  hInvMassAntiProtonXi_6  =new TH1F("hInvMassAntiProtonXi_6","FG AntiProtonXi(21,24)",10000,0,10);
  hInvMassAntiProtonXi_7  =new TH1F("hInvMassAntiProtonXi_7","FG AntiProtonXi(25,28)",10000,0,10);
  hInvMassAntiProtonXi_8  =new TH1F("hInvMassAntiProtonXi_8","FG AntiProtonXi(29,32)",10000,0,10);
  hInvMassAntiProtonXi_9  =new TH1F("hInvMassAntiProtonXi_9","FG AntiProtonXi(33,36)",10000,0,10);
  hInvMassAntiProtonXi_10 =new TH1F("hInvMassAntiProtonXi_10","FG AntiProtonXi(37,40)",10000,0,10);
  hInvMassAntiProtonXi_11 =new TH1F("hInvMassAntiProtonXi_11","FG AntiProtonXi(41,50)",10000,0,10);
  hInvMassAntiProtonXi_12 =new TH1F("hInvMassAntiProtonXi_12","FG AntiProtonXi(51,60)",10000,0,10);
  hInvMassAntiProtonXi_13 =new TH1F("hInvMassAntiProtonXi_13","FG AntiProtonXi(61,70)",10000,0,10);
  hInvMassAntiProtonXi_14 =new TH1F("hInvMassAntiProtonXi_14","FG AntiProtonXi(71,80)",10000,0,10);
  hInvMassAntiProtonXi_15 =new TH1F("hInvMassAntiProtonXi_15","FG AntiProtonXi(81,90)",10000,0,10);
  hInvMassAntiProtonXi_16 =new TH1F("hInvMassAntiProtonXi_16","FG AntiProtonXi(91,100)",10000,0,10);
  hInvMassAntiProtonXi_17 =new TH1F("hInvMassAntiProtonXi_17","FG AntiProtonXi(100>)",10000,0,10);                                                                    
  //FG LambdaXi
  hInvMassLambdaXi_1  =new TH1F("hInvMassLambdaXi_1","FG LambdaXi(1,5)",10000,0,10);
  hInvMassLambdaXi_2  =new TH1F("hInvMassLambdaXi_2","FG LambdaXi(6,10)",10000,0,10);
  hInvMassLambdaXi_3  =new TH1F("hInvMassLambdaXi_3","FG LambdaXi(11,15)",10000,0,10);
  hInvMassLambdaXi_4  =new TH1F("hInvMassLambdaXi_4","FG LambdaXi(16,20)",10000,0,10);
  hInvMassLambdaXi_5  =new TH1F("hInvMassLambdaXi_5","FG LambdaXi(21,25)",10000,0,10);
  hInvMassLambdaXi_6  =new TH1F("hInvMassLambdaXi_6","FG LambdaXi(26,30)",10000,0,10);
  hInvMassLambdaXi_7  =new TH1F("hInvMassLambdaXi_7","FG LambdaXi(31,35)",10000,0,10);
  hInvMassLambdaXi_8  =new TH1F("hInvMassLambdaXi_8","FG LambdaXi(36,40)",10000,0,10);
  hInvMassLambdaXi_9  =new TH1F("hInvMassLambdaXi_9","FG LambdaXi(41,50)",10000,0,10);
  hInvMassLambdaXi_10 =new TH1F("hInvMassLambdaXi_10","FG LambdaXi(51,60)",10000,0,10);
  hInvMassLambdaXi_11 =new TH1F("hInvMassLambdaXi_11","FG LambdaXi(61,70)",10000,0,10);
  hInvMassLambdaXi_12 =new TH1F("hInvMassLambdaXi_12","FG LambdaXi(71,80)",10000,0,10);
  hInvMassLambdaXi_13 =new TH1F("hInvMassLambdaXi_13","FG LambdaXi(81,90)",10000,0,10);
  hInvMassLambdaXi_14 =new TH1F("hInvMassLambdaXi_14","FG LambdaXi(91,100)",10000,0,10);
  hInvMassLambdaXi_15 =new TH1F("hInvMassLambdaXi_15","FG LambdaXi(100>)",10000,0,10);
  hInvMassAntiLambdaXi_1  =new TH1F("hInvMassAntiLambdaXi_1","FG AntiLambdaXi(1,5)",10000,0,10);
  hInvMassAntiLambdaXi_2  =new TH1F("hInvMassAntiLambdaXi_2","FG AntiLambdaXi(6,10)",10000,0,10);
  hInvMassAntiLambdaXi_3  =new TH1F("hInvMassAntiLambdaXi_3","FG AntiLambdaXi(11,15)",10000,0,10);
  hInvMassAntiLambdaXi_4  =new TH1F("hInvMassAntiLambdaXi_4","FG AntiLambdaXi(16,20)",10000,0,10);
  hInvMassAntiLambdaXi_5  =new TH1F("hInvMassAntiLambdaXi_5","FG AntiLambdaXi(21,25)",10000,0,10);
  hInvMassAntiLambdaXi_6  =new TH1F("hInvMassAntiLambdaXi_6","FG AntiLambdaXi(26,30)",10000,0,10);
  hInvMassAntiLambdaXi_7  =new TH1F("hInvMassAntiLambdaXi_7","FG AntiLambdaXi(31,35)",10000,0,10);
  hInvMassAntiLambdaXi_8  =new TH1F("hInvMassAntiLambdaXi_8","FG AntiLambdaXi(36,40)",10000,0,10);
  hInvMassAntiLambdaXi_9  =new TH1F("hInvMassAntiLambdaXi_9","FG AntiLambdaXi(41,50)",10000,0,10);
  hInvMassAntiLambdaXi_10 =new TH1F("hInvMassAntiLambdaXi_10","FG AntiLambdaXi(51,60)",10000,0,10);
  hInvMassAntiLambdaXi_11 =new TH1F("hInvMassAntiLambdaXi_11","FG AntiLambdaXi(61,70)",10000,0,10);
  hInvMassAntiLambdaXi_12 =new TH1F("hInvMassAntiLambdaXi_12","FG AntiLambdaXi(71,80)",10000,0,10);
  hInvMassAntiLambdaXi_13 =new TH1F("hInvMassAntiLambdaXi_13","FG AntiLambdaXi(81,90)",10000,0,10);
  hInvMassAntiLambdaXi_14 =new TH1F("hInvMassAntiLambdaXi_14","FG AntiLambdaXi(91,100)",10000,0,10);
  hInvMassAntiLambdaXi_15 =new TH1F("hInvMassAntiLambdaXi_15","FG AntiLambdaXi(100>)",10000,0,10);
  //BG LambdaLambda
  hInvMassLambdaLambda_evtpool1  =new TH1F("hInvMassLambdaLambda_evtpool1","LambdaLambda(1,4)",10000,0,10);
  hInvMassLambdaLambda_evtpool2  =new TH1F("hInvMassLambdaLambda_evtpool2","LambdaLambda(5,8)",10000,0,10);
  hInvMassLambdaLambda_evtpool3  =new TH1F("hInvMassLambdaLambda_evtpool3","LambdaLambda(9,12)",10000,0,10);
  hInvMassLambdaLambda_evtpool4  =new TH1F("hInvMassLambdaLambda_evtpool4","LambdaLambda(13,16)",10000,0,10);
  hInvMassLambdaLambda_evtpool5  =new TH1F("hInvMassLambdaLambda_evtpool5","LambdaLambda(17,20)",10000,0,10);
  hInvMassLambdaLambda_evtpool6  =new TH1F("hInvMassLambdaLambda_evtpool6","LambdaLambda(21,24)",10000,0,10);
  hInvMassLambdaLambda_evtpool7  =new TH1F("hInvMassLambdaLambda_evtpool7","LambdaLambda(25,28)",10000,0,10);
  hInvMassLambdaLambda_evtpool8  =new TH1F("hInvMassLambdaLambda_evtpool8","LambdaLambda(29,32)",10000,0,10);
  hInvMassLambdaLambda_evtpool9  =new TH1F("hInvMassLambdaLambda_evtpool9","LambdaLambda(33,36)",10000,0,10);
  hInvMassLambdaLambda_evtpool10 =new TH1F("hInvMassLambdaLambda_evtpool10","LambdaLambda(37,40)",10000,0,10);
  hInvMassLambdaLambda_evtpool11 =new TH1F("hInvMassLambdaLambda_evtpool11","LambdaLambda(41,50)",10000,0,10);
  hInvMassLambdaLambda_evtpool12 =new TH1F("hInvMassLambdaLambda_evtpool12","LambdaLambda(51,60)",10000,0,10);
  hInvMassLambdaLambda_evtpool13 =new TH1F("hInvMassLambdaLambda_evtpool13","LambdaLambda(61,70)",10000,0,10);
  hInvMassLambdaLambda_evtpool14 =new TH1F("hInvMassLambdaLambda_evtpool14","LambdaLambda(71,80)",10000,0,10);
  hInvMassLambdaLambda_evtpool15 =new TH1F("hInvMassLambdaLambda_evtpool15","LambdaLambda(81,90)",10000,0,10);
  hInvMassLambdaLambda_evtpool16 =new TH1F("hInvMassLambdaLambda_evtpool16","LambdaLambda(91,100)",10000,0,10);
  hInvMassLambdaLambda_evtpool17 =new TH1F("hInvMassLambdaLambda_evtpool17","LambdaLambda(100>)",10000,0,10);
  hInvMassAntiLambdaAntiLambda_evtpool1  =new TH1F("hInvMassAntiLambdaAntiLambda_evtpool1","AntiLambdaAntiLambda(1,4)",10000,0,10);
  hInvMassAntiLambdaAntiLambda_evtpool2  =new TH1F("hInvMassAntiLambdaAntiLambda_evtpool2","AntiLambdaAntiLambda(5,8)",10000,0,10);
  hInvMassAntiLambdaAntiLambda_evtpool3  =new TH1F("hInvMassAntiLambdaAntiLambda_evtpool3","AntiLambdaAntiLambda(9,12)",10000,0,10);
  hInvMassAntiLambdaAntiLambda_evtpool4  =new TH1F("hInvMassAntiLambdaAntiLambda_evtpool4","AntiLambdaAntiLambda(13,16)",10000,0,10);
  hInvMassAntiLambdaAntiLambda_evtpool5  =new TH1F("hInvMassAntiLambdaAntiLambda_evtpool5","AntiLambdaAntiLambda(17,20)",10000,0,10);
  hInvMassAntiLambdaAntiLambda_evtpool6  =new TH1F("hInvMassAntiLambdaAntiLambda_evtpool6","AntiLambdaAntiLambda(21,24)",10000,0,10);
  hInvMassAntiLambdaAntiLambda_evtpool7  =new TH1F("hInvMassAntiLambdaAntiLambda_evtpool7","AntiLambdaAntiLambda(25,28)",10000,0,10);
  hInvMassAntiLambdaAntiLambda_evtpool8  =new TH1F("hInvMassAntiLambdaAntiLambda_evtpool8","AntiLambdaAntiLambda(29,32)",10000,0,10);
  hInvMassAntiLambdaAntiLambda_evtpool9  =new TH1F("hInvMassAntiLambdaAntiLambda_evtpool9","AntiLambdaAntiLambda(33,36)",10000,0,10);
  hInvMassAntiLambdaAntiLambda_evtpool10 =new TH1F("hInvMassAntiLambdaAntiLambda_evtpool10","AntiLambdaAntiLambda(37,40)",10000,0,10);
  hInvMassAntiLambdaAntiLambda_evtpool11 =new TH1F("hInvMassAntiLambdaAntiLambda_evtpool11","AntiLambdaAntiLambda(41,50)",10000,0,10);
  hInvMassAntiLambdaAntiLambda_evtpool12 =new TH1F("hInvMassAntiLambdaAntiLambda_evtpool12","AntiLambdaAntiLambda(51,60)",10000,0,10);
  hInvMassAntiLambdaAntiLambda_evtpool13 =new TH1F("hInvMassAntiLambdaAntiLambda_evtpool13","AntiLambdaAntiLambda(61,70)",10000,0,10);
  hInvMassAntiLambdaAntiLambda_evtpool14 =new TH1F("hInvMassAntiLambdaAntiLambda_evtpool14","AntiLambdaAntiLambda(71,80)",10000,0,10);
  hInvMassAntiLambdaAntiLambda_evtpool15 =new TH1F("hInvMassAntiLambdaAntiLambda_evtpool15","AntiLambdaAntiLambda(81,90)",10000,0,10);
  hInvMassAntiLambdaAntiLambda_evtpool16 =new TH1F("hInvMassAntiLambdaAntiLambda_evtpool16","AntiLambdaAntiLambda(91,100)",10000,0,10);
  hInvMassAntiLambdaAntiLambda_evtpool17 =new TH1F("hInvMassAntiLambdaAntiLambda_evtpool17","AntiLambdaAntiLambda(100>)",10000,0,10);
  //BG ProtonXi
  hInvMassProtonXi_evtpool1  =new TH1F("hInvMassProtonXi_evtpool1","ProtonXi(1,4)",10000,0,10);
  hInvMassProtonXi_evtpool2  =new TH1F("hInvMassProtonXi_evtpool2","ProtonXi(5,8)",10000,0,10);
  hInvMassProtonXi_evtpool3  =new TH1F("hInvMassProtonXi_evtpool3","ProtonXi(9,12)",10000,0,10);
  hInvMassProtonXi_evtpool4  =new TH1F("hInvMassProtonXi_evtpool4","ProtonXi(13,16)",10000,0,10);
  hInvMassProtonXi_evtpool5  =new TH1F("hInvMassProtonXi_evtpool5","ProtonXi(17,20)",10000,0,10);
  hInvMassProtonXi_evtpool6  =new TH1F("hInvMassProtonXi_evtpool6","ProtonXi(21,24)",10000,0,10);
  hInvMassProtonXi_evtpool7  =new TH1F("hInvMassProtonXi_evtpool7","ProtonXi(25,28)",10000,0,10);
  hInvMassProtonXi_evtpool8  =new TH1F("hInvMassProtonXi_evtpool8","ProtonXi(29,32)",10000,0,10);
  hInvMassProtonXi_evtpool9  =new TH1F("hInvMassProtonXi_evtpool9","ProtonXi(33,36)",10000,0,10);
  hInvMassProtonXi_evtpool10 =new TH1F("hInvMassProtonXi_evtpool10","ProtonXi(37,40)",10000,0,10);
  hInvMassProtonXi_evtpool11 =new TH1F("hInvMassProtonXi_evtpool11","ProtonXi(41,50)",10000,0,10);
  hInvMassProtonXi_evtpool12 =new TH1F("hInvMassProtonXi_evtpool12","ProtonXi(51,60)",10000,0,10);
  hInvMassProtonXi_evtpool13 =new TH1F("hInvMassProtonXi_evtpool13","ProtonXi(61,70)",10000,0,10);
  hInvMassProtonXi_evtpool14 =new TH1F("hInvMassProtonXi_evtpool14","ProtonXi(71,80)",10000,0,10);
  hInvMassProtonXi_evtpool15 =new TH1F("hInvMassProtonXi_evtpool15","ProtonXi(81,90)",10000,0,10);
  hInvMassProtonXi_evtpool16 =new TH1F("hInvMassProtonXi_evtpool16","ProtonXi(91,100)",10000,0,10);
  hInvMassProtonXi_evtpool17 =new TH1F("hInvMassProtonXi_evtpool17","ProtonXi(100>)",10000,0,10);
  hInvMassAntiProtonXi_evtpool1  =new TH1F("hInvMassAntiProtonXi_evtpool1","AntiProtonXi(1,4)",10000,0,10);
  hInvMassAntiProtonXi_evtpool2  =new TH1F("hInvMassAntiProtonXi_evtpool2","AntiProtonXi(5,8)",10000,0,10);
  hInvMassAntiProtonXi_evtpool3  =new TH1F("hInvMassAntiProtonXi_evtpool3","AntiProtonXi(9,12)",10000,0,10);
  hInvMassAntiProtonXi_evtpool4  =new TH1F("hInvMassAntiProtonXi_evtpool4","AntiProtonXi(13,16)",10000,0,10);
  hInvMassAntiProtonXi_evtpool5  =new TH1F("hInvMassAntiProtonXi_evtpool5","AntiProtonXi(17,20)",10000,0,10);
  hInvMassAntiProtonXi_evtpool6  =new TH1F("hInvMassAntiProtonXi_evtpool6","AntiProtonXi(21,24)",10000,0,10);
  hInvMassAntiProtonXi_evtpool7  =new TH1F("hInvMassAntiProtonXi_evtpool7","AntiProtonXi(25,28)",10000,0,10);
  hInvMassAntiProtonXi_evtpool8  =new TH1F("hInvMassAntiProtonXi_evtpool8","AntiProtonXi(29,32)",10000,0,10);
  hInvMassAntiProtonXi_evtpool9  =new TH1F("hInvMassAntiProtonXi_evtpool9","AntiProtonXi(33,36)",10000,0,10);
  hInvMassAntiProtonXi_evtpool10 =new TH1F("hInvMassAntiProtonXi_evtpool10","AntiProtonXi(37,40)",10000,0,10);
  hInvMassAntiProtonXi_evtpool11 =new TH1F("hInvMassAntiProtonXi_evtpool11","AntiProtonXi(41,50)",10000,0,10);
  hInvMassAntiProtonXi_evtpool12 =new TH1F("hInvMassAntiProtonXi_evtpool12","AntiProtonXi(51,60)",10000,0,10);
  hInvMassAntiProtonXi_evtpool13 =new TH1F("hInvMassAntiProtonXi_evtpool13","AntiProtonXi(61,70)",10000,0,10);
  hInvMassAntiProtonXi_evtpool14 =new TH1F("hInvMassAntiProtonXi_evtpool14","AntiProtonXi(71,80)",10000,0,10);
  hInvMassAntiProtonXi_evtpool15 =new TH1F("hInvMassAntiProtonXi_evtpool15","AntiProtonXi(81,90)",10000,0,10);
  hInvMassAntiProtonXi_evtpool16 =new TH1F("hInvMassAntiProtonXi_evtpool16","AntiProtonXi(91,100)",10000,0,10);
  hInvMassAntiProtonXi_evtpool17 =new TH1F("hInvMassAntiProtonXi_evtpool17","AntiProtonXi(100>)",10000,0,10);
  //BG LambdaXi
  hInvMassLambdaXi_evtpool1  =new TH1F("hInvMassLambdaXi_evtpool1","LambdaXi(1,5)",10000,0,10);
  hInvMassLambdaXi_evtpool2  =new TH1F("hInvMassLambdaXi_evtpool2","LambdaXi(6,10)",10000,0,10);
  hInvMassLambdaXi_evtpool3  =new TH1F("hInvMassLambdaXi_evtpool3","LambdaXi(11,15)",10000,0,10);
  hInvMassLambdaXi_evtpool4  =new TH1F("hInvMassLambdaXi_evtpool4","LambdaXi(16,20)",10000,0,10);
  hInvMassLambdaXi_evtpool5  =new TH1F("hInvMassLambdaXi_evtpool5","LambdaXi(21,25)",10000,0,10);
  hInvMassLambdaXi_evtpool6  =new TH1F("hInvMassLambdaXi_evtpool6","LambdaXi(26,30)",10000,0,10);
  hInvMassLambdaXi_evtpool7  =new TH1F("hInvMassLambdaXi_evtpool7","LambdaXi(31,35)",10000,0,10);
  hInvMassLambdaXi_evtpool8  =new TH1F("hInvMassLambdaXi_evtpool8","LambdaXi(36,40)",10000,0,10);
  hInvMassLambdaXi_evtpool9  =new TH1F("hInvMassLambdaXi_evtpool9","LambdaXi(41,50)",10000,0,10);
  hInvMassLambdaXi_evtpool10 =new TH1F("hInvMassLambdaXi_evtpool10","LambdaXi(51,60)",10000,0,10);
  hInvMassLambdaXi_evtpool11 =new TH1F("hInvMassLambdaXi_evtpool11","LambdaXi(61,70)",10000,0,10);
  hInvMassLambdaXi_evtpool12 =new TH1F("hInvMassLambdaXi_evtpool12","LambdaXi(71,80)",10000,0,10);
  hInvMassLambdaXi_evtpool13 =new TH1F("hInvMassLambdaXi_evtpool13","LambdaXi(81,90)",10000,0,10);
  hInvMassLambdaXi_evtpool14 =new TH1F("hInvMassLambdaXi_evtpool14","LambdaXi(91,100)",10000,0,10);
  hInvMassLambdaXi_evtpool15 =new TH1F("hInvMassLambdaXi_evtpool15","LambdaXi(100>)",10000,0,10);
  hInvMassAntiLambdaXi_evtpool1  =new TH1F("hInvMassAntiLambdaXi_evtpool1","AntiLambdaXi(1,5)",10000,0,10);
  hInvMassAntiLambdaXi_evtpool2  =new TH1F("hInvMassAntiLambdaXi_evtpool2","AntiLambdaXi(6,10)",10000,0,10);
  hInvMassAntiLambdaXi_evtpool3  =new TH1F("hInvMassAntiLambdaXi_evtpool3","AntiLambdaXi(11,15)",10000,0,10);
  hInvMassAntiLambdaXi_evtpool4  =new TH1F("hInvMassAntiLambdaXi_evtpool4","AntiLambdaXi(16,20)",10000,0,10);
  hInvMassAntiLambdaXi_evtpool5  =new TH1F("hInvMassAntiLambdaXi_evtpool5","AntiLambdaXi(21,25)",10000,0,10);
  hInvMassAntiLambdaXi_evtpool6  =new TH1F("hInvMassAntiLambdaXi_evtpool6","AntiLambdaXi(26,30)",10000,0,10);
  hInvMassAntiLambdaXi_evtpool7  =new TH1F("hInvMassAntiLambdaXi_evtpool7","AntiLambdaXi(31,35)",10000,0,10);
  hInvMassAntiLambdaXi_evtpool8  =new TH1F("hInvMassAntiLambdaXi_evtpool8","AntiLambdaXi(36,40)",10000,0,10);
  hInvMassAntiLambdaXi_evtpool9  =new TH1F("hInvMassAntiLambdaXi_evtpool9","AntiLambdaXi(41,50)",10000,0,10);
  hInvMassAntiLambdaXi_evtpool10 =new TH1F("hInvMassAntiLambdaXi_evtpool10","AntiLambdaXi(51,60)",10000,0,10);
  hInvMassAntiLambdaXi_evtpool11 =new TH1F("hInvMassAntiLambdaXi_evtpool11","AntiLambdaXi(61,70)",10000,0,10);
  hInvMassAntiLambdaXi_evtpool12 =new TH1F("hInvMassAntiLambdaXi_evtpool12","AntiLambdaXi(71,80)",10000,0,10);
  hInvMassAntiLambdaXi_evtpool13 =new TH1F("hInvMassAntiLambdaXi_evtpool13","AntiLambdaXi(81,90)",10000,0,10);
  hInvMassAntiLambdaXi_evtpool14 =new TH1F("hInvMassAntiLambdaXi_evtpool14","AntiLambdaXi(91,100)",10000,0,10);
  hInvMassAntiLambdaXi_evtpool15 =new TH1F("hInvMassAntiLambdaXi_evtpool15","AntiLambdaXi(100>)",10000,0,10);

  //===============
  fOutputList->Add(fEvtCounter);
  fOutputList->Add(fEvtPassCut);
  fOutputList->Add(fEvtVtxX);
  fOutputList->Add(fEvtVtxY);
  fOutputList->Add(fEvtVtxZ);
  fOutputList->Add(fEvtVtxTrk);
  fOutputList->Add(fSPDVtxZ);
  fOutputList->Add(fSPDVtxTrk);
  fOutputList->Add(fSPDVtxCor);
  fOutputList->Add(fSPDVtxDisp);
  fOutputList->Add(fTrkEta);
  //multiplicity                                                                                                      
  fOutputList->Add(fmultiplicity);
  fOutputList->Add(fmultiplicity_lambda);
  fOutputList->Add(fmultiplicity_antilambda);
  fOutputList->Add(fmultiplicity_xim);
  fOutputList->Add(fmultiplicity_xip);
  fOutputList->Add(fmultiplicity_proton);
  fOutputList->Add(fmultiplicity_antiproton);
  fOutputList->Add(fmultiplicity_lambdalambda);
  fOutputList->Add(fmultiplicity_antilambdaantilambda);
  fOutputList->Add(fmultiplicity_protonxi);            
  fOutputList->Add(fmultiplicity_antiprotonxi);        
  fOutputList->Add(fmultiplicity_lambdaxi);            
  fOutputList->Add(fmultiplicity_antilambdaxi);        
  //proton
  fOutputList->Add(fTrkPt);
  fOutputList->Add(fTPCclus);
  fOutputList->Add(fTPCcrsR);
  fOutputList->Add(fTPCclusF);
  fOutputList->Add(fTrkDCAxy);
  fOutputList->Add(fTrkDCAz);
  fOutputList->Add(fPvsTPCsignal);
  fOutputList->Add(fPvsTOFsignal);
  //proton select
  fOutputList->Add(fPt_proton);
  fOutputList->Add(fPt_antiproton);
  fOutputList->Add(PID_cut);
  //v0 selection                                                                                    
  fOutputList->Add(fPtv0);
  fOutputList->Add(fPtv0_lambda);
  fOutputList->Add(fTransRadius);
  fOutputList->Add(fDCAdaugTov0Vtx);
  fOutputList->Add(fCPA);
  fOutputList->Add(hInvMassLambda);
  fOutputList->Add(fPtv0_antilambda);
  fOutputList->Add(fPtv0_antilambda_mass_masscut);
  fOutputList->Add(hInvMassAntiLambda);
  //invariant mass lambda
  fOutputList->Add(hInvMassXidecayLambda);
  fOutputList->Add(hInvMassXidecayantiLambda);
  //invariant mass xi          
  fOutputList->Add(hInvMassXi);
  fOutputList->Add(hInvMassposXi);
  fOutputList->Add(XiCPA);                  
  fOutputList->Add(XiTranverseradius);
  fOutputList->Add(dcaXiDaughters);
  fOutputList->Add(fPtcascade);
  fOutputList->Add(fPtv0casc);
  fOutputList->Add(fPtBachelor);
  fOutputList->Add(fPhiXidecayLambda);
  fOutputList->Add(fPtcascade_xip);
  fOutputList->Add(fPtv0casc_xip); 
  fOutputList->Add(fPtBachelor_xip);
  fOutputList->Add(fPhiXidecayLambda_xip);
  //invariant mass dibaryon
  fOutputList->Add(hInvMassProtonXi);
  fOutputList->Add(hInvMassantiProtonXi);
  fOutputList->Add(hInvMassLambdaXi);
  fOutputList->Add(hInvMassAntiLambdaXip);
  //decay rejection single
  fOutputList->Add(fPt_allLambda);
  fOutputList->Add(fPt_xidecaylambda);
  fOutputList->Add(fPt_allAntiLambda);
  fOutputList->Add(fPt_xidecayAntiLambda);
  //decay rejection pair                                                                       
  fOutputList->Add(hInvMassLambdaLambda_all);
  fOutputList->Add(hInvMassLambdaLambda_onlyprompt);
  fOutputList->Add(hInvMassAntiLambdaAntiLambda_all);
  fOutputList->Add(hInvMassAntiLambdaAntiLambda_onlyprompt);
  //eventmixing
  //fOutputList->Add(hInvMassLambdaLambda_combinatorial);
  fOutputList->Add(hInvMassLambdaLambda_evtpool);
  fOutputList->Add(hInvMassAntiLambdaAntiLambda_evtpool);
  fOutputList->Add(hInvMassProtonXi_evtpool);
  fOutputList->Add(hInvMassAntiProtonXi_evtpool);
  fOutputList->Add(hInvMassLambdaXi_evtpool);
  fOutputList->Add(hInvMassAntiLambdaXi_evtpool);
  //vs relative momentum                                                                                    
  fOutputList->Add(hInvMassLambdaLambda_relmomentum);
  fOutputList->Add(hInvMassLambdaLambda_openangle);
  fOutputList->Add(hInvMassAntiLambdaAntiLambda_relmomentum);
  fOutputList->Add(hInvMassAntiLambdaAntiLambda_openangle);
  fOutputList->Add(hAngle_relmom);
  fOutputList->Add(hInvMassProtonXi_relmomentum);
  fOutputList->Add(hInvMassProtonXi_openangle);
  fOutputList->Add(hInvMassantiProtonXi_relmomentum);
  fOutputList->Add(hInvMassantiProtonXi_openangle);
  fOutputList->Add(hInvMassLambdaXi_relmomentum);
  fOutputList->Add(hInvMassLambdaXi_openangle);
  fOutputList->Add(hInvMassAntiLambdaXi_relmomentum);
  fOutputList->Add(hInvMassAntiLambdaXi_openangle);
  //pt correlation
  fOutputList->Add(hPtcorrelation_lamlam);
  fOutputList->Add(hPtcorrelation_pxi);
  fOutputList->Add(hPtcorrelation_lamxi);
  //lowmass                                                                                             
  fOutputList->Add(hInvMassProtonXi_lowmassx);
  fOutputList->Add(hInvMassProtonXi_lowmassy);
  fOutputList->Add(hInvMassProtonXi_lowmassz);
  fOutputList->Add(hInvMassantiProtonXi_lowmassx);
  fOutputList->Add(hInvMassantiProtonXi_lowmassy);
  fOutputList->Add(hInvMassantiProtonXi_lowmassz);
  //FG LambdaLambda                                                                                                                     
  fOutputList->Add(hInvMassLambdaLambda_1);
  fOutputList->Add(hInvMassLambdaLambda_2);
  fOutputList->Add(hInvMassLambdaLambda_3);
  fOutputList->Add(hInvMassLambdaLambda_4);
  fOutputList->Add(hInvMassLambdaLambda_5);
  fOutputList->Add(hInvMassLambdaLambda_6);
  fOutputList->Add(hInvMassLambdaLambda_7);
  fOutputList->Add(hInvMassLambdaLambda_8);
  fOutputList->Add(hInvMassLambdaLambda_9);
  fOutputList->Add(hInvMassLambdaLambda_10);
  fOutputList->Add(hInvMassLambdaLambda_11);
  fOutputList->Add(hInvMassLambdaLambda_12);
  fOutputList->Add(hInvMassLambdaLambda_13);
  fOutputList->Add(hInvMassLambdaLambda_14);
  fOutputList->Add(hInvMassLambdaLambda_15);
  fOutputList->Add(hInvMassLambdaLambda_16);
  fOutputList->Add(hInvMassLambdaLambda_17);
  fOutputList->Add(hInvMassAntiLambdaAntiLambda_1);
  fOutputList->Add(hInvMassAntiLambdaAntiLambda_2);
  fOutputList->Add(hInvMassAntiLambdaAntiLambda_3);
  fOutputList->Add(hInvMassAntiLambdaAntiLambda_4);
  fOutputList->Add(hInvMassAntiLambdaAntiLambda_5);
  fOutputList->Add(hInvMassAntiLambdaAntiLambda_6);
  fOutputList->Add(hInvMassAntiLambdaAntiLambda_7);
  fOutputList->Add(hInvMassAntiLambdaAntiLambda_8);
  fOutputList->Add(hInvMassAntiLambdaAntiLambda_9);
  fOutputList->Add(hInvMassAntiLambdaAntiLambda_10);
  fOutputList->Add(hInvMassAntiLambdaAntiLambda_11);
  fOutputList->Add(hInvMassAntiLambdaAntiLambda_12);
  fOutputList->Add(hInvMassAntiLambdaAntiLambda_13);
  fOutputList->Add(hInvMassAntiLambdaAntiLambda_14);
  fOutputList->Add(hInvMassAntiLambdaAntiLambda_15);
  fOutputList->Add(hInvMassAntiLambdaAntiLambda_16);
  fOutputList->Add(hInvMassAntiLambdaAntiLambda_17);
  //FG ProtonXi                                                                                                                     
  fOutputList->Add(hInvMassProtonXi_1);
  fOutputList->Add(hInvMassProtonXi_2);
  fOutputList->Add(hInvMassProtonXi_3);
  fOutputList->Add(hInvMassProtonXi_4);
  fOutputList->Add(hInvMassProtonXi_5);
  fOutputList->Add(hInvMassProtonXi_6);
  fOutputList->Add(hInvMassProtonXi_7);
  fOutputList->Add(hInvMassProtonXi_8);
  fOutputList->Add(hInvMassProtonXi_9);
  fOutputList->Add(hInvMassProtonXi_10);
  fOutputList->Add(hInvMassProtonXi_11);
  fOutputList->Add(hInvMassProtonXi_12);
  fOutputList->Add(hInvMassProtonXi_13);
  fOutputList->Add(hInvMassProtonXi_14);
  fOutputList->Add(hInvMassProtonXi_15);
  fOutputList->Add(hInvMassProtonXi_16);
  fOutputList->Add(hInvMassProtonXi_17);
  fOutputList->Add(hInvMassAntiProtonXi_1);
  fOutputList->Add(hInvMassAntiProtonXi_2);
  fOutputList->Add(hInvMassAntiProtonXi_3);
  fOutputList->Add(hInvMassAntiProtonXi_4);
  fOutputList->Add(hInvMassAntiProtonXi_5);
  fOutputList->Add(hInvMassAntiProtonXi_6);
  fOutputList->Add(hInvMassAntiProtonXi_7);
  fOutputList->Add(hInvMassAntiProtonXi_8);
  fOutputList->Add(hInvMassAntiProtonXi_9);
  fOutputList->Add(hInvMassAntiProtonXi_10);
  fOutputList->Add(hInvMassAntiProtonXi_11);
  fOutputList->Add(hInvMassAntiProtonXi_12);
  fOutputList->Add(hInvMassAntiProtonXi_13);
  fOutputList->Add(hInvMassAntiProtonXi_14);
  fOutputList->Add(hInvMassAntiProtonXi_15);
  fOutputList->Add(hInvMassAntiProtonXi_16);
  fOutputList->Add(hInvMassAntiProtonXi_17);
  //FG LambdaXi                                                                                                                     
  fOutputList->Add(hInvMassLambdaXi_1);
  fOutputList->Add(hInvMassLambdaXi_2);
  fOutputList->Add(hInvMassLambdaXi_3);
  fOutputList->Add(hInvMassLambdaXi_4);
  fOutputList->Add(hInvMassLambdaXi_5);
  fOutputList->Add(hInvMassLambdaXi_6);
  fOutputList->Add(hInvMassLambdaXi_7);
  fOutputList->Add(hInvMassLambdaXi_8);
  fOutputList->Add(hInvMassLambdaXi_9);
  fOutputList->Add(hInvMassLambdaXi_10);
  fOutputList->Add(hInvMassLambdaXi_11);
  fOutputList->Add(hInvMassLambdaXi_12);
  fOutputList->Add(hInvMassLambdaXi_13);
  fOutputList->Add(hInvMassLambdaXi_14);
  fOutputList->Add(hInvMassLambdaXi_15);
  fOutputList->Add(hInvMassAntiLambdaXi_1);
  fOutputList->Add(hInvMassAntiLambdaXi_2);
  fOutputList->Add(hInvMassAntiLambdaXi_3);
  fOutputList->Add(hInvMassAntiLambdaXi_4);
  fOutputList->Add(hInvMassAntiLambdaXi_5);
  fOutputList->Add(hInvMassAntiLambdaXi_6);
  fOutputList->Add(hInvMassAntiLambdaXi_7);
  fOutputList->Add(hInvMassAntiLambdaXi_8);
  fOutputList->Add(hInvMassAntiLambdaXi_9);
  fOutputList->Add(hInvMassAntiLambdaXi_10);
  fOutputList->Add(hInvMassAntiLambdaXi_11);
  fOutputList->Add(hInvMassAntiLambdaXi_12);
  fOutputList->Add(hInvMassAntiLambdaXi_13);
  fOutputList->Add(hInvMassAntiLambdaXi_14);
  fOutputList->Add(hInvMassAntiLambdaXi_15);
  //BG LambdaLambda
  fOutputList->Add(hInvMassLambdaLambda_evtpool1);
  fOutputList->Add(hInvMassLambdaLambda_evtpool2);
  fOutputList->Add(hInvMassLambdaLambda_evtpool3);
  fOutputList->Add(hInvMassLambdaLambda_evtpool4);
  fOutputList->Add(hInvMassLambdaLambda_evtpool5);
  fOutputList->Add(hInvMassLambdaLambda_evtpool6);
  fOutputList->Add(hInvMassLambdaLambda_evtpool7);
  fOutputList->Add(hInvMassLambdaLambda_evtpool8);
  fOutputList->Add(hInvMassLambdaLambda_evtpool9);
  fOutputList->Add(hInvMassLambdaLambda_evtpool10);
  fOutputList->Add(hInvMassLambdaLambda_evtpool11);
  fOutputList->Add(hInvMassLambdaLambda_evtpool12);
  fOutputList->Add(hInvMassLambdaLambda_evtpool13);
  fOutputList->Add(hInvMassLambdaLambda_evtpool14);
  fOutputList->Add(hInvMassLambdaLambda_evtpool15);
  fOutputList->Add(hInvMassLambdaLambda_evtpool16);
  fOutputList->Add(hInvMassLambdaLambda_evtpool17);
  fOutputList->Add(hInvMassAntiLambdaAntiLambda_evtpool1);
  fOutputList->Add(hInvMassAntiLambdaAntiLambda_evtpool2);
  fOutputList->Add(hInvMassAntiLambdaAntiLambda_evtpool3);
  fOutputList->Add(hInvMassAntiLambdaAntiLambda_evtpool4);
  fOutputList->Add(hInvMassAntiLambdaAntiLambda_evtpool5);
  fOutputList->Add(hInvMassAntiLambdaAntiLambda_evtpool6);
  fOutputList->Add(hInvMassAntiLambdaAntiLambda_evtpool7);
  fOutputList->Add(hInvMassAntiLambdaAntiLambda_evtpool8);
  fOutputList->Add(hInvMassAntiLambdaAntiLambda_evtpool9);
  fOutputList->Add(hInvMassAntiLambdaAntiLambda_evtpool10);
  fOutputList->Add(hInvMassAntiLambdaAntiLambda_evtpool11);
  fOutputList->Add(hInvMassAntiLambdaAntiLambda_evtpool12);
  fOutputList->Add(hInvMassAntiLambdaAntiLambda_evtpool13);
  fOutputList->Add(hInvMassAntiLambdaAntiLambda_evtpool14);
  fOutputList->Add(hInvMassAntiLambdaAntiLambda_evtpool15);
  fOutputList->Add(hInvMassAntiLambdaAntiLambda_evtpool16);
  fOutputList->Add(hInvMassAntiLambdaAntiLambda_evtpool17);
  //BG ProtonXi
  fOutputList->Add(hInvMassProtonXi_evtpool1);
  fOutputList->Add(hInvMassProtonXi_evtpool2);
  fOutputList->Add(hInvMassProtonXi_evtpool3);
  fOutputList->Add(hInvMassProtonXi_evtpool4);
  fOutputList->Add(hInvMassProtonXi_evtpool5);
  fOutputList->Add(hInvMassProtonXi_evtpool6);
  fOutputList->Add(hInvMassProtonXi_evtpool7);
  fOutputList->Add(hInvMassProtonXi_evtpool8);
  fOutputList->Add(hInvMassProtonXi_evtpool9);
  fOutputList->Add(hInvMassProtonXi_evtpool10);
  fOutputList->Add(hInvMassProtonXi_evtpool11);
  fOutputList->Add(hInvMassProtonXi_evtpool12);
  fOutputList->Add(hInvMassProtonXi_evtpool13);
  fOutputList->Add(hInvMassProtonXi_evtpool14);
  fOutputList->Add(hInvMassProtonXi_evtpool15);
  fOutputList->Add(hInvMassProtonXi_evtpool16);
  fOutputList->Add(hInvMassProtonXi_evtpool17);
  fOutputList->Add(hInvMassAntiProtonXi_evtpool1);
  fOutputList->Add(hInvMassAntiProtonXi_evtpool2);
  fOutputList->Add(hInvMassAntiProtonXi_evtpool3);
  fOutputList->Add(hInvMassAntiProtonXi_evtpool4);
  fOutputList->Add(hInvMassAntiProtonXi_evtpool5);
  fOutputList->Add(hInvMassAntiProtonXi_evtpool6);
  fOutputList->Add(hInvMassAntiProtonXi_evtpool7);
  fOutputList->Add(hInvMassAntiProtonXi_evtpool8);
  fOutputList->Add(hInvMassAntiProtonXi_evtpool9);
  fOutputList->Add(hInvMassAntiProtonXi_evtpool10);
  fOutputList->Add(hInvMassAntiProtonXi_evtpool11);
  fOutputList->Add(hInvMassAntiProtonXi_evtpool12);
  fOutputList->Add(hInvMassAntiProtonXi_evtpool13);
  fOutputList->Add(hInvMassAntiProtonXi_evtpool14);
  fOutputList->Add(hInvMassAntiProtonXi_evtpool15);
  fOutputList->Add(hInvMassAntiProtonXi_evtpool16);
  fOutputList->Add(hInvMassAntiProtonXi_evtpool17);
  //BG LambdaXi
  fOutputList->Add(hInvMassLambdaXi_evtpool1);
  fOutputList->Add(hInvMassLambdaXi_evtpool2);
  fOutputList->Add(hInvMassLambdaXi_evtpool3);
  fOutputList->Add(hInvMassLambdaXi_evtpool4);
  fOutputList->Add(hInvMassLambdaXi_evtpool5);
  fOutputList->Add(hInvMassLambdaXi_evtpool6);
  fOutputList->Add(hInvMassLambdaXi_evtpool7);
  fOutputList->Add(hInvMassLambdaXi_evtpool8);
  fOutputList->Add(hInvMassLambdaXi_evtpool9);
  fOutputList->Add(hInvMassLambdaXi_evtpool10);
  fOutputList->Add(hInvMassLambdaXi_evtpool11);
  fOutputList->Add(hInvMassLambdaXi_evtpool12);
  fOutputList->Add(hInvMassLambdaXi_evtpool13);
  fOutputList->Add(hInvMassLambdaXi_evtpool14);
  fOutputList->Add(hInvMassLambdaXi_evtpool15);
  fOutputList->Add(hInvMassAntiLambdaXi_evtpool1);
  fOutputList->Add(hInvMassAntiLambdaXi_evtpool2);
  fOutputList->Add(hInvMassAntiLambdaXi_evtpool3);
  fOutputList->Add(hInvMassAntiLambdaXi_evtpool4);
  fOutputList->Add(hInvMassAntiLambdaXi_evtpool5);
  fOutputList->Add(hInvMassAntiLambdaXi_evtpool6);
  fOutputList->Add(hInvMassAntiLambdaXi_evtpool7);
  fOutputList->Add(hInvMassAntiLambdaXi_evtpool8);
  fOutputList->Add(hInvMassAntiLambdaXi_evtpool9);
  fOutputList->Add(hInvMassAntiLambdaXi_evtpool10);
  fOutputList->Add(hInvMassAntiLambdaXi_evtpool11);
  fOutputList->Add(hInvMassAntiLambdaXi_evtpool12);
  fOutputList->Add(hInvMassAntiLambdaXi_evtpool13);
  fOutputList->Add(hInvMassAntiLambdaXi_evtpool14);
  fOutputList->Add(hInvMassAntiLambdaXi_evtpool15);

  //===============

  PostData(1,fOutputList);
}

//________________________________________________________________________
void AliAnalysisTaskPPiLambda::UserExec(Option_t *)
{
  //cout<<"++++++++++++++++++++"<<endl;
  //cout<<"++ Analysis Start ++"<<endl;
  //cout<<"++++++++++++++++++++"<<endl;

  // Event loop
  //if(nevt!=0) return;
  
  // Input event
  AliAnalysisManager   *man         =AliAnalysisManager::GetAnalysisManager();
  AliInputEventHandler *inputHandler=(AliInputEventHandler*)(man->GetInputEventHandler());
  
  if(!inputHandler){
    printf("ERROR: AliInputEventHandler not available\n");
    return;
  }
  fAOD=dynamic_cast<AliAODEvent*>(InputEvent());
  if(!fAOD){
    printf("ERROR: fAOD not available\n");
    return;
  }
  fHeader=(AliAODHeader*)fAOD->GetHeader();
  if(!fHeader){
    printf("ERROR: fHeader not available\n");
    return;
  }
  fUtils=new AliAnalysisUtils();
  if(!fUtils){
    printf("ERROR: fUtils not available\n");
    return;
  }
  if(inputHandler){
    fPIDResponse=inputHandler->GetPIDResponse();
    if(!fPIDResponse){
      printf("ERROR: fPIDResponse not available\n");
      return;
    }
  }

  //=============== Event selection
  Bool_t isEvtMB=((((AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler()))->IsEventSelected()) & ftrigBit);
  if(!isEvtMB) return;
  fEvtCounter->Fill(0.5);
  // Event pile-up
  Bool_t isPileupEvt=kFALSE;
  isPileupEvt=(fAOD->IsPileupFromSPD() || fUtils->IsSPDClusterVsTrackletBG(fAOD));
  Bool_t isGoodEvt=EventSelection(fAOD);
  if(!isGoodEvt || isPileupEvt){ PostData(1,fOutputList); return; }
  fEvtPassCut->Fill(0.5);
  //if(nevt%1000==0) cout<<nevt<<"-th event"<<endl;
  nevt++; 
  
  //============= variable definition
  Int_t    nTracks   =fAOD->GetNumberOfTracks();
  Int_t    nV0s      =fAOD->GetNumberOfV0s();
  Int_t    nCascades =fAOD->GetNumberOfCascades();
  Double_t lamntrack[6][100]    ={{0},{0},{0},{0},{0},{0}};
  Double_t lambtrack[6][100]    ={{0},{0},{0},{0},{0},{0}};
  Double_t Lambda[10000]={0};
  Double_t AntiLambda[10000]={0};

  Float_t  vecTarget[3]={0.};
  Int_t multiplicity =0;
  
  vecTarget[0]=fAOD->GetPrimaryVertex()->GetX(); // primary vertex
  vecTarget[1]=fAOD->GetPrimaryVertex()->GetY();
  vecTarget[2]=fAOD->GetPrimaryVertex()->GetZ();
  
  AliAODTracklets *tracklets=(AliAODTracklets*)fAOD->GetTracklets();
  Int_t nTracklet = tracklets->GetNumberOfTracklets();
  
  for(Int_t itracklet = 0; itracklet<nTracklet; itracklet++){
    Double_t tracklet_eta   = tracklets->GetEta(itracklet);    
    if(fabs(tracklet_eta) > 0.8) continue;
    multiplicity++;    
  }
  
  //============== cascade loop
  TObjArray* fCascadeArraym=new TObjArray();
  TObjArray* fCascadeArrayp=new TObjArray();

  Int_t    nXim=0; // Xi-/event
  Int_t    nXip=0; // Xi+/event

  for (int iCasc=0; iCasc< nCascades; iCasc++){
    AliAODcascade *casc= (AliAODcascade*)fAOD->GetCascade(iCasc);
    if(!casc) continue;

    AliAODTrack* pTrackXi   = dynamic_cast<AliAODTrack*>(casc->GetDaughter(0));                     // daughter track +
    AliAODTrack* nTrackXi   = dynamic_cast<AliAODTrack*>(casc->GetDaughter(1));                     // daughter track -
    AliAODTrack* bachTrackXi= dynamic_cast<AliAODTrack*>(casc->GetDecayVertexXi()->GetDaughter(0)); // bachlor track
    double_t dcaxyp=0.,dcazp=0.,dcap=0.;  // daughter positive track
    double_t dcaxyn=0.,dcazn=0.,dcan=0.;  // daughter negative track
    double_t dcaxyb=0.,dcazb=0.,dcab=0.;  // bachelor track
    Float_t  dDCAp[2]    ={0.};        
    Float_t  cDCAp[3]    ={0.};        
    Float_t  dDCAn[2]    ={0.}; 
    Float_t  cDCAn[3]    ={0.}; 
    Float_t  dDCAb[2]    ={0.}; 
    Float_t  cDCAb[3]    ={0.}; 
    double_t v0Vtx[3]    ={0.};
    double_t xivertex[3] ={0.};
    Double_t lambdacpa=0.,lambdatransradius=0.;
    Double_t xicpa=0.,xitransradius=0.;
    Double_t v0vtxx=0.,v0vtxy=0.,v0vtxz=0.,xivtxx=0.,xivtxy=0.,xivtxz=0.;

    if((pTrackXi->Charge()) < (nTrackXi->Charge())) continue;
    
    //========== Daughter track selection 
    if(fabs(pTrackXi   ->Eta()) > 0.8)       continue;
    if(fabs(nTrackXi   ->Eta()) > 0.8)       continue;
    if(fabs(bachTrackXi->Eta()) > 0.8)       continue;
    if(pTrackXi   ->GetTPCNcls() < 70) continue;
    if(nTrackXi   ->GetTPCNcls() < 70) continue;
    if(bachTrackXi->GetTPCNcls() < 70) continue;
    if(pTrackXi   ->Pt() < 0.3)        continue;
    if(nTrackXi   ->Pt() < 0.3)        continue;
    if(bachTrackXi->Pt() < 0.3)        continue;

    pTrackXi->GetImpactParameters(dDCAp,cDCAp);
    dcaxyp          =dDCAp[0];
    dcazp           =dDCAp[1];
    dcap            =sqrt(dcaxyp*dcaxyp+dcazp*dcazp);
    nTrackXi->GetImpactParameters(dDCAn,cDCAn);
    dcaxyn          =dDCAn[0];
    dcazn           =dDCAn[1];
    dcan            =sqrt(dcaxyn*dcaxyn+dcazn*dcazn);
    bachTrackXi->GetImpactParameters(dDCAb,cDCAb);
    dcaxyb          =dDCAb[0];
    dcazb           =dDCAb[1];
    dcab            =sqrt(dcaxyb*dcaxyb+dcazb*dcazb);

    if(dcap < 0.05) continue;
    if(dcan < 0.05) continue;
    if(dcab < 0.05) continue;

    //========== PID
    Bool_t isProton     =(fabs(fPIDResponse->NumberOfSigmasTPC(pTrackXi, AliPID::kProton)) < 4);   // proton
    Bool_t isAntiproton =(fabs(fPIDResponse->NumberOfSigmasTPC(nTrackXi, AliPID::kProton)) < 4);   // proton- 
    Bool_t isPospion    =(fabs(fPIDResponse->NumberOfSigmasTPC(pTrackXi, AliPID::kPion)) < 4);     // pion+
    Bool_t isNegpion    =(fabs(fPIDResponse->NumberOfSigmasTPC(nTrackXi, AliPID::kPion)) < 4);     // pion-
    Bool_t Bachpion     =(fabs(fPIDResponse->NumberOfSigmasTPC(bachTrackXi, AliPID::kPion)) < 4);  // bacelor

    Bool_t isXi         =((isProton) && (isNegpion) && (bachTrackXi->Charge() < 0));
    Bool_t isXibar      =((isAntiproton) && (isPospion) && (bachTrackXi->Charge() > 0));

    //========== lambda selection
    v0vtxx               =casc->DecayVertexV0X();                      
    v0vtxy               =casc->DecayVertexV0Y(); 
    v0vtxz               =casc->DecayVertexV0Z();
    v0Vtx[0]             =v0vtxx;
    v0Vtx[1]             =v0vtxy;
    v0Vtx[2]             =v0vtxz;
    lambdacpa            =LambdaCosPointingAngle(casc,v0Vtx,vecTarget);
    lambdatransradius    =DecayLengthXY(v0Vtx,vecTarget);
        
    if(lambdacpa < 0.97)                 continue;
    if(lambdatransradius < 1.4)          continue;
    if(lambdatransradius > 200)          continue;
    if(casc->DcaV0Daughters() > 1.5)     continue;
    if(casc->DcaV0ToPrimVertex() < 0.07) continue;
    if(casc->GetOnFlyStatus())           continue; // select offline

    if(isXi){
      hInvMassXidecayLambda     ->Fill(InvMassLambda(casc));               
    }
    if(isXibar){
      hInvMassXidecayantiLambda ->Fill(InvMassAntiLambda(casc));      
    }

   //===== xi selection ===== 
    xivtxx        =casc->DecayVertexXiX();                                     
    xivtxy        =casc->DecayVertexXiY();
    xivtxz        =casc->DecayVertexXiZ();
    xivertex[0]   =xivtxx;
    xivertex[1]   =xivtxy;
    xivertex[2]   =xivtxz;
    xicpa         =casc->CosPointingAngleXi(vecTarget[0],vecTarget[1],vecTarget[2]);
    xitransradius =xiDecayLengthXY(xivertex,vecTarget);

    if(xicpa < 0.98)                 continue;
    if(xitransradius < 0.8)          continue;
    if(xitransradius > 200)          continue;
    if(casc->DcaXiDaughters() > 1.6) continue;
    if((InvMassOmega(casc) > 1.667) && (InvMassOmega(casc) < 1.677)) continue;

    //=========== lambda + pion- -> xi-
    if(isXi){

      if(InvMassLambda(casc) < 1.109) continue;
      if(InvMassLambda(casc) > 1.121) continue;

      hInvMassXi            ->Fill(InvMassXi(casc));

      if(InvMassXi(casc) < 1.317)     continue;
      if(InvMassXi(casc) > 1.327)     continue;

      XiCPA                 ->Fill(xicpa);        
      XiTranverseradius     ->Fill(xitransradius);
      dcaXiDaughters        ->Fill(casc->DcaXiDaughters());    
      fPtcascade            ->Fill(sqrt(casc->Pt2Xi()));
      fPtv0casc             ->Fill(sqrt(casc->Pt2V0()));
      fPtBachelor           ->Fill(sqrt(pow(casc->MomBachX(),2)+pow(casc->MomBachY(),2)));
      fPhiXidecayLambda     ->Fill(casc->OpenAngleV0());

      nXim++;

      fCascadeArraym->Add(casc);
    }
    //========== antilambda + pion+ -> xi+
    if(isXibar){
     
      if(InvMassAntiLambda(casc) < 1.109) continue;
      if(InvMassAntiLambda(casc) > 1.121) continue;

      hInvMassposXi             ->Fill(InvMassXi(casc));
      
      if(InvMassXi(casc) < 1.317)         continue;
      if(InvMassXi(casc) > 1.327)         continue;

      fPtcascade_xip        ->Fill(sqrt(casc->Pt2Xi()));
      fPtv0casc_xip         ->Fill(sqrt(casc->Pt2V0()));
      fPtBachelor_xip       ->Fill(sqrt(pow(casc->MomBachX(),2)+pow(casc->MomBachY(),2)));
      fPhiXidecayLambda_xip ->Fill(casc->OpenAngleV0());

      nXip++;

      fCascadeArrayp->Add(casc);
    }    
  }
  //============== cascade loop end

  //========== v0 loop
  TObjArray* fV0Arrayn=new TObjArray();
  TObjArray* fV0Arrayb=new TObjArray();

  int nLamn=0;
  int nLamb=0;

  for (Int_t iV0=0; iV0<nV0s; iV0++){
    AliAODv0 *v0= (AliAODv0*)fAOD->GetV0(iV0);
    if(!v0) continue; 

    AliAODTrack* ptrack=dynamic_cast<AliAODTrack*>(v0->GetDaughter(0));    
    AliAODTrack* ntrack=dynamic_cast<AliAODTrack*>(v0->GetDaughter(1));    
    Double_t dcaxyp=0.,dcazp=0.,dcap=0.; // proton
    Double_t dcaxyn=0.,dcazn=0.,dcan=0.; // pion
    Float_t  dDCAp[2]={0.};              // DCA to the vertex xy(d) and z
    Float_t  dDCAn[2]={0.};  
    Float_t  cDCAp[3]={0.};              // convariance of impact parameters
    Float_t  cDCAn[3]={0.}; 
    Double_t v0Vtx[3]={0.};
    Double_t energy=0.; 
    Double_t fcpa=0.,ftransradius=0.;

    if((ptrack->Charge()) < (ntrack->Charge())) continue;  

    ptrack->GetImpactParameters(dDCAp,cDCAp);
    dcaxyp          =dDCAp[0];
    dcazp           =dDCAp[1];
    dcap            =sqrt(dcaxyp*dcaxyp+dcazp*dcazp);
    ntrack->GetImpactParameters(dDCAn,cDCAn);
    dcaxyn          =dDCAn[0];
    dcazn           =dDCAn[1];
    dcan            =sqrt(dcaxyn*dcaxyn+dcazn*dcazn);

    //========== v0 daughter tracks selection
    if(fabs(ptrack->Eta()) > 0.8) continue;
    if(fabs(ntrack->Eta()) > 0.8) continue;
    if(ptrack->GetTPCNcls() < 70) continue;
    if(ntrack->GetTPCNcls() < 70) continue;
    if(dcap < 0.05) continue;
    if(dcan < 0.05) continue;

    //========== Out-of-bunch pile-up removal
    if(!(ptrack->HasPointOnITSLayer(0)) && !(ptrack->HasPointOnITSLayer(1))&&
       !(ptrack->HasPointOnITSLayer(4)) && !(ptrack->HasPointOnITSLayer(5))&&
       !(ptrack->GetTOFBunchCrossing()==0)) continue;
    if(!(ntrack->HasPointOnITSLayer(0)) && !(ntrack->HasPointOnITSLayer(1))&&
       !(ntrack->HasPointOnITSLayer(4)) && !(ntrack->HasPointOnITSLayer(5))&&
       !(ntrack->GetTOFBunchCrossing()==0)) continue;

    //========== daughter track PID
    Bool_t Daugproton     =(fabs(fPIDResponse->NumberOfSigmasTPC(ptrack, AliPID::kProton)) < 5); // proton +
    Bool_t Daugantiproton =(fabs(fPIDResponse->NumberOfSigmasTPC(ntrack, AliPID::kProton)) < 5); // proton -
    Bool_t Daugpion       =(fabs(fPIDResponse->NumberOfSigmasTPC(ntrack, AliPID::kPion)) < 5);   // pion -
    Bool_t Daugnpion      =(fabs(fPIDResponse->NumberOfSigmasTPC(ptrack, AliPID::kPion)) < 5);   // pion +
    Bool_t lambda         =(Daugproton && Daugpion);
    Bool_t antilambda     =(Daugnpion && Daugantiproton);

    v0Vtx[0]     =v0->DecayVertexV0X();
    v0Vtx[1]     =v0->DecayVertexV0Y();
    v0Vtx[2]     =v0->DecayVertexV0Z();
    fcpa         =CosPointingAngle(v0,v0Vtx,vecTarget);
    ftransradius =DecayLengthXY(v0Vtx,vecTarget);

    
    if(lambda){
      fPtv0->Fill(InvMasslambda(v0),sqrt(v0->Pt2V0()));
    }
    if(antilambda){
      fPtv0->Fill(InvMassAntilambda(v0),sqrt(v0->Pt2V0()));
    }

    //========== v0 selection
    if(sqrt(v0->Pt2V0()) < 0.3) continue;
    if(fabs(v0Vtx[0]) > 100) continue;
    if(fabs(v0Vtx[1]) > 100) continue;
    if(fabs(v0Vtx[2]) > 100) continue;
    if(ftransradius < 0.2) continue;
    if(ftransradius > 100) continue;
    if(v0->DcaV0Daughters() > 1.5) continue;
    if(fcpa < 0.99) continue;
    if((InvMassK0(v0) > 0.48) && (InvMassK0(v0) < 0.515)) continue;

    if(v0->GetOnFlyStatus()) continue; // select offline v0

    //========== lambda invariant mass
    if(lambda){
      
      fTransRadius    ->Fill(ftransradius);
      fDCAdaugTov0Vtx ->Fill(v0->DcaV0Daughters());
      fCPA            ->Fill(fcpa);
      hInvMassLambda  ->Fill(InvMasslambda(v0));

      if(InvMasslambda(v0) < 1.111) continue; 
      if(InvMasslambda(v0) > 1.119) continue; // lambda selection ±4 MeV/c^2

      fPtv0_lambda   ->Fill(InvMasslambda(v0),sqrt(v0->Pt2V0()));

  /*
      lamntrack[0][nLamn]=sqrt(pow(1.115,2)+pow(v0->MomV0X(),2)+pow(v0->MomV0Y(),2)+pow(v0->MomV0Z(),2));
      lamntrack[1][nLamn]=v0->MomV0X();
      lamntrack[2][nLamn]=v0->MomV0Y();
      lamntrack[3][nLamn]=v0->MomV0Z();
      lamntrack[4][nLamn]=v0->GetPosID();
      lamntrack[5][nLamn]=v0->GetNegID();
  */
      nLamn++;
      lambdacounter++;

      fV0Arrayn->Add(v0);
    }

    //========== antilambda invariant mass
    if(antilambda){    

      fPtv0_antilambda   ->Fill(sqrt(v0->Pt2V0()));
      hInvMassAntiLambda ->Fill(InvMassAntilambda(v0));
      
      if(InvMassAntilambda(v0) < 1.111) continue; 
      if(InvMassAntilambda(v0) > 1.119) continue; // lambda selection ±4 MeV/c^2
 
      fPtv0_antilambda_mass_masscut ->Fill(InvMassAntilambda(v0),sqrt(v0->Pt2V0()));

  /*
      lambtrack[0][nLamb]=sqrt(pow(1.115,2)+pow(v0->MomV0X(),2)+pow(v0->MomV0Y(),2)+pow(v0->MomV0Z(),2));
      lambtrack[1][nLamb]=v0->MomV0X();
      lambtrack[2][nLamb]=v0->MomV0Y();
      lambtrack[3][nLamb]=v0->MomV0Z();
      lambtrack[4][nLamb]=v0->GetPosID();
      lambtrack[5][nLamb]=v0->GetNegID();
  */ 
     nLamb++; 

      fV0Arrayb->Add(v0);
    }    
  }
  //=============== V0 loop end
  
  //=============== Track loop
  TObjArray* fProtonArray=new TObjArray();
  TObjArray* fProtonArrayb=new TObjArray();

  int np =0;
  int npb=0;

  Int_t pid=0;
  Int_t antipid=0;

  for(Int_t i=0; i<nTracks; i++){
    AliAODTrack *track=static_cast<AliAODTrack*>(fAOD->GetTrack(i));
    if(!track) continue;

    Double_t tpcratio=0.,dcaxy=0.,dcaz=0.,tofsig=0.,trklength=0.,beta=0.;
    UShort_t crsR=0,clsF=0; 
    Float_t  sigmatpc_pr=0.,sigmatof=0.,sigmacomb=0.;
    Float_t  dDCA[2]    ={0.};        
    Float_t  cDCA[3]    ={0.};        

    //========== proton selection 
    crsR         =track->GetTPCNCrossedRows();
    clsF         =track->GetTPCNclsF();
    tpcratio     =(clsF>0)?(double)crsR/clsF:0;
    track->GetImpactParameters(dDCA,cDCA); 
    dcaxy        =dDCA[0];
    dcaz         =dDCA[1];
    sigmatpc_pr  =fPIDResponse->NumberOfSigmasTPC(track, AliPID::kProton);
    sigmatof     =fPIDResponse->NumberOfSigmasTOF(track, AliPID::kProton);
    sigmacomb    =sqrt(pow(sigmatpc_pr,2)+pow(sigmatof,2));    
    tofsig       =track->GetTOFsignal();
    trklength    =track->GetIntegratedLength();
    beta         =(tofsig>0)?(double)trklength/ (2.99792457999999984e-02 * tofsig):0;

    if(fabs(track->Eta()) > 0.8) continue;
    if(track->Pt() < 0.5) continue;
    if(track->Pt() > 4.05) continue;
    if(track->GetTPCNcls() < 80) continue;
    if(crsR < 70) continue;
    if(tpcratio < 0.83) continue;
    if(!(track->GetTPCSharedMap()==0)) continue;
    if(fabs(dcaxy) > 0.1) continue;
    if(fabs(dcaz) > 0.2) continue;
    Bool_t isProton  =((track->P() < 0.75 && fabs(sigmatpc_pr) < 3) || (track->P() > 0.75 && sigmacomb < 3));
    //Bool_t isProton  =((track->P() < 0.75 && fabs(sigmatpc_pr) < 3) || (track->P() > 0.75 && ((fabs(sigmatpc_pr) < 3)||(fabs(sigmatof) < 3)) ));

    fTrkEta         ->Fill(track->Eta());   
    fTrkPt          ->Fill(track->Pt());
    fTPCclus        ->Fill(track->GetTPCNcls());
    fTPCcrsR        ->Fill(crsR);   
    fTPCclusF       ->Fill(clsF);
    fTrkDCAxy       ->Fill(dcaxy);
    fTrkDCAz        ->Fill(dcaz);
    fPvsTPCsignal   ->Fill(track->P(),track->GetTPCsignal());
    fPvsTOFsignal   ->Fill(track->P(),beta);  

    if(!isProton) continue;

    //========== proton
    if(track->Charge() > 0){
      fPt_proton               ->Fill(track->Pt());
      PID_cut                  ->Fill(track->P(),track->GetTPCsignal());  
 
      np++;

      fProtonArray->Add(track);
    }

    //========== anti proton
    if(track->Charge() < 0){    

      fPt_antiproton            ->Fill(track->Pt());

      npb++;
      
      fProtonArrayb->Add(track);
    }

  }
  //=================== track loop end

  // remove same track id 
  for(Int_t i=0; i < fProtonArray->GetEntriesFast(); i++){
    TObject *objProton1=fProtonArray->At(i);    
    if(!objProton1) continue;
    AliAODTrack *proton1=dynamic_cast<AliAODTrack*>(objProton1); 
    Int_t id1 =proton1->GetID();
    if(id1 < 0) id1 =-id1-1;

    for(Int_t j=i+1; j < fProtonArray->GetEntriesFast();){
      TObject *objProton2=fProtonArray->At(j);    
      if(!objProton2) { j++; continue;}
      AliAODTrack *proton2=dynamic_cast<AliAODTrack*>(objProton2); 
      Int_t id2 =proton2->GetID();
      if(id2 < 0) id2 =-id2-1;

      if(id1==id2) fProtonArray->RemoveAt(j);
      else j++;
    }
  }
 
  for(Int_t i=0; i < fProtonArrayb->GetEntriesFast(); i++){
    TObject *objProton1=fProtonArrayb->At(i);    
    if(!objProton1) continue;
    AliAODTrack *proton1=dynamic_cast<AliAODTrack*>(objProton1); 
    Int_t id1 =proton1->GetID();
    if(id1 < 0) id1 =-id1-1;

    for(Int_t j=i+1; j < fProtonArrayb->GetEntriesFast();){
      TObject *objProton2=fProtonArrayb->At(j);    
      if(!objProton2) { j++; continue;}
      AliAODTrack *proton2=dynamic_cast<AliAODTrack*>(objProton2); 
      Int_t id2 =proton2->GetID();
      if(id2 < 0) id2 =-id2-1;

      if(id1==id2) fProtonArrayb->RemoveAt(j);
      else j++;
    }
  }

  // decay rejection
  // Lambda
  if(nLamn>0){
    for(Int_t i=0; i < fV0Arrayn->GetEntriesFast(); i++){
      TObject *objV0=fV0Arrayn->At(i);    
      if(!objV0) continue;
      AliAODv0 *lambda=dynamic_cast<AliAODv0*>(objV0); 

      if(nXim > 0){
	for(Int_t j=0; j < fCascadeArraym->GetEntriesFast(); j++){
	  TObject *objCascade=fCascadeArraym->At(j);    
	  if(!objCascade) continue;
	  AliAODcascade *xi=dynamic_cast<AliAODcascade*>(objCascade); 
	  Bool_t isXidecayLambda = ((lambda->GetPosID()) == xi->GetPosID()) && ((lambda->GetNegID()) == (xi->GetNegID()));
	  // xi decay lambda → 1 + 0 + 0 = 1, prompt lambda → 0 + 0 + 0 = 0
	  Lambda[i] += isXidecayLambda;
	}
      }
      if(nXim == 0){
        Lambda[i] = 0;
      }
    }
  }

  // Lambda pT
  if(nLamn > 0){
    for(Int_t i=0; i < fV0Arrayn->GetEntriesFast(); i++){
      TObject *objV0=fV0Arrayn->At(i);    
      if(!objV0) continue;
      AliAODv0 *lambda=dynamic_cast<AliAODv0*>(objV0); 

      fPt_allLambda->Fill(sqrt(pow(lambda->MomV0X(),2)+pow(lambda->MomV0Y(),2)));
	
	if(Lambda[i]==1){
	  fPt_xidecaylambda->Fill(sqrt(pow(lambda->MomV0X(),2)+pow(lambda->MomV0Y(),2)));	  
	}	
    }
  }

  //antiLambda
  if(nLamb>0){
    for(Int_t i=0; i < fV0Arrayb->GetEntriesFast(); i++){
      TObject *objV0=fV0Arrayb->At(i);    
      if(!objV0) continue;
      AliAODv0 *lambda=dynamic_cast<AliAODv0*>(objV0); 
      if(nXip > 0){
	for(Int_t j=0; j < fCascadeArrayp->GetEntriesFast(); j++){
	  TObject *objCascade=fCascadeArrayp->At(j);    
	  if(!objCascade) continue;
	  AliAODcascade *xi=dynamic_cast<AliAODcascade*>(objCascade); 
	  
	  Bool_t isXidecayAntiLambda = ((lambda->GetPosID()) == xi->GetPosID()) && ((lambda->GetNegID()) == (xi->GetNegID()));
	  AntiLambda[i] += isXidecayAntiLambda;
	}
      }
      if(nXip == 0){
        AntiLambda[i] = 0;
      }

    }
  }

  // antiLambda pT
  if(nLamb > 0){
    for(Int_t i=0; i < fV0Arrayb->GetEntriesFast(); i++){
      TObject *objV0=fV0Arrayb->At(i);    
      if(!objV0) continue;
      AliAODv0 *lambda=dynamic_cast<AliAODv0*>(objV0); 

      fPt_allAntiLambda->Fill(sqrt(pow(lambda->MomV0X(),2)+pow(lambda->MomV0Y(),2)));

      if(AntiLambda[i]==1){
	fPt_xidecayAntiLambda->Fill(sqrt(pow(lambda->MomV0X(),2)+pow(lambda->MomV0Y(),2)));	  
      }	
    }
  }

  //==================
  Double_t energysum=0.,energy1=0.,energy2=0.,pt=0.,l1px=0.,l1py=0.,l1pz=0.,l2px=0.,l2py=0.,l2pz=0.,
    invMass=0.,relmom=0.,openangle=0.;
  Double_t isXi1=0.,isXi2=0.;
  Short_t posID1=0,posID2=0,negID1=0,negID2=0;
  Short_t posID=0,negID=0;

  Double_t energyxi=0.,energypr=0.,energylam=0.,prpx=0.,prpy=0.,prpz=0.,xipx=0.,xipy=0.,xipz=0.,lampx=0.,lampy=0.,lampz=0.;
  Short_t bachid=0,ptrackid=0,ntrackid=0,protonid=0;

  Double_t protonpt1=0.,protonpt2=0.,xidecaypt=0.,protonpt=0.;

  Bool_t islambdalambda=false;
  Bool_t isantilambdaantilambda=false;
  Bool_t isprotonxi=false;
  Bool_t isantiprotonxi=false;
  Bool_t islambdaxi=false; 
  Bool_t isantilambdaxi=false; 

  // p + Xi-
  if((np > 0) && (nXim > 0)){
    for(Int_t i=0; i < fCascadeArraym->GetEntriesFast(); i++){
      TObject *objCascade=fCascadeArraym->At(i);    
      if(!objCascade) continue;
      AliAODcascade *xi=dynamic_cast<AliAODcascade*>(objCascade); 
      AliAODTrack* ptrkxi   = dynamic_cast<AliAODTrack*>(xi->GetDaughter(0)); 

      for(Int_t j=0; j < fProtonArray->GetEntriesFast(); j++){
	TObject *objProton=fProtonArray->At(j);    
	if(!objProton) continue;
	AliAODTrack *proton=dynamic_cast<AliAODTrack*>(objProton); 

	xidecaypt   =ptrkxi->Pt();
	protonpt    =proton->Pt();

	xipx        =xi->MomXiX();
	xipy        =xi->MomXiY();
	xipz        =xi->MomXiZ();
	energyxi    =sqrt(pow(1.322,2)+pow(xipx,2)+pow(xipy,2)+pow(xipz,2));
	bachid      =xi->GetBachID();
	ptrackid    =xi->GetPosID();
	ntrackid    =xi->GetNegID();
	prpx        =proton->Px();
	prpy        =proton->Py();
	prpz        =proton->Pz();
	energypr    =sqrt(pow(0.938,2)+pow(prpx,2)+pow(prpy,2)+pow(prpz,2));
	protonid    =proton->GetID();
	if(protonid < 0) protonid =-protonid-1;
	energysum   =energyxi+energypr;

	relmom    =sqrt(pow((xipx-prpx),2)+pow((xipy-prpy),2)+pow((xipz-prpz),2))/2;
	openangle =OpenAngle(xipx,xipy,xipz,prpx,prpy,prpz)*57.2957796;
	pt        =sqrt((prpx+xipx)*(prpx+xipx)+(prpy+xipy)*(prpy+xipy));
	invMass   =InvariantMass(xipx,xipy,xipz,prpx,prpy,prpz,energysum);

	if(protonid == ptrackid) continue;
	if(protonid == ntrackid) continue;
	if(protonid == bachid) continue;
	if(bachid == ptrackid) continue;
	if(bachid == ntrackid) continue;
	if(ptrackid == ntrackid) continue;

	hInvMassProtonXi            ->Fill(pt,invMass); 
	hPtcorrelation_pxi          ->Fill(xidecaypt,protonpt);
	hInvMassProtonXi_relmomentum->Fill(relmom,invMass);
	hInvMassProtonXi_openangle  ->Fill(openangle,invMass);
	fmultiplicity_protonxi      ->Fill(multiplicity);

	isprotonxi=true; 

	//multiplicity
	if((multiplicity>=1) && (multiplicity<=4)){
          hInvMassProtonXi_1 ->Fill(invMass);
        }
        if((multiplicity>=5) && (multiplicity<=8)){
          hInvMassProtonXi_2 ->Fill(invMass);
        }
        if((multiplicity>=9) && (multiplicity<=12)){
          hInvMassProtonXi_3 ->Fill(invMass);
        }
        if((multiplicity>=13) && (multiplicity<=16)){
          hInvMassProtonXi_4 ->Fill(invMass);
        }
        if((multiplicity>=17) && (multiplicity<=20)){
          hInvMassProtonXi_5 ->Fill(invMass);
        }
        if((multiplicity>=21) && (multiplicity<=24)){
          hInvMassProtonXi_6 ->Fill(invMass);
        }
        if((multiplicity>=25) && (multiplicity<=28)){
          hInvMassProtonXi_7 ->Fill(invMass);
        }
        if((multiplicity>=29) && (multiplicity<=32)){
          hInvMassProtonXi_8 ->Fill(invMass);
        }
        if((multiplicity>=33) && (multiplicity<=36)){
          hInvMassProtonXi_9 ->Fill(invMass);
        }
        if((multiplicity>=37) && (multiplicity<=40)){
          hInvMassProtonXi_10 ->Fill(invMass);
        }
        if((multiplicity>40) && (multiplicity<=50)){
          hInvMassProtonXi_11 ->Fill(invMass);
        }
	if((multiplicity>50) && (multiplicity<=60)){
          hInvMassProtonXi_12 ->Fill(invMass);
        }
        if((multiplicity>60) && (multiplicity<=70)){
          hInvMassProtonXi_13 ->Fill(invMass);
        }
        if((multiplicity>70) && (multiplicity<=80)){
          hInvMassProtonXi_14 ->Fill(invMass);
        }
        if((multiplicity>80) && (multiplicity<=90)){
          hInvMassProtonXi_15 ->Fill(invMass);
        }
        if((multiplicity>90) && (multiplicity<=100)){
          hInvMassProtonXi_16 ->Fill(invMass);
        }
        if(multiplicity>100){
          hInvMassProtonXi_17 ->Fill(invMass);
        }


	if((invMass > 2.25) && (invMass < 2.27)){
	  hInvMassProtonXi_lowmassx->Fill(xipx,prpx); 
	  hInvMassProtonXi_lowmassy->Fill(xipy,prpy); 
	  hInvMassProtonXi_lowmassz->Fill(xipz,prpz); 
	}

      }
    }
  }

  // pbar + Xi+
  if((npb > 0) && (nXip > 0)){
    for(Int_t i=0; i < fCascadeArrayp->GetEntriesFast(); i++){
      TObject *objCascade=fCascadeArrayp->At(i);    
      if(!objCascade) continue;
      AliAODcascade *xi=dynamic_cast<AliAODcascade*>(objCascade); 

      for(Int_t j=0; j < fProtonArrayb->GetEntriesFast(); j++){
	TObject *objProton=fProtonArrayb->At(j);    
	if(!objProton) continue;
	AliAODTrack *proton=dynamic_cast<AliAODTrack*>(objProton); 

	xipx        =xi->MomXiX();
	xipy        =xi->MomXiY();
	xipz        =xi->MomXiZ();
	energyxi    =sqrt(pow(1.322,2)+pow(xipx,2)+pow(xipy,2)+pow(xipz,2));
	bachid      =xi->GetBachID();
	ptrackid    =xi->GetPosID();
	ntrackid    =xi->GetNegID();
	prpx        =proton->Px();
	prpy        =proton->Py();
	prpz        =proton->Pz();
	energypr    =sqrt(pow(0.938,2)+pow(prpx,2)+pow(prpy,2)+pow(prpz,2));
	protonid    =proton->GetID();
	if(protonid < 0) protonid =-protonid-1;
	energysum   =energyxi+energypr;

	relmom    =sqrt(pow((xipx-prpx),2)+pow((xipy-prpy),2)+pow((xipz-prpz),2))/2;
	openangle =OpenAngle(xipx,xipy,xipz,prpx,prpy,prpz)*57.2957796;
	pt        =sqrt((prpx+xipx)*(prpx+xipx)+(prpy+xipy)*(prpy+xipy));
	invMass   =InvariantMass(xipx,xipy,xipz,prpx,prpy,prpz,energysum);

	if(protonid == ptrackid) continue;
	if(protonid == ntrackid) continue;
	if(protonid == bachid) continue;
	if(bachid == ptrackid) continue;
	if(bachid == ntrackid) continue;
	if(ptrackid == ntrackid) continue;

	hInvMassantiProtonXi             ->Fill(pt,invMass); 
	hInvMassantiProtonXi_relmomentum ->Fill(relmom,invMass);
	hInvMassantiProtonXi_openangle   ->Fill(openangle,invMass);
	fmultiplicity_antiprotonxi       ->Fill(multiplicity);

	isantiprotonxi=true; 

	//multiplicity
	if((multiplicity>=1) && (multiplicity<=4)){
          hInvMassAntiProtonXi_1 ->Fill(invMass);
        }
        if((multiplicity>=5) && (multiplicity<=8)){
          hInvMassAntiProtonXi_2 ->Fill(invMass);
        }
        if((multiplicity>=9) && (multiplicity<=12)){
          hInvMassAntiProtonXi_3 ->Fill(invMass);
        }
        if((multiplicity>=13) && (multiplicity<=16)){
          hInvMassAntiProtonXi_4 ->Fill(invMass);
        }
        if((multiplicity>=17) && (multiplicity<=20)){
          hInvMassAntiProtonXi_5 ->Fill(invMass);
        }
        if((multiplicity>=21) && (multiplicity<=24)){
          hInvMassAntiProtonXi_6 ->Fill(invMass);
        }
        if((multiplicity>=25) && (multiplicity<=28)){
          hInvMassAntiProtonXi_7 ->Fill(invMass);
        }
        if((multiplicity>=29) && (multiplicity<=32)){
          hInvMassAntiProtonXi_8 ->Fill(invMass);
        }
        if((multiplicity>=33) && (multiplicity<=36)){
          hInvMassAntiProtonXi_9 ->Fill(invMass);
        }
	if((multiplicity>=37) && (multiplicity<=40)){
          hInvMassAntiProtonXi_10 ->Fill(invMass);
        }
        if((multiplicity>40) && (multiplicity<=50)){
          hInvMassAntiProtonXi_11 ->Fill(invMass);
        }
        if((multiplicity>50) && (multiplicity<=60)){
          hInvMassAntiProtonXi_12 ->Fill(invMass);
        }
        if((multiplicity>60) && (multiplicity<=70)){
          hInvMassAntiProtonXi_13 ->Fill(invMass);
        }
        if((multiplicity>70) && (multiplicity<=80)){
          hInvMassAntiProtonXi_14 ->Fill(invMass);
        }
        if((multiplicity>80) && (multiplicity<=90)){
          hInvMassAntiProtonXi_15 ->Fill(invMass);
        }
        if((multiplicity>90) && (multiplicity<=100)){
          hInvMassAntiProtonXi_16 ->Fill(invMass);
        }
        if(multiplicity>100){
          hInvMassAntiProtonXi_17 ->Fill(invMass);
        }

	if((invMass > 2.25) && (invMass < 2.27)){
	  hInvMassantiProtonXi_lowmassx->Fill(xipx,prpx); 
	  hInvMassantiProtonXi_lowmassy->Fill(xipy,prpy); 
	  hInvMassantiProtonXi_lowmassz->Fill(xipz,prpz); 
	}
      }
    }
  }

  //lambda+lambda
  if(nLamn>0){
    for(Int_t i=0; i < fV0Arrayn->GetEntriesFast()-1; i++){
      TObject *objV01=fV0Arrayn->At(i);    
      if(!objV01) continue;
      AliAODv0 *lambda1=dynamic_cast<AliAODv0*>(objV01); 
      AliAODTrack* ptrack1=dynamic_cast<AliAODTrack*>(lambda1->GetDaughter(0));    
      for(Int_t j=i+1; j < fV0Arrayn->GetEntriesFast(); j++){
	TObject *objV02=fV0Arrayn->At(j);    
	if(!objV02) continue;
	AliAODv0 *lambda2=dynamic_cast<AliAODv0*>(objV02); 
	AliAODTrack* ptrack2=dynamic_cast<AliAODTrack*>(lambda2->GetDaughter(0));    

	protonpt1 = ptrack1->Pt();
	protonpt2 = ptrack2->Pt();

	l1px     =lambda1->MomV0X();
	l1py     =lambda1->MomV0Y();
	l1pz     =lambda1->MomV0Z();
	l2px     =lambda2->MomV0X();
	l2py     =lambda2->MomV0Y();
	l2pz     =lambda2->MomV0Z();
	energy1  =sqrt(pow(1.115,2)+pow(l1px,2)+pow(l1py,2)+pow(l1pz,2));  
	energy2  =sqrt(pow(1.115,2)+pow(l2px,2)+pow(l2py,2)+pow(l2pz,2));  
	posID1   =lambda1->GetPosID();
	negID1   =lambda1->GetNegID();
	posID2   =lambda2->GetPosID();
	negID2   =lambda2->GetNegID();
	energysum=energy1+energy2;
	pt        =sqrt(pow(l1px+l2px,2)+pow(l1py+l2py,2));
	relmom    =sqrt(pow((l1px-l2px),2)+pow((l1py-l2py),2)+pow((l1pz-l2pz),2))/2;
	openangle =OpenAngle(l1px,l1py,l1pz,l2px,l2py,l2pz)*57.2957796;
	invMass   =InvariantMass(l1px,l1py,l1pz,l2px,l2py,l2pz,energysum);
	isXi1     =Lambda[i];
	isXi2     =Lambda[j];	

	if(posID1==posID2) continue;  
        if(negID1==negID2) continue;    
	
	// all lambda + all lambda
	hInvMassLambdaLambda_all        ->Fill(pt,invMass);
	hInvMassLambdaLambda_relmomentum->Fill(relmom,invMass);
	hInvMassLambdaLambda_openangle  ->Fill(openangle,invMass);
	hAngle_relmom                   ->Fill(relmom,openangle);
	hPtcorrelation_lamlam           ->Fill(protonpt1,protonpt2);
	fmultiplicity_lambdalambda      ->Fill(multiplicity);
	
	islambdalambda=true; 

	//multiplicity
	if((multiplicity>=1) && (multiplicity<=4)){
          hInvMassLambdaLambda_1 ->Fill(invMass);
        }
        if((multiplicity>=5) && (multiplicity<=8)){
          hInvMassLambdaLambda_2 ->Fill(invMass);
        }
        if((multiplicity>=9) && (multiplicity<=12)){
          hInvMassLambdaLambda_3 ->Fill(invMass);
        }
        if((multiplicity>=13) && (multiplicity<=16)){
          hInvMassLambdaLambda_4 ->Fill(invMass);
        }
        if((multiplicity>=17) && (multiplicity<=20)){
          hInvMassLambdaLambda_5 ->Fill(invMass);
        }
        if((multiplicity>=21) && (multiplicity<=24)){
          hInvMassLambdaLambda_6 ->Fill(invMass);
        }
        if((multiplicity>=25) && (multiplicity<=28)){
          hInvMassLambdaLambda_7 ->Fill(invMass);
        }
        if((multiplicity>=29) && (multiplicity<=32)){
          hInvMassLambdaLambda_8 ->Fill(invMass);
        }
        if((multiplicity>=33) && (multiplicity<=36)){
          hInvMassLambdaLambda_9 ->Fill(invMass);
        }
        if((multiplicity>=37) && (multiplicity<=40)){
          hInvMassLambdaLambda_10 ->Fill(invMass);
        }
        if((multiplicity>40) && (multiplicity<=50)){
          hInvMassLambdaLambda_11 ->Fill(invMass);
        }
        if((multiplicity>50) && (multiplicity<=60)){
          hInvMassLambdaLambda_12 ->Fill(invMass);
        }
        if((multiplicity>60) && (multiplicity<=70)){
          hInvMassLambdaLambda_13 ->Fill(invMass);
        }
	if((multiplicity>70) && (multiplicity<=80)){
          hInvMassLambdaLambda_14 ->Fill(invMass);
        }
        if((multiplicity>80) && (multiplicity<=90)){
          hInvMassLambdaLambda_15 ->Fill(invMass);
        }
        if((multiplicity>90) && (multiplicity<=100)){
          hInvMassLambdaLambda_16 ->Fill(invMass);
        }
        if(multiplicity>100){
          hInvMassLambdaLambda_17 ->Fill(invMass);
        }

	// prompt lambda + prompt lambda
	if((isXi1==0) && (isXi2==0)){
	  hInvMassLambdaLambda_onlyprompt->Fill(pt,invMass);
	}	
      }
    }
  }

  // antilambda + antilambda
  if(nLamb>0){
    for(Int_t i=0; i < fV0Arrayb->GetEntriesFast()-1; i++){
      TObject *objV01=fV0Arrayb->At(i);    
      if(!objV01) continue;
      AliAODv0 *lambda1=dynamic_cast<AliAODv0*>(objV01); 
      for(Int_t j=i+1; j < fV0Arrayb->GetEntriesFast(); j++){
	TObject *objV02=fV0Arrayb->At(j);    
	if(!objV02) continue;
	AliAODv0 *lambda2=dynamic_cast<AliAODv0*>(objV02); 
	
	l1px     =lambda1->MomV0X();
	l1py     =lambda1->MomV0Y();
	l1pz     =lambda1->MomV0Z();
	l2px     =lambda2->MomV0X();
	l2py     =lambda2->MomV0Y();
	l2pz     =lambda2->MomV0Z();
	energy1  =sqrt(pow(1.115,2)+pow(l1px,2)+pow(l1py,2)+pow(l1pz,2));  
	energy2  =sqrt(pow(1.115,2)+pow(l2px,2)+pow(l2py,2)+pow(l2pz,2));  
	posID1   =lambda1->GetPosID();
	negID1   =lambda1->GetNegID();
	posID2   =lambda2->GetPosID();
	negID2   =lambda2->GetNegID();
	energysum=energy1+energy2;
	pt        =sqrt(pow(l1px+l2px,2)+pow(l1py+l2py,2));
	relmom    =sqrt(pow((l1px-l2px),2)+pow((l1py-l2py),2)+pow((l1pz-l2pz),2))/2;
	openangle =OpenAngle(l1px,l1py,l1pz,l2px,l2py,l2pz)*57.2957796;
	invMass   =InvariantMass(l1px,l1py,l1pz,l2px,l2py,l2pz,energysum);
	isXi1     =AntiLambda[i];
	isXi2     =AntiLambda[j];	

	if(posID1==posID2) continue;  
        if(negID1==negID2) continue;    
	
	// all lambda + all lambda
	hInvMassAntiLambdaAntiLambda_all        ->Fill(pt,invMass);
	hInvMassAntiLambdaAntiLambda_relmomentum->Fill(relmom,invMass);
	hInvMassAntiLambdaAntiLambda_openangle  ->Fill(openangle,invMass);
	fmultiplicity_antilambdaantilambda      ->Fill(multiplicity);

	isantilambdaantilambda=true; 

	//multiplicity
	if((multiplicity>=1) && (multiplicity<=4)){
          hInvMassAntiLambdaAntiLambda_1 ->Fill(invMass);
        }
        if((multiplicity>=5) && (multiplicity<=8)){
          hInvMassAntiLambdaAntiLambda_2 ->Fill(invMass);
        }
        if((multiplicity>=9) && (multiplicity<=12)){
          hInvMassAntiLambdaAntiLambda_3 ->Fill(invMass);
        }
        if((multiplicity>=13) && (multiplicity<=16)){
          hInvMassAntiLambdaAntiLambda_4 ->Fill(invMass);
        }
        if((multiplicity>=17) && (multiplicity<=20)){
          hInvMassAntiLambdaAntiLambda_5 ->Fill(invMass);
        }
        if((multiplicity>=21) && (multiplicity<=24)){
          hInvMassAntiLambdaAntiLambda_6 ->Fill(invMass);
        }
        if((multiplicity>=25) && (multiplicity<=28)){
          hInvMassAntiLambdaAntiLambda_7 ->Fill(invMass);
        }
        if((multiplicity>=29) && (multiplicity<=32)){
          hInvMassAntiLambdaAntiLambda_8 ->Fill(invMass);
        }
        if((multiplicity>=33) && (multiplicity<=36)){
          hInvMassAntiLambdaAntiLambda_9 ->Fill(invMass);
        }
        if((multiplicity>=37) && (multiplicity<=40)){
          hInvMassAntiLambdaAntiLambda_10 ->Fill(invMass);
        }
        if((multiplicity>40) && (multiplicity<=50)){
          hInvMassAntiLambdaAntiLambda_11 ->Fill(invMass);
        }
        if((multiplicity>50) && (multiplicity<=60)){
          hInvMassAntiLambdaAntiLambda_12 ->Fill(invMass);
        }
        if((multiplicity>60) && (multiplicity<=70)){
          hInvMassAntiLambdaAntiLambda_13 ->Fill(invMass);
        }
        if((multiplicity>70) && (multiplicity<=80)){
          hInvMassAntiLambdaAntiLambda_14 ->Fill(invMass);
        }
	if((multiplicity>80) && (multiplicity<=90)){
          hInvMassAntiLambdaAntiLambda_15 ->Fill(invMass);
        }
        if((multiplicity>90) && (multiplicity<=100)){
          hInvMassAntiLambdaAntiLambda_16 ->Fill(invMass);
        }
        if(multiplicity>100){
          hInvMassAntiLambdaAntiLambda_17 ->Fill(invMass);
        }

	// prompt lambda + prompt lambda
	if((isXi1==0) && (isXi2==0)){
	  hInvMassAntiLambdaAntiLambda_onlyprompt->Fill(pt,invMass);
	}
      }
    }
  }

  // lambda + xi-
  if((nLamn > 0) && (nXim > 0)){
    for(Int_t i=0; i < fCascadeArraym->GetEntriesFast(); i++){
      TObject *objCascade=fCascadeArraym->At(i);    
      if(!objCascade) continue;
      AliAODcascade *xi=dynamic_cast<AliAODcascade*>(objCascade); 
      AliAODTrack *ptrkxi=dynamic_cast<AliAODTrack*>(xi->GetDaughter(0)); 

      for(Int_t j=0; j < fV0Arrayn->GetEntriesFast(); j++){
	TObject *objV0=fV0Arrayn->At(j);    
	if(!objV0) continue;
	AliAODv0 *lambda=dynamic_cast<AliAODv0*>(objV0); 
	AliAODTrack *ptrack=dynamic_cast<AliAODTrack*>(lambda->GetDaughter(0));    

	xidecaypt   =ptrkxi->Pt();
	protonpt1   =ptrack->Pt();

	xipx        =xi->MomXiX();
	xipy        =xi->MomXiY();
	xipz        =xi->MomXiZ();
	energyxi    =sqrt(pow(1.322,2)+pow(xipx,2)+pow(xipy,2)+pow(xipz,2));
	bachid      =xi->GetBachID();
	ptrackid    =xi->GetPosID();
	ntrackid    =xi->GetNegID();
	lampx       =lambda->MomV0X();
	lampy       =lambda->MomV0Y();
	lampz       =lambda->MomV0Z();
	posID       =lambda->GetPosID();
	negID       =lambda->GetNegID();
	energylam   =sqrt(pow(1.115,2)+pow(lampx,2)+pow(lampy,2)+pow(lampz,2));  
	energysum   =energyxi+energylam;

	relmom    =sqrt(pow((xipx-lampx),2)+pow((xipy-lampy),2)+pow((xipz-lampz),2))/2;
	openangle =OpenAngle(xipx,xipy,xipz,lampx,lampy,lampz)*57.2957796;
	pt          =sqrt(pow((lampx+xipx),2)+pow((lampy+xipy),2));
      	invMass     =InvariantMass(xipx,xipy,xipz,lampx,lampy,lampz,energysum);	

      	if(posID == negID) continue;
	if(posID == ptrackid) continue;
	if(posID == ntrackid) continue;
	if(posID == bachid) continue;
	if(negID == ptrackid) continue;
	if(negID == ntrackid) continue;
	if(negID == bachid) continue;
	if(ptrackid == ntrackid) continue;
	if(ptrackid == bachid) continue;
	if(ntrackid == bachid) continue;

	hInvMassLambdaXi            ->Fill(pt,invMass);
	hPtcorrelation_lamxi        ->Fill(xidecaypt,protonpt1);
	hInvMassLambdaXi_relmomentum->Fill(relmom,invMass);
	hInvMassLambdaXi_openangle  ->Fill(openangle,invMass);
	fmultiplicity_lambdaxi      ->Fill(multiplicity);

	islambdaxi=true; 

	//multiplicity
	if((multiplicity>=1) && (multiplicity<=5)){
          hInvMassLambdaXi_1 ->Fill(invMass);
        }
        if((multiplicity>=6) && (multiplicity<=10)){
          hInvMassLambdaXi_2 ->Fill(invMass);
        }
        if((multiplicity>=11) && (multiplicity<=15)){
          hInvMassLambdaXi_3 ->Fill(invMass);
        }
        if((multiplicity>=16) && (multiplicity<=20)){
          hInvMassLambdaXi_4 ->Fill(invMass);
        }
        if((multiplicity>=21) && (multiplicity<=25)){
          hInvMassLambdaXi_5 ->Fill(invMass);
        }
        if((multiplicity>=26) && (multiplicity<=30)){
          hInvMassLambdaXi_6 ->Fill(invMass);
        }
        if((multiplicity>=31) && (multiplicity<=35)){
          hInvMassLambdaXi_7 ->Fill(invMass);
        }
        if((multiplicity>=36) && (multiplicity<=40)){
          hInvMassLambdaXi_8 ->Fill(invMass);
        }
        if((multiplicity>40) && (multiplicity<=50)){
          hInvMassLambdaXi_9 ->Fill(invMass);
        }
        if((multiplicity>50) && (multiplicity<=60)){
          hInvMassLambdaXi_10 ->Fill(invMass);
        }
	if((multiplicity>60) && (multiplicity<=70)){
          hInvMassLambdaXi_11 ->Fill(invMass);
        }
        if((multiplicity>70) && (multiplicity<=80)){
          hInvMassLambdaXi_12 ->Fill(invMass);
        }
        if((multiplicity>80) && (multiplicity<=90)){
          hInvMassLambdaXi_13 ->Fill(invMass);
        }
        if((multiplicity>90) && (multiplicity<=100)){
          hInvMassLambdaXi_14 ->Fill(invMass);
        }
        if(multiplicity>100){
          hInvMassLambdaXi_15 ->Fill(invMass);
        }

      }
    }
  }

  // antilambda + xi+
  if((nLamb > 0) && (nXip > 0)){
    for(Int_t i=0; i < fCascadeArrayp->GetEntriesFast(); i++){
      TObject *objCascade=fCascadeArrayp->At(i);    
      if(!objCascade) continue;
      AliAODcascade *xi=dynamic_cast<AliAODcascade*>(objCascade); 

      for(Int_t j=0; j < fV0Arrayb->GetEntriesFast(); j++){
	TObject *objV0=fV0Arrayb->At(j);    
	if(!objV0) continue;
	AliAODv0 *lambda=dynamic_cast<AliAODv0*>(objV0); 

	xipx        =xi->MomXiX();
	xipy        =xi->MomXiY();
	xipz        =xi->MomXiZ();
	energyxi    =sqrt(pow(1.322,2)+pow(xipx,2)+pow(xipy,2)+pow(xipz,2));
	bachid      =xi->GetBachID();
	ptrackid    =xi->GetPosID();
	ntrackid    =xi->GetNegID();
	lampx       =lambda->MomV0X();
	lampy       =lambda->MomV0Y();
	lampz       =lambda->MomV0Z();
	posID       =lambda->GetPosID();
	negID       =lambda->GetNegID();
	energylam   =sqrt(pow(1.115,2)+pow(lampx,2)+pow(lampy,2)+pow(lampz,2));  
	energysum   =energyxi+energylam;

	relmom    =sqrt(pow((xipx-lampx),2)+pow((xipy-lampy),2)+pow((xipz-lampz),2))/2;
	openangle =OpenAngle(xipx,xipy,xipz,lampx,lampy,lampz)*57.2957796;
	pt          =sqrt(pow((lampx+xipx),2)+pow((lampy+xipy),2));
      	invMass     =InvariantMass(xipx,xipy,xipz,lampx,lampy,lampz,energysum);	

      	if(posID == negID) continue;
	if(posID == ptrackid) continue;
	if(posID == ntrackid) continue;
	if(posID == bachid) continue;
	if(negID == ptrackid) continue;
	if(negID == ntrackid) continue;
	if(negID == bachid) continue;
	if(ptrackid == ntrackid) continue;
	if(ptrackid == bachid) continue;
	if(ntrackid == bachid) continue;

	hInvMassAntiLambdaXip           ->Fill(pt,invMass);
	hInvMassAntiLambdaXi_relmomentum->Fill(relmom,invMass);
	hInvMassAntiLambdaXi_openangle  ->Fill(openangle,invMass);
	fmultiplicity_antilambdaxi      ->Fill(multiplicity);

	isantilambdaxi=true; 

	//multiplicity
	if((multiplicity>=1) && (multiplicity<=5)){
          hInvMassAntiLambdaXi_1 ->Fill(invMass);
        }
        if((multiplicity>=6) && (multiplicity<=10)){
          hInvMassAntiLambdaXi_2 ->Fill(invMass);
        }
        if((multiplicity>=11) && (multiplicity<=15)){
          hInvMassAntiLambdaXi_3 ->Fill(invMass);
        }
        if((multiplicity>=16) && (multiplicity<=20)){
          hInvMassAntiLambdaXi_4 ->Fill(invMass);
        }
        if((multiplicity>=21) && (multiplicity<=25)){
          hInvMassAntiLambdaXi_5 ->Fill(invMass);
        }
        if((multiplicity>=26) && (multiplicity<=30)){
          hInvMassAntiLambdaXi_6 ->Fill(invMass);
        }
        if((multiplicity>=31) && (multiplicity<=35)){
          hInvMassAntiLambdaXi_7 ->Fill(invMass);
        }
        if((multiplicity>=36) && (multiplicity<=40)){
          hInvMassAntiLambdaXi_8 ->Fill(invMass);
        }
        if((multiplicity>40) && (multiplicity<=50)){
          hInvMassAntiLambdaXi_9 ->Fill(invMass);
        }
        if((multiplicity>50) && (multiplicity<=60)){
          hInvMassAntiLambdaXi_10 ->Fill(invMass);
        }
	if((multiplicity>60) && (multiplicity<=70)){
          hInvMassAntiLambdaXi_11 ->Fill(invMass);
        }
        if((multiplicity>70) && (multiplicity<=80)){
          hInvMassAntiLambdaXi_12 ->Fill(invMass);
        }
        if((multiplicity>80) && (multiplicity<=90)){
          hInvMassAntiLambdaXi_13 ->Fill(invMass);
        }
        if((multiplicity>90) && (multiplicity<=100)){
          hInvMassAntiLambdaXi_14 ->Fill(invMass);
        }
        if(multiplicity>100){
          hInvMassAntiLambdaXi_15 ->Fill(invMass);
        }

      }
    }
  }

  //========== Event mixing
  fmultiplicity->Fill(multiplicity,vecTarget[2]);
  /*
  if(0 < nLamn){

    Double_t z_min = -10; 
    Double_t z_max = 10; 
    Double_t m_min = 0;
    Double_t m_max = 200;
    Double_t z_divide = 4;
    Double_t m_divide = 20;
    Double_t z_width = ((z_max)-(z_min))/(z_divide); 
    Double_t m_width = ((m_max)-(m_min))/(m_divide); 
    
    Int_t iz    = (vecTarget[2]-(z_min))/z_width;
    Int_t im    = (multiplicity-(m_min))/m_width;

    for(Int_t i=0; i < nLamn; i++){ 
      
      for(Int_t j=0; j < eventdepth2[iz][im]; j++){
	for(Int_t k=0; k < mnLambda[iz][im][j]; k++){
	  
	  Double_t energysum=0.,l1px=0.,l1py=0.,l1pz=0.,l2px=0.,l2py=0.,l2pz=0.,pt=0.,invMass=0.;
	  
	  energysum = lamntrack[0][i]+mEnergy[iz][im][j][k];
	  l1px      = lamntrack[1][i];
	  l1py      = lamntrack[2][i];
	  l1pz      = lamntrack[3][i];
	  l2px      = mPx[iz][im][j][k];
	  l2py      = mPy[iz][im][j][k];
	  l2pz      = mPz[iz][im][j][k];
	  
	  pt        =sqrt(pow(l1px+l2px,2)+pow(l1py+l2py,2));
	  invMass   =InvariantMass(l1px,l1py,l1pz,l2px,l2py,l2pz,energysum);       
	  
	  hInvMassLambdaLambda_combinatorial->Fill(pt,invMass);
	  
	}
      }
    }
    
    mnLambda[iz][im][eventdepth[iz][im]] = nLamn;
        
    for(Int_t i=0; i < nLamn; i++){
      
      mEnergy[iz][im][eventdepth[iz][im]][i] = lamntrack[0][i];
      mPx[iz][im][eventdepth[iz][im]][i] = lamntrack[1][i];
      mPy[iz][im][eventdepth[iz][im]][i] = lamntrack[2][i];
      mPz[iz][im][eventdepth[iz][im]][i] = lamntrack[3][i];
      
    }
    
    if(eventdepth[iz][im] < 5){
      eventdepth[iz][im]++;
    }
    if(eventdepth[iz][im] == 5){       
      eventdepth[iz][im] = 0;
    }    
    if(eventdepth2[iz][im] < 5){
      eventdepth2[iz][im]++;
    }    
    if((eventdepth2[iz][im] == 5) || (eventdepth2[iz][im] > 5)){
      eventdepth2[iz][im] = 5;
    }
  }
  */

  //========== eventmixing with AliEventpool
  AliEventPool *lambdapool =fPoolManagerlambda->GetEventPool(multiplicity,vecTarget[2],0.);
  if(!lambdapool) return; 
  AliEventPool *antilambdapool =fPoolManagerantilambda->GetEventPool(multiplicity,vecTarget[2],0.);
  if(!antilambdapool) return; 

  AliEventPool *xipool =fPoolManagerxi->GetEventPool(multiplicity,vecTarget[2],0.);
  if(!xipool) return; 
  AliEventPool *xippool =fPoolManagerxip->GetEventPool(multiplicity,vecTarget[2],0.);
  if(!xippool) return; 
  AliEventPool *protonpool =fPoolManagerproton->GetEventPool(multiplicity,vecTarget[2],0.);
  if(!protonpool) return; 
  AliEventPool *antiprotonpool =fPoolManagerantiproton->GetEventPool(multiplicity,vecTarget[2],0.);
  if(!antiprotonpool) return; 

  AliEventPool *lampool =fPoolManagerlam->GetEventPool(multiplicity,vecTarget[2],0.);
  if(!lampool) return; 
  AliEventPool *antilampool =fPoolManagerantilam->GetEventPool(multiplicity,vecTarget[2],0.);
  if(!antilampool) return; 
  AliEventPool *xpool =fPoolManagerx->GetEventPool(multiplicity,vecTarget[2],0.);
  if(!xpool) return; 
  AliEventPool *xppool =fPoolManagerxp->GetEventPool(multiplicity,vecTarget[2],0.);
  if(!xppool) return; 

  // lambda + lambda
  //if(islambdalambda){
  for(Int_t iv0=0; iv0 < fV0Arrayn->GetEntriesFast(); iv0++){
    TObject *objV0=fV0Arrayn->At(iv0);    
    if(!objV0) continue;
    AliAODv0 *lambda1=dynamic_cast<AliAODv0*>(objV0); 
    l1px     =lambda1->MomV0X();
    l1py     =lambda1->MomV0Y();
    l1pz     =lambda1->MomV0Z();
    energy1  =sqrt(pow(1.115,2)+pow(l1px,2)+pow(l1py,2)+pow(l1pz,2));  
    
    //mixed event loop starts
    for(Int_t imevt=0; imevt < lambdapool->GetCurrentNEvents(); imevt++){
      TObjArray *poolv0=(TObjArray*)lambdapool->GetEvent(imevt);
      for(Int_t imv0=0; imv0<poolv0->GetEntriesFast(); imv0++){
	TObject *objpV0=poolv0->At(imv0);
	if(!objpV0) continue;
	AliAODv0 *lambda2=dynamic_cast<AliAODv0*>(objpV0);
	l2px     =lambda2->MomV0X();
	l2py     =lambda2->MomV0Y();
	l2pz     =lambda2->MomV0Z();
	energy2  =sqrt(pow(1.115,2)+pow(l2px,2)+pow(l2py,2)+pow(l2pz,2));
	energysum=energy1+energy2;
	pt        =sqrt(pow(l1px+l2px,2)+pow(l1py+l2py,2));
	invMass   =InvariantMass(l1px,l1py,l1pz,l2px,l2py,l2pz,energysum);
	hInvMassLambdaLambda_evtpool->Fill(pt,invMass);
	relmom    =sqrt(pow((l1px-l2px),2)+pow((l1py-l2py),2)+pow((l1pz-l2pz),2))/2;
	
	//multiplicity
	if((multiplicity>=1) && (multiplicity<=4)){
	  hInvMassLambdaLambda_evtpool1 ->Fill(invMass);
	}
	if((multiplicity>=5) && (multiplicity<=8)){
	  hInvMassLambdaLambda_evtpool2 ->Fill(invMass);
	}
	if((multiplicity>=9) && (multiplicity<=12)){
	  hInvMassLambdaLambda_evtpool3 ->Fill(invMass);
	}
	if((multiplicity>=13) && (multiplicity<=16)){
	  hInvMassLambdaLambda_evtpool4 ->Fill(invMass);
	}
	if((multiplicity>=17) && (multiplicity<=20)){
	  hInvMassLambdaLambda_evtpool5 ->Fill(invMass);
	}
	if((multiplicity>=21) && (multiplicity<=24)){
	  hInvMassLambdaLambda_evtpool6 ->Fill(invMass);
	}
	if((multiplicity>=25) && (multiplicity<=28)){
	  hInvMassLambdaLambda_evtpool7 ->Fill(invMass);
	}
	if((multiplicity>=29) && (multiplicity<=32)){
	  hInvMassLambdaLambda_evtpool8 ->Fill(invMass);
	}
	if((multiplicity>=33) && (multiplicity<=36)){
	  hInvMassLambdaLambda_evtpool9 ->Fill(invMass);
	}
	if((multiplicity>=37) && (multiplicity<=40)){
	  hInvMassLambdaLambda_evtpool10 ->Fill(invMass);
	}
	if((multiplicity>40) && (multiplicity<=50)){
	  hInvMassLambdaLambda_evtpool11 ->Fill(invMass);
	}
	if((multiplicity>50) && (multiplicity<=60)){
	  hInvMassLambdaLambda_evtpool12 ->Fill(invMass);
	}
	if((multiplicity>60) && (multiplicity<=70)){
	  hInvMassLambdaLambda_evtpool13 ->Fill(invMass);
	}
	if((multiplicity>70) && (multiplicity<=80)){
	  hInvMassLambdaLambda_evtpool14 ->Fill(invMass);
	}
	if((multiplicity>80) && (multiplicity<=90)){
	  hInvMassLambdaLambda_evtpool15 ->Fill(invMass);
	}
	if((multiplicity>90) && (multiplicity<=100)){
	  hInvMassLambdaLambda_evtpool16 ->Fill(invMass);
	}
	if(multiplicity>100){
	  hInvMassLambdaLambda_evtpool17 ->Fill(invMass);
	}
	
      }
    }   
  }  
  //}

  // antilambda + antilambda
  //if(isantilambdaantilambda){
  for(Int_t iv0=0; iv0 < fV0Arrayb->GetEntriesFast(); iv0++){
    TObject *objV0=fV0Arrayb->At(iv0);    
    if(!objV0) continue;
    AliAODv0 *lambda1=dynamic_cast<AliAODv0*>(objV0); 
    l1px     =lambda1->MomV0X();
    l1py     =lambda1->MomV0Y();
    l1pz     =lambda1->MomV0Z();
    energy1  =sqrt(pow(1.115,2)+pow(l1px,2)+pow(l1py,2)+pow(l1pz,2));  
    
    //mixed event loop starts
    for(Int_t imevt=0; imevt < antilambdapool->GetCurrentNEvents(); imevt++){
      TObjArray *poolv0=(TObjArray*)antilambdapool->GetEvent(imevt);
      for(Int_t imv0=0; imv0<poolv0->GetEntriesFast(); imv0++){
	TObject *objpV0=poolv0->At(imv0);
	if(!objpV0) continue;
	AliAODv0 *lambda2=dynamic_cast<AliAODv0*>(objpV0);
	l2px     =lambda2->MomV0X();
	l2py     =lambda2->MomV0Y();
	l2pz     =lambda2->MomV0Z();
	energy2  =sqrt(pow(1.115,2)+pow(l2px,2)+pow(l2py,2)+pow(l2pz,2));
	energysum=energy1+energy2;
	pt        =sqrt(pow(l1px+l2px,2)+pow(l1py+l2py,2));
	invMass   =InvariantMass(l1px,l1py,l1pz,l2px,l2py,l2pz,energysum);
	hInvMassAntiLambdaAntiLambda_evtpool->Fill(pt,invMass);       
	relmom    =sqrt(pow((l1px-l2px),2)+pow((l1py-l2py),2)+pow((l1pz-l2pz),2))/2;	
	
	//multiplicity
	if((multiplicity>=1) && (multiplicity<=4)){
	  hInvMassAntiLambdaAntiLambda_evtpool1 ->Fill(invMass);
	}
	if((multiplicity>=5) && (multiplicity<=8)){
	  hInvMassAntiLambdaAntiLambda_evtpool2 ->Fill(invMass);
	}
	if((multiplicity>=9) && (multiplicity<=12)){
	  hInvMassAntiLambdaAntiLambda_evtpool3 ->Fill(invMass);
	}
	if((multiplicity>=13) && (multiplicity<=16)){
	  hInvMassAntiLambdaAntiLambda_evtpool4 ->Fill(invMass);
	}
	if((multiplicity>=17) && (multiplicity<=20)){
	  hInvMassAntiLambdaAntiLambda_evtpool5 ->Fill(invMass);
	}
	if((multiplicity>=21) && (multiplicity<=24)){
	  hInvMassAntiLambdaAntiLambda_evtpool6 ->Fill(invMass);
	}
	if((multiplicity>=25) && (multiplicity<=28)){
          hInvMassAntiLambdaAntiLambda_evtpool7 ->Fill(invMass);
	}
	if((multiplicity>=29) && (multiplicity<=32)){
	  hInvMassAntiLambdaAntiLambda_evtpool8 ->Fill(invMass);
	}
	if((multiplicity>=33) && (multiplicity<=36)){
	  hInvMassAntiLambdaAntiLambda_evtpool9 ->Fill(invMass);
	}
	if((multiplicity>=37) && (multiplicity<=40)){
	  hInvMassAntiLambdaAntiLambda_evtpool10 ->Fill(invMass);
	}
	if((multiplicity>40) && (multiplicity<=50)){
	  hInvMassAntiLambdaAntiLambda_evtpool11 ->Fill(invMass);
	}
	if((multiplicity>50) && (multiplicity<=60)){
	  hInvMassAntiLambdaAntiLambda_evtpool12 ->Fill(invMass);
	}
	if((multiplicity>60) && (multiplicity<=70)){
	  hInvMassAntiLambdaAntiLambda_evtpool13 ->Fill(invMass);
	}
	if((multiplicity>70) && (multiplicity<=80)){
	  hInvMassAntiLambdaAntiLambda_evtpool14 ->Fill(invMass);
	}
	if((multiplicity>80) && (multiplicity<=90)){
	  hInvMassAntiLambdaAntiLambda_evtpool15 ->Fill(invMass);
	}
	if((multiplicity>90) && (multiplicity<=100)){
	  hInvMassAntiLambdaAntiLambda_evtpool16 ->Fill(invMass);
	}
	if(multiplicity>100){
	  hInvMassAntiLambdaAntiLambda_evtpool17 ->Fill(invMass);
	}
	
      }
    }   
  }   
  //}

  // p + xi-
  // proton
  //if(isprotonxi){
  for(Int_t itrk=0; itrk < fProtonArray->GetEntriesFast(); itrk++){
    TObject *objTrack=fProtonArray->At(itrk);    
    if(!objTrack) continue;
    AliAODTrack *proton=dynamic_cast<AliAODTrack*>(objTrack);    
    prpx        =proton->Px();
    prpy        =proton->Py();
    prpz        =proton->Pz();
    energypr    =sqrt(pow(0.938,2)+pow(prpx,2)+pow(prpy,2)+pow(prpz,2));
    
    //mixed event loop starts
    for(Int_t imevt=0; imevt < xipool->GetCurrentNEvents(); imevt++){
      TObjArray *poolcascade=(TObjArray*)xipool->GetEvent(imevt);
      for(Int_t imcascade=0; imcascade < poolcascade->GetEntriesFast(); imcascade++){
	TObject *objCascade=poolcascade->At(imcascade);
	if(!objCascade) continue;
	AliAODcascade *xi=dynamic_cast<AliAODcascade*>(objCascade); 
	xipx        =xi->MomXiX();
	xipy        =xi->MomXiY();
	xipz        =xi->MomXiZ();
	energyxi    =sqrt(pow(1.322,2)+pow(xipx,2)+pow(xipy,2)+pow(xipz,2));	
	
	pt        =sqrt((prpx+xipx)*(prpx+xipx)+(prpy+xipy)*(prpy+xipy));
	energysum = energypr+energyxi;
	invMass   =InvariantMass(xipx,xipy,xipz,prpx,prpy,prpz,energysum);
	hInvMassProtonXi_evtpool->Fill(pt,invMass);       
	relmom    =sqrt(pow((xipx-prpx),2)+pow((xipy-prpy),2)+pow((xipz-prpz),2))/2;
	
	//multiplicity
	if((multiplicity>=1) && (multiplicity<=4)){
	  hInvMassProtonXi_evtpool1 ->Fill(invMass);
	}
	if((multiplicity>=5) && (multiplicity<=8)){
	  hInvMassProtonXi_evtpool2 ->Fill(invMass);
	}
	if((multiplicity>=9) && (multiplicity<=12)){
	  hInvMassProtonXi_evtpool3 ->Fill(invMass);
	}
	if((multiplicity>=13) && (multiplicity<=16)){
	  hInvMassProtonXi_evtpool4 ->Fill(invMass);
	}
	if((multiplicity>=17) && (multiplicity<=20)){
	  hInvMassProtonXi_evtpool5 ->Fill(invMass);
	}
	if((multiplicity>=21) && (multiplicity<=24)){
	  hInvMassProtonXi_evtpool6 ->Fill(invMass);
	}
	if((multiplicity>=25) && (multiplicity<=28)){
	  hInvMassProtonXi_evtpool7 ->Fill(invMass);
	}
	if((multiplicity>=29) && (multiplicity<=32)){
	  hInvMassProtonXi_evtpool8 ->Fill(invMass);
	}
	if((multiplicity>=33) && (multiplicity<=36)){
	  hInvMassProtonXi_evtpool9 ->Fill(invMass);
	}
	if((multiplicity>=37) && (multiplicity<=40)){
	  hInvMassProtonXi_evtpool10 ->Fill(invMass);
	}
	if((multiplicity>40) && (multiplicity<=50)){
	  hInvMassProtonXi_evtpool11 ->Fill(invMass);
	}
	if((multiplicity>50) && (multiplicity<=60)){
	  hInvMassProtonXi_evtpool12 ->Fill(invMass);
	}
	if((multiplicity>60) && (multiplicity<=70)){
	  hInvMassProtonXi_evtpool13 ->Fill(invMass);
	}
	if((multiplicity>70) && (multiplicity<=80)){
	  hInvMassProtonXi_evtpool14 ->Fill(invMass);
	}
	if((multiplicity>80) && (multiplicity<=90)){
	  hInvMassProtonXi_evtpool15 ->Fill(invMass);
	}
	if((multiplicity>90) && (multiplicity<=100)){
	  hInvMassProtonXi_evtpool16 ->Fill(invMass);
	}
	if(multiplicity>100){
	  hInvMassProtonXi_evtpool17 ->Fill(invMass);
	}
	
      }
    }    
  } 
  // Xi
  for(Int_t icascade=0; icascade < fCascadeArraym->GetEntriesFast(); icascade++){
    TObject *objCascade=fCascadeArraym->At(icascade);    
    if(!objCascade) continue;
    AliAODcascade *xi=dynamic_cast<AliAODcascade*>(objCascade); 
    xipx        =xi->MomXiX();
    xipy        =xi->MomXiY();
    xipz        =xi->MomXiZ();
    energyxi    =sqrt(pow(1.322,2)+pow(xipx,2)+pow(xipy,2)+pow(xipz,2));
    
    //mixed event loop starts
    for(Int_t imevt=0; imevt < protonpool->GetCurrentNEvents(); imevt++){
      TObjArray *pooltrk=(TObjArray*)protonpool->GetEvent(imevt);      
      for(Int_t imtrk=0; imtrk < pooltrk->GetEntriesFast(); imtrk++){
	TObject *objProton=pooltrk->At(imtrk);
	if(!objProton) continue;
	AliAODTrack *proton=dynamic_cast<AliAODTrack*>(objProton); 	
	prpx        =proton->Px();
	prpy        =proton->Py();
	prpz        =proton->Pz();
	energypr    =sqrt(pow(0.938,2)+pow(prpx,2)+pow(prpy,2)+pow(prpz,2));
	
	energysum   =energyxi+energypr;
	pt        =sqrt((prpx+xipx)*(prpx+xipx)+(prpy+xipy)*(prpy+xipy));
	invMass   =InvariantMass(xipx,xipy,xipz,prpx,prpy,prpz,energysum);
	hInvMassProtonXi_evtpool->Fill(pt,invMass);       
	relmom    =sqrt(pow((xipx-prpx),2)+pow((xipy-prpy),2)+pow((xipz-prpz),2))/2;
	
	//multiplicity 
	if((multiplicity>=1) && (multiplicity<=4)){
	  hInvMassProtonXi_evtpool1 ->Fill(invMass);
	}
	if((multiplicity>=5) && (multiplicity<=8)){
	  hInvMassProtonXi_evtpool2 ->Fill(invMass);
	}
	if((multiplicity>=9) && (multiplicity<=12)){
	  hInvMassProtonXi_evtpool3 ->Fill(invMass);
	}
	if((multiplicity>=13) && (multiplicity<=16)){
	  hInvMassProtonXi_evtpool4 ->Fill(invMass);
	}
	if((multiplicity>=17) && (multiplicity<=20)){
	  hInvMassProtonXi_evtpool5 ->Fill(invMass);
	}
	if((multiplicity>=21) && (multiplicity<=24)){
	  hInvMassProtonXi_evtpool6 ->Fill(invMass);
	}
	if((multiplicity>=25) && (multiplicity<=28)){
	  hInvMassProtonXi_evtpool7 ->Fill(invMass);
	}
	if((multiplicity>=29) && (multiplicity<=32)){
	  hInvMassProtonXi_evtpool8 ->Fill(invMass);
	}
	if((multiplicity>=33) && (multiplicity<=36)){
	  hInvMassProtonXi_evtpool9 ->Fill(invMass);
	}
	if((multiplicity>=37) && (multiplicity<=40)){
	  hInvMassProtonXi_evtpool10 ->Fill(invMass);
	}
	if((multiplicity>40) && (multiplicity<=50)){
	  hInvMassProtonXi_evtpool11 ->Fill(invMass);
	}
	if((multiplicity>50) && (multiplicity<=60)){
	  hInvMassProtonXi_evtpool12 ->Fill(invMass);
	}
	if((multiplicity>60) && (multiplicity<=70)){
	  hInvMassProtonXi_evtpool13 ->Fill(invMass);
	}
	if((multiplicity>70) && (multiplicity<=80)){
	  hInvMassProtonXi_evtpool14 ->Fill(invMass);
	}
	if((multiplicity>80) && (multiplicity<=90)){
	  hInvMassProtonXi_evtpool15 ->Fill(invMass);
	}
	if((multiplicity>90) && (multiplicity<=100)){
	  hInvMassProtonXi_evtpool16 ->Fill(invMass);
	}
	if(multiplicity>100){
	  hInvMassProtonXi_evtpool17 ->Fill(invMass);
	}
	
      }
    }   
  }          
  //}

  // pbar + xi+
  // protonbar
  //if(isantiprotonxi){
  for(Int_t itrk=0; itrk < fProtonArrayb->GetEntriesFast(); itrk++){
    TObject *objTrack=fProtonArrayb->At(itrk);    
    if(!objTrack) continue;
    AliAODTrack *proton=dynamic_cast<AliAODTrack*>(objTrack);    
    prpx        =proton->Px();
    prpy        =proton->Py();
    prpz        =proton->Pz();
    energypr    =sqrt(pow(0.938,2)+pow(prpx,2)+pow(prpy,2)+pow(prpz,2));
    
    //mixed event loop starts
    for(Int_t imevt=0; imevt < xippool->GetCurrentNEvents(); imevt++){
      TObjArray *poolcascade=(TObjArray*)xippool->GetEvent(imevt);
      for(Int_t imcascade=0; imcascade < poolcascade->GetEntriesFast(); imcascade++){
	TObject *objCascade=poolcascade->At(imcascade);
	if(!objCascade) continue;
	AliAODcascade *xi=dynamic_cast<AliAODcascade*>(objCascade); 
	xipx        =xi->MomXiX();
	xipy        =xi->MomXiY();
	xipz        =xi->MomXiZ();
	energyxi    =sqrt(pow(1.322,2)+pow(xipx,2)+pow(xipy,2)+pow(xipz,2));	
	
	pt        =sqrt((prpx+xipx)*(prpx+xipx)+(prpy+xipy)*(prpy+xipy));
	energysum = energypr+energyxi;
	invMass   =InvariantMass(xipx,xipy,xipz,prpx,prpy,prpz,energysum);
	hInvMassAntiProtonXi_evtpool->Fill(pt,invMass);       
	relmom    =sqrt(pow((xipx-prpx),2)+pow((xipy-prpy),2)+pow((xipz-prpz),2))/2;
	
	//multiplicity 
	if((multiplicity>=1) && (multiplicity<=4)){
	  hInvMassAntiProtonXi_evtpool1 ->Fill(invMass);
	}
	if((multiplicity>=5) && (multiplicity<=8)){
	  hInvMassAntiProtonXi_evtpool2 ->Fill(invMass);
	}
	if((multiplicity>=9) && (multiplicity<=12)){
	  hInvMassAntiProtonXi_evtpool3 ->Fill(invMass);
	}
	if((multiplicity>=13) && (multiplicity<=16)){
	  hInvMassAntiProtonXi_evtpool4 ->Fill(invMass);
	}
	if((multiplicity>=17) && (multiplicity<=20)){
	  hInvMassAntiProtonXi_evtpool5 ->Fill(invMass);
	}
	if((multiplicity>=21) && (multiplicity<=24)){
	  hInvMassAntiProtonXi_evtpool6 ->Fill(invMass);
	}
	if((multiplicity>=25) && (multiplicity<=28)){
	  hInvMassAntiProtonXi_evtpool7 ->Fill(invMass);
	}
	if((multiplicity>=29) && (multiplicity<=32)){
	  hInvMassAntiProtonXi_evtpool8 ->Fill(invMass);
	}
	if((multiplicity>=33) && (multiplicity<=36)){
	  hInvMassAntiProtonXi_evtpool9 ->Fill(invMass);
	}
	if((multiplicity>=37) && (multiplicity<=40)){
	  hInvMassAntiProtonXi_evtpool10 ->Fill(invMass);
	}
	if((multiplicity>40) && (multiplicity<=50)){
	  hInvMassAntiProtonXi_evtpool11 ->Fill(invMass);
	}
	if((multiplicity>50) && (multiplicity<=60)){
	  hInvMassAntiProtonXi_evtpool12 ->Fill(invMass);
	}
	if((multiplicity>60) && (multiplicity<=70)){
	  hInvMassAntiProtonXi_evtpool13 ->Fill(invMass);
	}
	if((multiplicity>70) && (multiplicity<=80)){
	  hInvMassAntiProtonXi_evtpool14 ->Fill(invMass);
	}
	if((multiplicity>80) && (multiplicity<=90)){
	  hInvMassAntiProtonXi_evtpool15 ->Fill(invMass);
	}
	if((multiplicity>90) && (multiplicity<=100)){
	  hInvMassAntiProtonXi_evtpool16 ->Fill(invMass);
	}
	if(multiplicity>100){
	  hInvMassAntiProtonXi_evtpool17 ->Fill(invMass);
	}
	
      }
    }    
  } 
  // Xi+
  for(Int_t icascade=0; icascade < fCascadeArrayp->GetEntriesFast(); icascade++){
    TObject *objCascade=fCascadeArrayp->At(icascade);    
    if(!objCascade) continue;
    AliAODcascade *xi=dynamic_cast<AliAODcascade*>(objCascade); 
    xipx        =xi->MomXiX();
    xipy        =xi->MomXiY();
    xipz        =xi->MomXiZ();
    energyxi    =sqrt(pow(1.322,2)+pow(xipx,2)+pow(xipy,2)+pow(xipz,2));
    
    //mixed event loop starts
    for(Int_t imevt=0; imevt < antiprotonpool->GetCurrentNEvents(); imevt++){
      TObjArray *pooltrk=(TObjArray*)antiprotonpool->GetEvent(imevt);      
      for(Int_t imtrk=0; imtrk < pooltrk->GetEntriesFast(); imtrk++){
	TObject *objProton=pooltrk->At(imtrk);
	if(!objProton) continue;
	AliAODTrack *proton=dynamic_cast<AliAODTrack*>(objProton); 	
	prpx        =proton->Px();
	prpy        =proton->Py();
	prpz        =proton->Pz();
	energypr    =sqrt(pow(0.938,2)+pow(prpx,2)+pow(prpy,2)+pow(prpz,2));
	
	energysum   =energyxi+energypr;
	pt        =sqrt((prpx+xipx)*(prpx+xipx)+(prpy+xipy)*(prpy+xipy));
	invMass   =InvariantMass(xipx,xipy,xipz,prpx,prpy,prpz,energysum);
	relmom    =sqrt(pow((xipx-prpx),2)+pow((xipy-prpy),2)+pow((xipz-prpz),2))/2;
	hInvMassAntiProtonXi_evtpool->Fill(pt,invMass);       
	
	//multiplicity 
	if((multiplicity>=1) && (multiplicity<=4)){
	  hInvMassAntiProtonXi_evtpool1 ->Fill(invMass);
	}
	if((multiplicity>=5) && (multiplicity<=8)){
	  hInvMassAntiProtonXi_evtpool2 ->Fill(invMass);
	}
	if((multiplicity>=9) && (multiplicity<=12)){
	  hInvMassAntiProtonXi_evtpool3 ->Fill(invMass);
	}
	if((multiplicity>=13) && (multiplicity<=16)){
	  hInvMassAntiProtonXi_evtpool4 ->Fill(invMass);
	}
	if((multiplicity>=17) && (multiplicity<=20)){
	  hInvMassAntiProtonXi_evtpool5 ->Fill(invMass);
	}
	if((multiplicity>=21) && (multiplicity<=24)){
	  hInvMassAntiProtonXi_evtpool6 ->Fill(invMass);
	}
	if((multiplicity>=25) && (multiplicity<=28)){
	  hInvMassAntiProtonXi_evtpool7 ->Fill(invMass);
	}
	if((multiplicity>=29) && (multiplicity<=32)){
	  hInvMassAntiProtonXi_evtpool8 ->Fill(invMass);
	}
	if((multiplicity>=33) && (multiplicity<=36)){
	  hInvMassAntiProtonXi_evtpool9 ->Fill(invMass);
	}
	if((multiplicity>=37) && (multiplicity<=40)){
	  hInvMassAntiProtonXi_evtpool10 ->Fill(invMass);
	}
	if((multiplicity>40) && (multiplicity<=50)){
	  hInvMassAntiProtonXi_evtpool11 ->Fill(invMass);
	}
	if((multiplicity>50) && (multiplicity<=60)){
	  hInvMassAntiProtonXi_evtpool12 ->Fill(invMass);
	}
	if((multiplicity>60) && (multiplicity<=70)){
	  hInvMassAntiProtonXi_evtpool13 ->Fill(invMass);
	}
	if((multiplicity>70) && (multiplicity<=80)){
	  hInvMassAntiProtonXi_evtpool14 ->Fill(invMass);
	}
	if((multiplicity>80) && (multiplicity<=90)){
	  hInvMassAntiProtonXi_evtpool15 ->Fill(invMass);
	}
	if((multiplicity>90) && (multiplicity<=100)){
	  hInvMassAntiProtonXi_evtpool16 ->Fill(invMass);
	}
	if(multiplicity>100){
	  hInvMassAntiProtonXi_evtpool17 ->Fill(invMass);
	}
	
      }
    }   
  }    
  //}

  // lambda + Xi-
  // lambda
  //if(islambdaxi){
  for(Int_t iv0=0; iv0 < fV0Arrayn->GetEntriesFast(); iv0++){
    TObject *objV0=fV0Arrayn->At(iv0);    
    if(!objV0) continue;
    AliAODv0 *lambda=dynamic_cast<AliAODv0*>(objV0); 
    lampx       =lambda->MomV0X();
    lampy       =lambda->MomV0Y();
    lampz       =lambda->MomV0Z();
    energylam   =sqrt(pow(1.115,2)+pow(lampx,2)+pow(lampy,2)+pow(lampz,2));  
    
    //mixed event loop starts
    for(Int_t imevt=0; imevt < xpool->GetCurrentNEvents(); imevt++){
      TObjArray *poolcascade=(TObjArray*)xpool->GetEvent(imevt);
      for(Int_t imcascade=0; imcascade < poolcascade->GetEntriesFast(); imcascade++){
	TObject *objCascade=poolcascade->At(imcascade);
	if(!objCascade) continue;
	AliAODcascade *xi=dynamic_cast<AliAODcascade*>(objCascade); 
	xipx        =xi->MomXiX();
	xipy        =xi->MomXiY();
	xipz        =xi->MomXiZ();
	energyxi    =sqrt(pow(1.322,2)+pow(xipx,2)+pow(xipy,2)+pow(xipz,2));	
	
	energysum   =energyxi+energylam;
	pt          =sqrt(pow((lampx+xipx),2)+pow((lampy+xipy),2));
	invMass     =InvariantMass(xipx,xipy,xipz,lampx,lampy,lampz,energysum);	
	relmom    =sqrt(pow((xipx-lampx),2)+pow((xipy-lampy),2)+pow((xipz-lampz),2))/2;
	hInvMassLambdaXi_evtpool->Fill(pt,invMass);       
	
	//multiplicity //else ifで
	if((multiplicity>=1) && (multiplicity<=5)){
	  hInvMassLambdaXi_evtpool1 ->Fill(invMass);
	}
	if((multiplicity>=6) && (multiplicity<=10)){
	  hInvMassLambdaXi_evtpool2 ->Fill(invMass);
	}
	if((multiplicity>=11) && (multiplicity<=15)){
	  hInvMassLambdaXi_evtpool3 ->Fill(invMass);
	}
	if((multiplicity>=16) && (multiplicity<=20)){
	  hInvMassLambdaXi_evtpool4 ->Fill(invMass);
	}
	if((multiplicity>=21) && (multiplicity<=25)){
	  hInvMassLambdaXi_evtpool5 ->Fill(invMass);
	}
	if((multiplicity>=26) && (multiplicity<=30)){
	  hInvMassLambdaXi_evtpool6 ->Fill(invMass);
	}
	if((multiplicity>=31) && (multiplicity<=35)){
	  hInvMassLambdaXi_evtpool7 ->Fill(invMass);
	}
	if((multiplicity>=36) && (multiplicity<=40)){
	  hInvMassLambdaXi_evtpool8 ->Fill(invMass);
	}
	if((multiplicity>40) && (multiplicity<=50)){
	  hInvMassLambdaXi_evtpool9 ->Fill(invMass);
	}
	if((multiplicity>50) && (multiplicity<=60)){
	  hInvMassLambdaXi_evtpool10 ->Fill(invMass);
	}
	if((multiplicity>60) && (multiplicity<=70)){
	  hInvMassLambdaXi_evtpool11 ->Fill(invMass);
	}
	if((multiplicity>70) && (multiplicity<=80)){
	  hInvMassLambdaXi_evtpool12 ->Fill(invMass);
	}
	if((multiplicity>80) && (multiplicity<=90)){
	  hInvMassLambdaXi_evtpool13 ->Fill(invMass);
	}
	if((multiplicity>90) && (multiplicity<=100)){
	  hInvMassLambdaXi_evtpool14 ->Fill(invMass);
	}
	if(multiplicity>100){
	  hInvMassLambdaXi_evtpool15 ->Fill(invMass);
	}	  
	
      }
    }    
  } 
  // Xi-
  for(Int_t icascade=0; icascade < fCascadeArraym->GetEntriesFast(); icascade++){
    TObject *objCascade=fCascadeArraym->At(icascade);    
    if(!objCascade) continue;
    AliAODcascade *xi=dynamic_cast<AliAODcascade*>(objCascade); 
    xipx        =xi->MomXiX();
    xipy        =xi->MomXiY();
    xipz        =xi->MomXiZ();
    energyxi    =sqrt(pow(1.322,2)+pow(xipx,2)+pow(xipy,2)+pow(xipz,2));
    
    //mixed event loop starts
    for(Int_t imevt=0; imevt < lampool->GetCurrentNEvents(); imevt++){
      TObjArray *poolv0=(TObjArray*)lampool->GetEvent(imevt);
      for(Int_t imv0=0; imv0<poolv0->GetEntriesFast(); imv0++){
	TObject *objpV0=poolv0->At(imv0);
	if(!objpV0) continue;
	AliAODv0 *lambda=dynamic_cast<AliAODv0*>(objpV0);
	lampx       =lambda->MomV0X();
	lampy       =lambda->MomV0Y();
	lampz       =lambda->MomV0Z();
	energylam   =sqrt(pow(1.115,2)+pow(lampx,2)+pow(lampy,2)+pow(lampz,2));  
	
	energysum   =energyxi+energylam;
	pt          =sqrt(pow((lampx+xipx),2)+pow((lampy+xipy),2));
	invMass     =InvariantMass(xipx,xipy,xipz,lampx,lampy,lampz,energysum);	
	relmom    =sqrt(pow((xipx-lampx),2)+pow((xipy-lampy),2)+pow((xipz-lampz),2))/2;
	hInvMassLambdaXi_evtpool->Fill(pt,invMass);       
	
	//multiplicity 
	if((multiplicity>=1) && (multiplicity<=5)){
	  hInvMassLambdaXi_evtpool1 ->Fill(invMass);
	}
	if((multiplicity>=6) && (multiplicity<=10)){
	  hInvMassLambdaXi_evtpool2 ->Fill(invMass);
	}
	if((multiplicity>=11) && (multiplicity<=15)){
	  hInvMassLambdaXi_evtpool3 ->Fill(invMass);
	}
	if((multiplicity>=16) && (multiplicity<=20)){
	  hInvMassLambdaXi_evtpool4 ->Fill(invMass);
	}
	if((multiplicity>=21) && (multiplicity<=25)){
	  hInvMassLambdaXi_evtpool5 ->Fill(invMass);
	}
	if((multiplicity>=26) && (multiplicity<=30)){
	  hInvMassLambdaXi_evtpool6 ->Fill(invMass);
	}
	if((multiplicity>=31) && (multiplicity<=35)){
	  hInvMassLambdaXi_evtpool7 ->Fill(invMass);
	}
	if((multiplicity>=36) && (multiplicity<=40)){
	  hInvMassLambdaXi_evtpool8 ->Fill(invMass);
	}
	if((multiplicity>40) && (multiplicity<=50)){
	  hInvMassLambdaXi_evtpool9 ->Fill(invMass);
	}
	if((multiplicity>50) && (multiplicity<=60)){
	  hInvMassLambdaXi_evtpool10 ->Fill(invMass);
	}
	if((multiplicity>60) && (multiplicity<=70)){
	  hInvMassLambdaXi_evtpool11 ->Fill(invMass);
	}
	if((multiplicity>70) && (multiplicity<=80)){
	  hInvMassLambdaXi_evtpool12 ->Fill(invMass);
	}
	if((multiplicity>80) && (multiplicity<=90)){
          hInvMassLambdaXi_evtpool13 ->Fill(invMass);
	}
	if((multiplicity>90) && (multiplicity<=100)){
	  hInvMassLambdaXi_evtpool14 ->Fill(invMass);
	}
	if(multiplicity>100){
	  hInvMassLambdaXi_evtpool15 ->Fill(invMass);
	}	  
	
      }
    }   
  }    
  //}

  // antilambda + Xi+
  // antilambda
  //if(isantilambdaxi){
  for(Int_t iv0=0; iv0 < fV0Arrayb->GetEntriesFast(); iv0++){
    TObject *objV0=fV0Arrayb->At(iv0);    
    if(!objV0) continue;
    AliAODv0 *lambda=dynamic_cast<AliAODv0*>(objV0); 
    lampx       =lambda->MomV0X();
    lampy       =lambda->MomV0Y();
    lampz       =lambda->MomV0Z();
    energylam   =sqrt(pow(1.115,2)+pow(lampx,2)+pow(lampy,2)+pow(lampz,2));  
    
    //mixed event loop starts
    for(Int_t imevt=0; imevt < xppool->GetCurrentNEvents(); imevt++){
      TObjArray *poolcascade=(TObjArray*)xppool->GetEvent(imevt);
      for(Int_t imcascade=0; imcascade < poolcascade->GetEntriesFast(); imcascade++){
	TObject *objCascade=poolcascade->At(imcascade);
	if(!objCascade) continue;
	AliAODcascade *xi=dynamic_cast<AliAODcascade*>(objCascade); 
	xipx        =xi->MomXiX();
	xipy        =xi->MomXiY();
	xipz        =xi->MomXiZ();
	energyxi    =sqrt(pow(1.322,2)+pow(xipx,2)+pow(xipy,2)+pow(xipz,2));	
	
	energysum   =energyxi+energylam;
	pt          =sqrt(pow((lampx+xipx),2)+pow((lampy+xipy),2));
	invMass     =InvariantMass(xipx,xipy,xipz,lampx,lampy,lampz,energysum);	
	relmom    =sqrt(pow((xipx-lampx),2)+pow((xipy-lampy),2)+pow((xipz-lampz),2))/2;
	hInvMassAntiLambdaXi_evtpool->Fill(pt,invMass);       
	
	//multiplicity 
	if((multiplicity>=1) && (multiplicity<=5)){
	  hInvMassAntiLambdaXi_evtpool1 ->Fill(invMass);
	}
	if((multiplicity>=6) && (multiplicity<=10)){
	  hInvMassAntiLambdaXi_evtpool2 ->Fill(invMass);
	}
	if((multiplicity>=11) && (multiplicity<=15)){
	  hInvMassAntiLambdaXi_evtpool3 ->Fill(invMass);
	}
	if((multiplicity>=16) && (multiplicity<=20)){
	  hInvMassAntiLambdaXi_evtpool4 ->Fill(invMass);
	}
	if((multiplicity>=21) && (multiplicity<=25)){
	  hInvMassAntiLambdaXi_evtpool5 ->Fill(invMass);
	}
	if((multiplicity>=26) && (multiplicity<=30)){
	  hInvMassAntiLambdaXi_evtpool6 ->Fill(invMass);
	}
	if((multiplicity>=31) && (multiplicity<=35)){
	  hInvMassAntiLambdaXi_evtpool7 ->Fill(invMass);
	}
	if((multiplicity>=36) && (multiplicity<=40)){
	  hInvMassAntiLambdaXi_evtpool8 ->Fill(invMass);
	}
	if((multiplicity>40) && (multiplicity<=50)){
	  hInvMassAntiLambdaXi_evtpool9 ->Fill(invMass);
	}
	if((multiplicity>50) && (multiplicity<=60)){
	  hInvMassAntiLambdaXi_evtpool10 ->Fill(invMass);
	}
	if((multiplicity>60) && (multiplicity<=70)){
	  hInvMassAntiLambdaXi_evtpool11 ->Fill(invMass);
	}
	if((multiplicity>70) && (multiplicity<=80)){
	  hInvMassAntiLambdaXi_evtpool12 ->Fill(invMass);
	}
	if((multiplicity>80) && (multiplicity<=90)){
          hInvMassAntiLambdaXi_evtpool13 ->Fill(invMass);
	}
	if((multiplicity>90) && (multiplicity<=100)){
	  hInvMassAntiLambdaXi_evtpool14 ->Fill(invMass);
	}
	if(multiplicity>100){
	  hInvMassAntiLambdaXi_evtpool15 ->Fill(invMass);
	}	  
	
      }
    }    
  } 
  // Xi+
  for(Int_t icascade=0; icascade < fCascadeArrayp->GetEntriesFast(); icascade++){
    TObject *objCascade=fCascadeArrayp->At(icascade);    
    if(!objCascade) continue;
    AliAODcascade *xi=dynamic_cast<AliAODcascade*>(objCascade); 
    xipx        =xi->MomXiX();
    xipy        =xi->MomXiY();
    xipz        =xi->MomXiZ();
    energyxi    =sqrt(pow(1.322,2)+pow(xipx,2)+pow(xipy,2)+pow(xipz,2));
    
    //mixed event loop starts
    for(Int_t imevt=0; imevt < antilampool->GetCurrentNEvents(); imevt++){
      TObjArray *poolv0=(TObjArray*)antilampool->GetEvent(imevt);
      for(Int_t imv0=0; imv0<poolv0->GetEntriesFast(); imv0++){
	TObject *objpV0=poolv0->At(imv0);
	if(!objpV0) continue;
	AliAODv0 *lambda=dynamic_cast<AliAODv0*>(objpV0);
	lampx       =lambda->MomV0X();
	lampy       =lambda->MomV0Y();
	lampz       =lambda->MomV0Z();
	energylam   =sqrt(pow(1.115,2)+pow(lampx,2)+pow(lampy,2)+pow(lampz,2));  
	
	energysum   =energyxi+energylam;
	pt          =sqrt(pow((lampx+xipx),2)+pow((lampy+xipy),2));
	invMass     =InvariantMass(xipx,xipy,xipz,lampx,lampy,lampz,energysum);	
	relmom    =sqrt(pow((xipx-lampx),2)+pow((xipy-lampy),2)+pow((xipz-lampz),2))/2;
	hInvMassAntiLambdaXi_evtpool->Fill(pt,invMass);       
	
	//multiplicity 
	if((multiplicity>=1) && (multiplicity<=5)){
	  hInvMassAntiLambdaXi_evtpool1 ->Fill(invMass);
	}
	if((multiplicity>=6) && (multiplicity<=10)){
	  hInvMassAntiLambdaXi_evtpool2 ->Fill(invMass);
	}
	if((multiplicity>=11) && (multiplicity<=15)){
	  hInvMassAntiLambdaXi_evtpool3 ->Fill(invMass);
	}
	if((multiplicity>=16) && (multiplicity<=20)){
	    hInvMassAntiLambdaXi_evtpool4 ->Fill(invMass);
	}
	if((multiplicity>=21) && (multiplicity<=25)){
	  hInvMassAntiLambdaXi_evtpool5 ->Fill(invMass);
	}
	if((multiplicity>=26) && (multiplicity<=30)){
	  hInvMassAntiLambdaXi_evtpool6 ->Fill(invMass);
	}
	if((multiplicity>=31) && (multiplicity<=35)){
	  hInvMassAntiLambdaXi_evtpool7 ->Fill(invMass);
	}
	if((multiplicity>=36) && (multiplicity<=40)){
	  hInvMassAntiLambdaXi_evtpool8 ->Fill(invMass);
	}
	if((multiplicity>40) && (multiplicity<=50)){
	  hInvMassAntiLambdaXi_evtpool9 ->Fill(invMass);
	}
	if((multiplicity>50) && (multiplicity<=60)){
	  hInvMassAntiLambdaXi_evtpool10 ->Fill(invMass);
	}
	if((multiplicity>60) && (multiplicity<=70)){
	  hInvMassAntiLambdaXi_evtpool11 ->Fill(invMass);
	}
	if((multiplicity>70) && (multiplicity<=80)){
	  hInvMassAntiLambdaXi_evtpool12 ->Fill(invMass);
	}
	if((multiplicity>80) && (multiplicity<=90)){
          hInvMassAntiLambdaXi_evtpool13 ->Fill(invMass);
	}
	if((multiplicity>90) && (multiplicity<=100)){
	  hInvMassAntiLambdaXi_evtpool14 ->Fill(invMass);
	}
	if(multiplicity>100){
	  hInvMassAntiLambdaXi_evtpool15 ->Fill(invMass);
	}	  
	
      }
    }   
  }
  //}    

  //lambda
  //if(islambdalambda){
  TObjArray* tracksClonelambda=(TObjArray*)fV0Arrayn->Clone();
  tracksClonelambda->SetOwner();
  lambdapool->UpdatePool(tracksClonelambda);
  //}
  //lambda bar
  //if(isantilambdaantilambda){
  TObjArray* tracksClonelambdab=(TObjArray*)fV0Arrayb->Clone();
  tracksClonelambdab->SetOwner();
  antilambdapool->UpdatePool(tracksClonelambdab);    
  //}

  //Xi-
  //if(isprotonxi){
  TObjArray* tracksClonexim=(TObjArray*)fCascadeArraym->Clone();
  tracksClonexim->SetOwner();
  xipool->UpdatePool(tracksClonexim);
  //proton
  TObjArray* tracksClonep=(TObjArray*)fProtonArray->Clone();
  tracksClonep->SetOwner();
  protonpool->UpdatePool(tracksClonep);
  //}
  //Xi+
  //if(isantiprotonxi){
  TObjArray* tracksClonexip=(TObjArray*)fCascadeArrayp->Clone();
  tracksClonexip->SetOwner();
  xippool->UpdatePool(tracksClonexip);
  //proton bar
  TObjArray* tracksClonepb=(TObjArray*)fProtonArrayb->Clone();
  tracksClonepb->SetOwner();
  antiprotonpool->UpdatePool(tracksClonepb);
  //}

  //lambda
  //if(islambdaxi){
  TObjArray* tracksClonelam=(TObjArray*)fV0Arrayn->Clone();
  tracksClonelam->SetOwner();
  lampool->UpdatePool(tracksClonelam);
  //Xi-
  TObjArray* tracksClonex=(TObjArray*)fCascadeArraym->Clone();
  tracksClonex->SetOwner();
  xpool->UpdatePool(tracksClonex);
  //}
  //lambda bar
  //if(isantilambdaxi){
  TObjArray* tracksClonelamb=(TObjArray*)fV0Arrayb->Clone();
  tracksClonelamb->SetOwner();
  antilampool->UpdatePool(tracksClonelamb);    
  TObjArray* tracksClonexp=(TObjArray*)fCascadeArrayp->Clone();
  tracksClonexp->SetOwner();
  xppool->UpdatePool(tracksClonexp); 
  //}

  // multiplicity
  if(nLamn > 0){
    fmultiplicity_lambda->Fill(multiplicity);
  }
  if(nLamb > 0){
    fmultiplicity_antilambda->Fill(multiplicity);
  }
  if(np > 0){
    fmultiplicity_proton->Fill(multiplicity);
    }
  if(npb > 0){
    fmultiplicity_antiproton->Fill(multiplicity);
    }
  if(nXim > 0){
    fmultiplicity_xim->Fill(multiplicity);
    }
  if(nXip > 0){
    fmultiplicity_xip->Fill(multiplicity);
    }

  // Post output data
  PostData(1,fOutputList);

  /*
    cout<<"+++++++++++++++++++++"<<endl;
    cout<<"++ Analysis Finish ++"<<endl;
    cout<<"++     bye bye..   ++"<<endl;
    cout<<"event number = "<<nevt<<endl;
    cout<<"+++++++++++++++++++++"<<endl; 
  */
  
}


//________________________________________________________________________
Bool_t AliAnalysisTaskPPiLambda::EventSelection(AliAODEvent* data)
{
  const AliVVertex *vtx   =data->GetPrimaryVertex();
  const AliVVertex *vtxSPD=data->GetPrimaryVertexSPD();  
  
  Double_t xvtx=0.,yvtx=0.,zvtx=0.;
  Double_t xvtxSPD=0.,yvtxSPD=0.,zvtxSPD=0.;
  Int_t    ncont=0,ncontSPD=0;
  Double_t vdisp=0.;
  // Event vertex
  xvtx =vtx->GetX();
  yvtx =vtx->GetY();
  zvtx =vtx->GetZ();
  ncont=vtx->GetNContributors();

  fEvtVtxX   ->Fill(xvtx);
  fEvtVtxY   ->Fill(yvtx);
  fEvtVtxZ   ->Fill(zvtx);
  fEvtVtxTrk ->Fill(ncont);
  // SPD vertex
  xvtxSPD =vtxSPD->GetX();
  yvtxSPD =vtxSPD->GetY();
  zvtxSPD =vtxSPD->GetZ();
  ncontSPD=vtxSPD->GetNContributors();
  vdisp=fabs(zvtx-zvtxSPD);

  fSPDVtxZ   ->Fill(zvtxSPD);
  fSPDVtxTrk ->Fill(ncontSPD);
  fSPDVtxCor ->Fill(zvtx,zvtxSPD);
  fSPDVtxDisp->Fill(vdisp);

  Bool_t zvtx_cut=kFALSE;
  zvtx_cut=(fabs(zvtx)<10. && ncont>1 && ncontSPD>0 && vdisp<0.5);
  return zvtx_cut;
}
//________________________________________________________________________

//invariant mass Lambda cascade class
Double_t AliAnalysisTaskPPiLambda::InvMassLambda(AliAODcascade *casc)
{  
  Double_t pPx=0.,pPy=0.,pPz=0.,pE=0.; //proton
  Double_t nPx=0.,nPy=0.,nPz=0.,nE=0.; //pion-
  Float_t  energysum=0.,psum2=0.;
  Double_t invMass =0.;  
  pPx=casc->MomPosX();
  pPy=casc->MomPosY();
  pPz=casc->MomPosZ();
  nPx=casc->MomNegX();
  nPy=casc->MomNegY();
  nPz=casc->MomNegZ();
  pE =sqrt(0.938*0.938+pPx*pPx+pPy*pPy+pPz*pPz); //proton
  nE =sqrt(0.14*0.14+nPx*nPx+nPy*nPy+nPz*nPz);   //pion-
  energysum =pE+nE;
  psum2     =(pPx+nPx)*(pPx+nPx)+(pPy+nPy)*(pPy+nPy)+(pPz+nPz)*(pPz+nPz);
  invMass   =((energysum*energysum-psum2)>0)?sqrt((energysum*energysum)-psum2):0; 
  return invMass;
}
//invariant mass anti lambda
Double_t AliAnalysisTaskPPiLambda::InvMassAntiLambda(AliAODcascade *casc)
{  
  Double_t nPx=0.,nPy=0.,nPz=0.,nE=0.; //anti proton
  Double_t pPx=0.,pPy=0.,pPz=0.,pE=0.; //pion+
  Float_t  energysum=0.,psum2=0.;
  Double_t invMass =0.;  
  nPx=casc->MomNegX();
  nPy=casc->MomNegY();
  nPz=casc->MomNegZ();
  pPx=casc->MomPosX();
  pPy=casc->MomPosY();
  pPz=casc->MomPosZ();
  nE =sqrt(0.938*0.938+nPx*nPx+nPy*nPy+nPz*nPz); //anti proton
  pE =sqrt(0.14*0.14+pPx*pPx+pPy*pPy+pPz*pPz);   //pion+
  energysum =nE+pE;
  psum2     =(nPx+pPx)*(nPx+pPx)+(nPy+pPy)*(nPy+pPy)+(nPz+pPz)*(nPz+pPz);
  invMass   =((energysum*energysum-psum2)>0)?sqrt((energysum*energysum)-psum2):0; 
  return invMass;
}
//invariant mass Xi-
Double_t AliAnalysisTaskPPiLambda::InvMassXi(AliAODcascade *casc)
{  
  Double_t lPx=0.,lPy=0.,lPz=0.,lE=0.; //lambda
  Double_t bPx=0.,bPy=0.,bPz=0.,bE=0.; //bachelor pion-
  Float_t  energysum=0.,psum2=0.;
  Double_t invMass =0.;  
  lPx=casc->MomV0X();
  lPy=casc->MomV0Y();
  lPz=casc->MomV0Z();
  bPx=casc->MomBachX();
  bPy=casc->MomBachY();
  bPz=casc->MomBachZ();
  lE =sqrt(1.115*1.115+lPx*lPx+lPy*lPy+lPz*lPz); //lambda
  bE =sqrt(0.14*0.14+bPx*bPx+bPy*bPy+bPz*bPz);   //pion-
  energysum =lE+bE;
  psum2     =(lPx+bPx)*(lPx+bPx)+(lPy+bPy)*(lPy+bPy)+(lPz+bPz)*(lPz+bPz);
  invMass   =((energysum*energysum-psum2)>0)?sqrt((energysum*energysum)-psum2):0; 
  return invMass;
}
//omega rejection
Double_t AliAnalysisTaskPPiLambda::InvMassOmega(AliAODcascade *casc)
{  
  Double_t lPx=0.,lPy=0.,lPz=0.,lE=0.; //lambda
  Double_t bPx=0.,bPy=0.,bPz=0.,bE=0.; //bachelor kaon
  Float_t  energysum=0.,psum2=0.;
  Double_t invMass =0.;  
  lPx=casc->MomV0X();
  lPy=casc->MomV0Y();
  lPz=casc->MomV0Z();
  bPx=casc->MomBachX();
  bPy=casc->MomBachY();
  bPz=casc->MomBachZ();
  lE =sqrt(1.116*1.116+lPx*lPx+lPy*lPy+lPz*lPz); //lambda
  bE =sqrt(0.494*0.494+bPx*bPx+bPy*bPy+bPz*bPz); //kaon-
  energysum =lE+bE;
  psum2     =(lPx+bPx)*(lPx+bPx)+(lPy+bPy)*(lPy+bPy)+(lPz+bPz)*(lPz+bPz);
  invMass   =((energysum*energysum-psum2)>0)?sqrt((energysum*energysum)-psum2):0; 
  return invMass;
}
//lambda cosine pointing angle
Double_t AliAnalysisTaskPPiLambda::LambdaCosPointingAngle(AliAODcascade *casc,const Double_t *DecayVtx,
						    const Float_t *point) const 
{  
  TVector3 v0Mom(casc->MomV0X(),casc->MomV0Y(),casc->MomV0Z());
  TVector3 fline(DecayVtx[0] - point[0], DecayVtx[1] - point[1],
                 DecayVtx[2] - point[2]);
    Double_t ptot2 = v0Mom.Mag2() * fline.Mag2();
  if (ptot2 <= 0) {
    return 0.0;
  }
  else {
    Double_t cos = v0Mom.Dot(fline) / TMath::Sqrt(ptot2);
    if (cos > 1.0)
      cos = 1.0;
    if (cos < -1.0)
      cos = -1.0;
    return cos;
  }
}
//transverse radius of the lambda decay vertex
Double_t AliAnalysisTaskPPiLambda::DecayLengthXY(const Double_t *DecayVtx,const Float_t *point) const 
{ return TMath::Sqrt( (point[0] - DecayVtx[0]) * (point[0] - DecayVtx[0])
		     + (point[1] - DecayVtx[1]) * (point[1] - DecayVtx[1]));
}
//transverse radius of the xi decay vertex
Double_t AliAnalysisTaskPPiLambda::xiDecayLengthXY(const Double_t *xiDecayVtx,const Float_t *point) const 
{ return TMath::Sqrt( (point[0] - xiDecayVtx[0]) * (point[0] - xiDecayVtx[0])
		     + (point[1] - xiDecayVtx[1]) * (point[1] - xiDecayVtx[1]));
}
//invariant mass Lambda v0 class
Double_t AliAnalysisTaskPPiLambda::InvMasslambda(AliAODv0 *v0)
{  
  Double_t pPx=0.,pPy=0.,pPz=0.,pE=0.; //proton
  Double_t nPx=0.,nPy=0.,nPz=0.,nE=0.; //pion
  Float_t  energysum=0.,psum2=0.;
  Double_t invMass =0.;
  pPx=v0->MomPosX();
  pPy=v0->MomPosY();
  pPz=v0->MomPosZ();
  nPx=v0->MomNegX();
  nPy=v0->MomNegY();
  nPz=v0->MomNegZ();
  pE =sqrt(0.938*0.938+pPx*pPx+pPy*pPy+pPz*pPz); //proton
  nE =sqrt(0.14*0.14+nPx*nPx+nPy*nPy+nPz*nPz);   //pion-
  energysum =pE+nE;
  psum2     =(pPx+nPx)*(pPx+nPx)+(pPy+nPy)*(pPy+nPy)+(pPz+nPz)*(pPz+nPz);
  invMass   =((energysum*energysum-psum2)>0)?sqrt((energysum*energysum)-psum2):0; 
  return invMass;
}
//K0 rejection                                                                                         
Double_t AliAnalysisTaskPPiLambda::InvMassK0(AliAODv0 *v0)
{
  Double_t pPx=0.,pPy=0.,pPz=0.,pE=0.; //pion+                                                         
  Double_t nPx=0.,nPy=0.,nPz=0.,nE=0.; //pion                                                          
  Float_t  energysum=0.,psum2=0.;
  Double_t invMass =0.;
  pPx=v0->MomPosX();
  pPy=v0->MomPosY();
  pPz=v0->MomPosZ();
  nPx=v0->MomNegX();
  nPy=v0->MomNegY();
  nPz=v0->MomNegZ();
  pE =sqrt(0.14*0.14+pPx*pPx+pPy*pPy+pPz*pPz); //pion+                                                 
  nE =sqrt(0.14*0.14+nPx*nPx+nPy*nPy+nPz*nPz); //pion-                                                 
  energysum =pE+nE;
  psum2     =(pPx+nPx)*(pPx+nPx)+(pPy+nPy)*(pPy+nPy)+(pPz+nPz)*(pPz+nPz);
  invMass   =((energysum*energysum-psum2)>0)?sqrt((energysum*energysum)-psum2):0;
  return invMass;
}
//invariant mass antiLambda v0 class
Double_t AliAnalysisTaskPPiLambda::InvMassAntilambda(AliAODv0 *v0)
{
  Double_t pPx=0.,pPy=0.,pPz=0.,pE=0.; //pion+                                                                                                            
  Double_t nPx=0.,nPy=0.,nPz=0.,nE=0.; //antiproton                                                                                                      
  Float_t  energysum=0.,psum2=0.;
  Double_t invMass =0.;
  pPx=v0->MomPosX();
  pPy=v0->MomPosY();
  pPz=v0->MomPosZ();
  nPx=v0->MomNegX();
  nPy=v0->MomNegY();
  nPz=v0->MomNegZ();
  pE =sqrt(0.14*0.14+pPx*pPx+pPy*pPy+pPz*pPz);     //pion+                                                                                                  
  nE =sqrt(0.938*0.938+nPx*nPx+nPy*nPy+nPz*nPz);   //antiproton                                                                                             
  energysum =pE+nE;
  psum2     =(pPx+nPx)*(pPx+nPx)+(pPy+nPy)*(pPy+nPy)+(pPz+nPz)*(pPz+nPz);
  invMass   =((energysum*energysum-psum2)>0)?sqrt((energysum*energysum)-psum2):0;
  return invMass;
}
Double_t AliAnalysisTaskPPiLambda::CosPointingAngle(AliAODv0 *v0,const Double_t *DecayVtx,
						    const Float_t *point) const {  
  /// Cosine of pointing angle in space assuming it is produced at "point"
  TVector3 v0Mom(v0->MomV0X(),v0->MomV0Y(),v0->MomV0Z());
  TVector3 fline(DecayVtx[0] - point[0], DecayVtx[1] - point[1],
                 DecayVtx[2] - point[2]);
    Double_t ptot2 = v0Mom.Mag2() * fline.Mag2();
  if (ptot2 <= 0) {
    return 0.0;
  }
  else {
    Double_t cos = v0Mom.Dot(fline) / TMath::Sqrt(ptot2);
    if (cos > 1.0)
      cos = 1.0;
    if (cos < -1.0)
      cos = -1.0;
    return cos;
  }
}

Double_t AliAnalysisTaskPPiLambda::OpenAngle(Double_t px1,Double_t py1,Double_t pz1,
					     Double_t px2,Double_t py2,Double_t pz2){
  Double_t lScalPtot1Ptot2=0.,lPtot1xPtot2=0.;

  lScalPtot1Ptot2 = (px1*px2)+(py1*py2)+(pz1*pz2);
  lPtot1xPtot2 = (sqrt(pow(px1,2)+pow(py1,2)+pow(pz1,2)))*(sqrt(pow(px2,2)+pow(py2,2)+pow(pz2,2)));
  
  return acos(lScalPtot1Ptot2/lPtot1xPtot2);
}

Double_t AliAnalysisTaskPPiLambda::InvariantMass(Double_t px1,Double_t py1,Double_t pz1,
						 Double_t px2,Double_t py2,Double_t pz2,Double_t energysum){

  Double_t psum2=0.,pt=0.,invMass=0.;

  psum2     =pow((px1+px2),2)+pow((py1+py2),2)+pow((pz1+pz2),2);
  pt        =sqrt(pow(px1+px2,2)+pow(py1+py2,2));
  invMass   =sqrt((energysum*energysum)-psum2);

  return invMass;

}
