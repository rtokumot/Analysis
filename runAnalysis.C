#include "AliAnalysisTaskPPiLambda.h"
#include "AliAnalysisTaskPIDResponse.h"

void runAnalysis(bool grid=1,int gmode=1){
  // if grid==true/false, run on grid/local
  // gmode: 0=test, 1=full, 2=terminate

  // Load common libraries
  //================================================================
  gSystem->Load("libCore.so");  
  gSystem->Load("libTree.so");
  gSystem->Load("libGeom.so");
  gSystem->Load("libVMC.so");
  gSystem->Load("libPhysics.so");
  gSystem->Load("libSTEERBase");
  gSystem->Load("libESD");
  gSystem->Load("libAOD");
  gSystem->Load("libANALYSIS");
  gSystem->Load("libANALYSISalice");   
  gSystem->Load("libANALYSISaliceBase");
  gSystem->Load("libOADB");

  // Use AliRoot includes
  //================================================================
#if !defined (__CINT__) || defined (__CLING__)
  gInterpreter->ProcessLine(".include $ROOTSYS/include");
  gInterpreter->ProcessLine(".include $ALICE_ROOT/include");
  gInterpreter->ProcessLine(".include $ALICE_PHYSICS/include");
#else
  gROOT->ProcessLine(".include $ROOTSYS/include");
  gROOT->ProcessLine(".include $ALICE_ROOT/include");
  gROOT->ProcessLine(".include $ALICE_PHYSICS/include");
#endif

  // Create the analysis manager
  //================================================================
  AliAnalysisManager *mgr =new AliAnalysisManager("AnalsysTaskExample");
  AliAODInputHandler *aodH=new AliAODInputHandler();
  mgr->SetInputEventHandler(aodH);

  // Load my analysis task
  //================================================================
#if !defined (__CINT__) || defined (__CLING__)
  gInterpreter->LoadMacro("AliAnalysisTaskPPiLambda.cxx++g");  
  AliAnalysisTaskPPiLambda *task=reinterpret_cast<AliAnalysisTaskPPiLambda*>(gInterpreter->ExecuteMacro("AddTaskPPiLambda.C()"));
  AliAnalysisTaskPIDResponse *pidtask=reinterpret_cast<AliAnalysisTaskPIDResponse*>(gInterpreter->ExecuteMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C"));
#else
  gROOT->LoadMacro("AliAnalysisTaskPPiLambda.cxx++g");  
  gROOT->LoadMacro("AddTaskPPiLambda.C");
  AliAnalysisTaskPPiLambda *task=AddTaskPPiLambda();
  gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C");
  AliAnalysisTaskPIDResponse *pidtask=AddTaskPIDResponse();
#endif

  // Debug
  //================================================================
  if(!mgr->InitAnalysis()) return;
  mgr->SetDebugLevel(2);
  mgr->PrintStatus();
  mgr->SetUseProgressBar(1,25);

  // Start analysis
  if(!grid){ // local test
    TChain *chain=new TChain("aodTree");
    // samples/AliAOD.root copied from 
    // /alice/data/2016/LHC16q/000265525/pass2_CENT_woSDD/16000265525019.100/
    //chain->Add("AliAOD.root");
    chain->Add("AliAOD.root.merged");
    mgr->StartAnalysis("local",chain);
  }
  else{ // analysis on grid.
    // Create the plugin
    //================================================================
    AliAnalysisAlien *alienHandler=new AliAnalysisAlien();
    
    // Specify the include paths on grid
    //================================================================
    alienHandler->AddIncludePath("-I$ROOTSYS/include");
    alienHandler->AddIncludePath("-I$ALICE_ROOT/include");
    alienHandler->AddIncludePath("-I$ALICE_PHYSICS/include");

    // Copy of source files 
    //================================================================
    alienHandler->SetAnalysisSource("AliAnalysisTaskPPiLambda.cxx");
    alienHandler->SetAdditionalLibs("AliAnalysisTaskPPiLambda.cxx AliAnalysisTaskPPiLambda.h");
    alienHandler->SetAliPhysicsVersion("vAN-20200904_ROOT6-1");

    // Input data
    //================================================================
    alienHandler->SetGridDataDir("/alice/data/2016/LHC16q");
    alienHandler->SetDataPattern("pass2_CENT_woSDD/*/AliAOD.root");
    alienHandler->SetRunPrefix("000");
    // selected runs
    /*
    int runnum=0;
    ifstream frun("runlist.dat");
    while(!frun.eof()){
      runnum=0;
      frun>>runnum;
      alienHandler->AddRunNumber(runnum);
    }
    */
    alienHandler->AddRunNumber(265525); // For test

    // Job setting
    //================================================================
    alienHandler->SetAPIVersion("V1.1x");
    alienHandler->SetNrunsPerMaster(1); // 1 output folder/input run
    alienHandler->SetMergeViaJDL(1); // merged output/input run
    alienHandler->SetSplitMaxInputFileNumber(20); // 20 input files/sub-job
    alienHandler->SetTTL(30000);
    alienHandler->SetJDLName("analysis_ryoka_PPiLambda.jdl");
    alienHandler->SetOutputToRunNo(kTRUE);
    alienHandler->SetKeepLogs(kTRUE);
    alienHandler->SetGridWorkingDir("test/AliAnalysisTaskPPiLambda");
    alienHandler->SetGridOutputDir("1029001");
    alienHandler->SetDefaultOutputs(kTRUE);
    alienHandler->SetCheckCopy(kFALSE);
    if     (gmode==0){ alienHandler->SetNtestFiles(1); alienHandler->SetRunMode("test"); }
    else if(gmode==1){ alienHandler->SetRunMode("full"); }
    else             { alienHandler->SetRunMode("terminate"); }

    // Executables
    //================================================================
    alienHandler->SetExecutableCommand("aliroot -b -q");
    //  alienHandler->SetAnalysisMacro("analysis_ryoka_PPiLambda.C");
    //  alienHandler->SetExecutable("analysis_ryoka_PPiLambda.sh");
    alienHandler->SetInputFormat("xml-single");
    alienHandler->SetSplitMode("se");
    mgr->SetGridHandler(alienHandler);
    mgr->StartAnalysis("grid");
  }


}
