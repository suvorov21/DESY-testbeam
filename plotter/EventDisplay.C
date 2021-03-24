#include <TGClient.h>
#include <TGButton.h>
#include <TGLabel.h>
#include <TGFrame.h>
#include <TGLayout.h>
#include <TGSplitter.h>
#include <TCanvas.h>
#include "TROOT.h"
#include "TApplication.h"
#include "TRootEmbeddedCanvas.h"
#include "TGDNDManager.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include "TTree.h"
#include "TStyle.h"
#include "TBox.h"
#include "TGNumberEntry.h"
#include "TGTextEntry.h"

#include <iostream>

class T2K {
 public:
  TStyle* SetT2KStyle(Int_t WhichStyle, TString styleName) {
    TStyle *t2kStyle= new TStyle(styleName, "T2K approved plots style");

    // -- WhichStyle --
    // 1 = presentation large fonts
    // 2 = presentation small fonts
    // 3 = publication/paper

    Int_t FontStyle = 22;
    Float_t FontSizeLabel = 0.035;
    Float_t FontSizeTitle = 0.05;
    Float_t YOffsetTitle = 1.3;

    switch(WhichStyle) {
    case 1:
      FontStyle = 42;
      FontSizeLabel = 0.05;
      FontSizeTitle = 0.065;
      YOffsetTitle = 1.19;
      break;
    case 2:
      FontStyle = 42;
      FontSizeLabel = 0.035;
      FontSizeTitle = 0.05;
      YOffsetTitle = 1.6;
      break;
    case 3:
      FontStyle = 132;
      FontSizeLabel = 0.035;
      FontSizeTitle = 0.05;
      YOffsetTitle = 1.6;
      break;
    }

    // use plain black on white colors
    t2kStyle->SetFrameBorderMode(0);
    t2kStyle->SetCanvasBorderMode(0);
    t2kStyle->SetPadBorderMode(0);
    t2kStyle->SetCanvasBorderSize(0);
    t2kStyle->SetFrameBorderSize(0);
    t2kStyle->SetDrawBorder(0);
    t2kStyle->SetTitleBorderSize(0);

    t2kStyle->SetPadColor(0);
    t2kStyle->SetCanvasColor(0);
    t2kStyle->SetStatColor(0);
    t2kStyle->SetFillColor(0);

    t2kStyle->SetEndErrorSize(4);
    t2kStyle->SetStripDecimals(kFALSE);

    t2kStyle->SetLegendBorderSize(0);
    t2kStyle->SetLegendFont(FontStyle);

    // set the paper & margin sizes
    t2kStyle->SetPaperSize(20, 26);
    t2kStyle->SetPadTopMargin(0.1);
    t2kStyle->SetPadBottomMargin(0.15);
    t2kStyle->SetPadRightMargin(0.15); // 0.075 -> 0.13 for colz option
    t2kStyle->SetPadLeftMargin(0.16);//to include both large/small font options

    // Fonts, sizes, offsets
    t2kStyle->SetTextFont(FontStyle);
    t2kStyle->SetTextSize(0.08);

    t2kStyle->SetLabelFont(FontStyle, "x");
    t2kStyle->SetLabelFont(FontStyle, "y");
    t2kStyle->SetLabelFont(FontStyle, "z");
    t2kStyle->SetLabelFont(FontStyle, "t");
    t2kStyle->SetLabelSize(FontSizeLabel, "x");
    t2kStyle->SetLabelSize(FontSizeLabel, "y");
    t2kStyle->SetLabelSize(FontSizeLabel, "z");
    t2kStyle->SetLabelOffset(0.015, "x");
    t2kStyle->SetLabelOffset(0.015, "y");
    t2kStyle->SetLabelOffset(0.015, "z");

    t2kStyle->SetTitleFont(FontStyle, "x");
    t2kStyle->SetTitleFont(FontStyle, "y");
    t2kStyle->SetTitleFont(FontStyle, "z");
    t2kStyle->SetTitleFont(FontStyle, "t");
    t2kStyle->SetTitleSize(FontSizeTitle, "y");
    t2kStyle->SetTitleSize(FontSizeTitle, "x");
    t2kStyle->SetTitleSize(FontSizeTitle, "z");
    t2kStyle->SetTitleOffset(1.14, "x");
    t2kStyle->SetTitleOffset(YOffsetTitle, "y");
    t2kStyle->SetTitleOffset(1.2, "z");

    t2kStyle->SetTitleStyle(0);
    t2kStyle->SetTitleFontSize(0.06);//0.08
    t2kStyle->SetTitleFont(FontStyle, "pad");
    t2kStyle->SetTitleBorderSize(0);
    t2kStyle->SetTitleX(0.1f);
    t2kStyle->SetTitleW(0.8f);

    // use bold lines and markers
    t2kStyle->SetMarkerStyle(20);
    t2kStyle->SetHistLineWidth( Width_t(2.5) );
    t2kStyle->SetLineStyleString(2, "[12 12]"); // postscript dashes

    // get rid of X error bars and y error bar caps
    t2kStyle->SetErrorX(0.001);

    // do not display any of the standard histogram decorations
    t2kStyle->SetOptTitle(1);
    t2kStyle->SetOptStat(0);
    t2kStyle->SetOptFit(1); // switch fitter pad on

    // put tick marks on top and RHS of plots
    t2kStyle->SetPadTickX(1);
    t2kStyle->SetPadTickY(1);

    // -- color --
    // functions blue
    t2kStyle->SetFuncColor(600-4);

    t2kStyle->SetFillColor(1); // make color fillings (not white)
    // - color setup for 2D -
    // - "cold"/ blue-ish -
    Double_t red[]   = { 0.00, 0.00, 0.00 };
    Double_t green[] = { 1.00, 0.00, 0.00 };
    Double_t blue[]  = { 1.00, 1.00, 0.25 };
    // - "warm" red-ish colors -
    //  Double_t red[]   = {1.00, 1.00, 0.25 };
    //  Double_t green[] = {1.00, 0.00, 0.00 };
    //  Double_t blue[]  = {0.00, 0.00, 0.00 };

    Double_t stops[] = { 0.25, 0.75, 1.00 };
    const Int_t NRGBs = 3;
    const Int_t NCont = 500;

    TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
    t2kStyle->SetNumberContours(NCont);

    // - Rainbow -
    //  t2kStyle->SetPalette(1);  // use the rainbow color set

    // Setup light blue - bright yellow pallete
    // gStyle->SetPalette(kBird);

    // -- axis --
    t2kStyle->SetStripDecimals(kFALSE); // don't do 1.0 -> 1
    //  TGaxis::SetMaxDigits(3); // doesn't have an effect
    // no supressed zeroes!
    t2kStyle->SetHistMinimumZero(kTRUE);


   return(t2kStyle);
  }


  void CenterHistoTitles(TH1 *thisHisto){
    thisHisto->GetXaxis()->CenterTitle();
    thisHisto->GetYaxis()->CenterTitle();
    thisHisto->GetZaxis()->CenterTitle();
  }


  int AddGridLinesToPad(TPad *thisPad) {
    thisPad->SetGridx();
    thisPad->SetGridy();
    return(0);
  }
};


class MyMainFrame : public TGMainFrame {

public:
  MyMainFrame(const TGWindow *p, UInt_t w, UInt_t h, TString name);
  virtual ~MyMainFrame();
  void     CloseWindow();

  void DoDraw();
  void NextEvent();
  void PrevEvent();
  void UpdateNumber();

  void EventInfo(Int_t event, Int_t px, Int_t py, TObject *selected);

  ClassDef(MyMainFrame, 0)

protected:
  // file
  TFile* f;
  TTree* t;
  bool saclay_cosmics;
  Int_t padAmpl[36][32][511];
  Int_t padAmpl_saclay[36][32][510];
  Int_t eventID = 0;

  TH2F* MM;
  TH1F* WF[9];
  Int_t WFstart = 100;
  Int_t WFend = 260;


  // GUI
  TCanvas *f_ED_canvas;
  TCanvas *f_WF_canvas;
  TGHorizontalFrame *fhf;
  TRootEmbeddedCanvas  *fED;
  TRootEmbeddedCanvas  *fWF;

  TGTextButton* fButtonExit;
  TGTextButton* fButtonDraw;
  TGTextButton* fNextEvent;
  TGTextButton* fPrevEvent;

  TGNumberEntry* fNumber;
  TGTextEntry* fEntry;
};

//______________________________________________________________________________
MyMainFrame::MyMainFrame(const TGWindow *p, UInt_t w, UInt_t h, TString name) :
TGMainFrame(p, w, h)
{
  TString localStyleName = "T2K";
  int localWhichStyle = 3;
  TStyle* t2kstyle = T2K().SetT2KStyle(localWhichStyle, localStyleName);
  gROOT->SetStyle(t2kstyle->GetName());
  gROOT->ForceStyle();
  // read file
  f = new TFile(name, "READ");
  t = (TTree*)f->Get("tree");

  saclay_cosmics = false;
  TString branch_name = t->GetBranch("PadAmpl")->GetTitle();
  if (branch_name.Contains("[510]")) {
    saclay_cosmics = true;
    t->SetBranchAddress("PadAmpl", padAmpl_saclay);
  } else if (branch_name.Contains("[511]")) {
    saclay_cosmics = false;
    t->SetBranchAddress("PadAmpl", padAmpl);
  } else {
    std::cerr << "Time binning is unknown." << std::endl;
    std::cerr << "Read from file " << branch_name  << std::endl;
    exit(1);
  }
  MM = new TH2F("h", "", 38, -1., 37., 34, -1., 33.);
  for (auto i = 0; i < 9; ++i)
    WF[i] = new TH1F(Form("WF_%i", i), "", 511, 0., 511);

  // ini GUI
  TGHorizontalFrame* fMain = new TGHorizontalFrame(this, w, h);
  fED = new TRootEmbeddedCanvas("glec1", fMain, 700, 500);
  fWF = new TRootEmbeddedCanvas("glec2", fMain, 700, 500);

  fMain->AddFrame(fED, new TGLayoutHints(kLHintsLeft | kLHintsCenterY));
  fMain->AddFrame(fWF, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY));
  TGHorizontalFrame *hfrm = new TGHorizontalFrame(this, 10, 10);

  AddFrame(fMain, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY));

  // Exit
  fButtonExit = new TGTextButton(hfrm, "        &Exit...        ", 3);
  fButtonExit->Connect("Clicked()" , "TApplication", gApplication, "Terminate()");
  hfrm->AddFrame(fButtonExit, new TGLayoutHints(kLHintsCenterX | kLHintsRight,
   10, 10, 10, 10));

  // next event
  fNextEvent = new TGTextButton(hfrm, "        &Next        ", 3);
  fNextEvent->Connect("Clicked()" , "MyMainFrame", this, "NextEvent()");
  hfrm->AddFrame(fNextEvent, new TGLayoutHints(kLHintsCenterX | kLHintsRight,
   10, 10, 10, 10));

  // previous event
  fPrevEvent = new TGTextButton(hfrm, "        &Previous        ", 3);
  fPrevEvent->Connect("Clicked()" , "MyMainFrame", this, "PrevEvent()");
  hfrm->AddFrame(fPrevEvent, new TGLayoutHints(kLHintsCenterX | kLHintsRight,
   10, 10, 10, 10));

  // do draw
  fButtonDraw = new TGTextButton(hfrm, "        &Draw        ", 3);
  fButtonDraw->Connect("Clicked()" , "MyMainFrame", this, "UpdateNumber()");
  hfrm->AddFrame(fButtonDraw, new TGLayoutHints(kLHintsCenterX | kLHintsRight,
   10, 10, 10, 10));

  fNumber = new TGNumberEntry(this, 0, 9,999, TGNumberFormat::kNESInteger,
                                               TGNumberFormat::kNEANonNegative,
                                               TGNumberFormat::kNELLimitMinMax,
                                               0, 99999);
  // fEntry = new TGTextEntry(this);
  // hfrm->AddFrame(fEntry, new TGLayoutHints(kLHintsCenterX | kLHintsRight, 10, 10, 10, 10));
  hfrm->AddFrame(fNumber, new TGLayoutHints(kLHintsBottom | kLHintsRight, 10, 10, 10, 10));

  AddFrame(hfrm, new TGLayoutHints(kLHintsTop | kLHintsExpandX, 5, 5, 5, 5));

  // What to clean up in destructor
  SetCleanup(kDeepCleanup);

   // Set a name to the main frame
  SetWindowName("Event display");
  MapSubwindows();
  Resize(GetDefaultSize());
  MapWindow();
}

void MyMainFrame::NextEvent() {
  ++eventID;
  fNumber->SetIntNumber(eventID);
  DoDraw();
}

void MyMainFrame::PrevEvent() {
  --eventID;
  fNumber->SetIntNumber(eventID);
  DoDraw();
}

void MyMainFrame::UpdateNumber() {
  eventID = fNumber->GetNumberEntry()->GetIntNumber();
  DoDraw();
}

void MyMainFrame::DoDraw() {
  if (eventID > t->GetEntries()) {
    std::cout << "EOF" << std::endl;
    return;
  }
  // get canvas and connect to monitor
  f_ED_canvas = fED->GetCanvas();
  f_ED_canvas->Connect("ProcessedEvent(Int_t,Int_t,Int_t,TObject*)","MyMainFrame",this,
   "EventInfo(Int_t,Int_t,Int_t,TObject*)");

  // read event
  t->GetEntry(eventID);
  std::cout << "Event\t" << eventID << std::endl;
  MM->Reset();
  for (auto x = 0; x < 36; ++x) {
    for (auto y = 0; y < 32; ++y) {
      auto max = 0;
      for (auto t = 0; t < 511; ++t) {
        int Q = 0;
        if (saclay_cosmics)
          Q = padAmpl_saclay[x][y][t] - 250;
        else
          Q = padAmpl[x][y][t] - 250;
        if (Q > max) {
          max = Q;
        }
      } // over t
      if (max)
        MM->Fill(x, y, max);
    }
  }
  f_ED_canvas->cd();
  gStyle->SetOptStat(0);
  MM->Draw("colz");
  // MM->GetXaxis()->SetNdivisions(38);
  // MM->GetXaxis()->SetLabelSize(0.025);
  // MM->GetYaxis()->SetNdivisions(36);
  gPad->SetGrid();
  f_ED_canvas->Update();
}

void MyMainFrame::EventInfo(Int_t event, Int_t px, Int_t py, TObject *selected)
{
  TString s = selected->GetObjectInfo(px,py);
  Ssiz_t first = 0, last;

  first = s.Index("=", 0);
  last  = s.Index(",", first);

  if (first == kNPOS || last == kNPOS)
    return;

  Int_t x = int(TString(s(first+1, last-1)).Atof());

  first = s.Index("=", last);
  last  = s.Index(",", first);

  if (first == kNPOS)
    return;

  Int_t y = int(TString(s(first+1, last-1)).Atof());

  //std::cout << x << "  " << y << std::endl;

  f_WF_canvas = fWF->GetCanvas();
  f_WF_canvas->Clear();
  f_WF_canvas->Divide(3, 3);
  for (auto i = 0; i < 9; ++i) {
    WF[i]->Reset();
    for (auto t_id = 0; t_id < 510; ++t_id) {
      if (x+1 > 35 || x-1 < 0 || y+1 > 31 || y-1 < 0)
        continue;
      int WF_signal = 0;
      if (saclay_cosmics)
        WF_signal = padAmpl_saclay[x-1+i%3][y+1-i/3][t_id] - 250;
      else
        WF_signal = padAmpl[x-1+i%3][y+1-i/3][t_id] - 250;
      if (WF_signal > -250)
        WF[i]->SetBinContent(t_id, WF_signal);
      else
        WF[i]->SetBinContent(t_id, 0);
    }

    f_WF_canvas->cd(i+1);
    WF[i]->GetYaxis()->SetRangeUser(-300, 3000);
    WF[i]->GetXaxis()->SetRangeUser(WFstart, WFend);
    WF[i]->Draw("hist");
    f_WF_canvas->Modified();
    f_WF_canvas->Update();
  }

  f_ED_canvas = fED->GetCanvas();
  f_ED_canvas->cd();

  TBox box(x-1, y-1, x+2, y+2);
  box.SetFillStyle(0);
  box.SetLineColor(kRed);
  box.SetLineWidth(3);
  box.Draw();
  f_ED_canvas->Update();

}

//______________________________________________________________________________
MyMainFrame::~MyMainFrame()
{
   // Clean up all widgets, frames and layouthints that were used
 Cleanup();
}

//______________________________________________________________________________
void MyMainFrame::CloseWindow()
{
   // Called when window is closed via the window manager.

 delete this;
}

//------------------------------------------------------------------------------
void EventDisplay(TString name="")
{
   // Main function (entry point)
  // TString name = "haha";

 new MyMainFrame(gClient->GetRoot(), 1000, 1000, name);
}