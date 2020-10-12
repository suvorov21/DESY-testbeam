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

#include <iostream>


class MyMainFrame : public TGMainFrame {

public:
  MyMainFrame(const TGWindow *p, UInt_t w, UInt_t h, TString name);
  virtual ~MyMainFrame();
  void     CloseWindow();

  void DoDraw();
  void NextEvent();
  void PrevEvent();

  void EventInfo(Int_t event, Int_t px, Int_t py, TObject *selected);

  ClassDef(MyMainFrame, 0)

protected:
  // file
  TFile* f;
  TTree* t;
  Int_t padAmpl[36][32][511];
  Int_t eventID = 0;

  TH2F* MM;
  TH1F* WF[9];

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
};

//______________________________________________________________________________
MyMainFrame::MyMainFrame(const TGWindow *p, UInt_t w, UInt_t h, TString name) :
TGMainFrame(p, w, h)
{
  // read file
  f = new TFile(name, "READ");
  t = (TTree*)f->Get("tree");
  t->SetBranchAddress("PadAmpl", padAmpl);
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


  // do draw
  fButtonDraw = new TGTextButton(hfrm, "        &Draw        ", 3);
  fButtonDraw->Connect("Clicked()" , "MyMainFrame", this, "DoDraw()");
  hfrm->AddFrame(fButtonDraw, new TGLayoutHints(kLHintsCenterX | kLHintsRight,
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
  DoDraw();
}

void MyMainFrame::PrevEvent() {
  --eventID;
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
  MM->Reset();
  for (auto x = 0; x < 36; ++x) {
    for (auto y = 0; y < 32; ++y) {
      auto max = 0;
      for (auto t = 0; t < 511; ++t) {
        if (padAmpl[x][y][t] > max) {
          max = padAmpl[x][y][t];
        }
      } // over t
      if (max)
        MM->Fill(x, y, max);
    }
  }
  f_ED_canvas->cd();
  gStyle->SetOptStat(0);
  MM->Draw("colz");
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
    for (auto t_id = 0; t_id < 511; ++t_id) {
      if (x+1 > 35 || x-1 < 0 || y+1 > 31 || y-1 < 0)
        continue;
      WF[i]->SetBinContent(t_id, padAmpl[x-1+i%3][y+1-i/3][t_id]);
    }

    f_WF_canvas->cd(i+1);
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