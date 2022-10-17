#include <algorithm>
#include <memory>
#include <cstdlib>
#include <iomanip>

#include "TROOT.h"

#include <GenericToolbox.h>

#include "AnalysisBase.hxx"

//******************************************************************************
AnalysisBase::AnalysisBase() :
    _clustering(std::make_unique<Clustering>(clusterType::kRowColumn, false)),
    _start_ID(-1),
    _end_ID(-1),
    _selected(0),
    _reconstructed{0},
    _reconstruction(nullptr),
    _max_mult(6),
    _max_mean_mult(5),
    _cut_gap(true),
    _min_clusters(30),
    _verbose(1),
    _to_store_wf(true) {
//******************************************************************************

    // CLI reader
    _clParser.setIsUnixGnuMode(true);
    _clParser.setIsFascist((true));

    _clParser.addOption("input_file", {"-i", "--input"}, "Input file name", 1);
    _clParser.addOption("output_file", {"-o", "--output"}, "Output file name", 1);
    _clParser.addOption("nthreads", {"--nthread"}, "Number of reader threads [1-4]");

    _clParser.addOption("param_file", {"-p", "--param"}, "Parameter file name", 1);
    _clParser.addOption("verbosity", {"-v", "--verbose"}, "Verbosity level", 1);
    _clParser.addOption("start_id", {"--start"}, "Start event ID", 1);
    _clParser.addOption("end_id", {"--end"}, "End event ID", 1);

    _clParser.addTriggerOption("batch", {"-b"}, "Batch mode");
    _clParser.addTriggerOption("debug", {"-d", "--debug"}, "Debug mode");
    _clParser.addTriggerOption("overwrite", {"-r", "--overwrite"}, "Overwrite output file");

    _clParser.addTriggerOption("help", {"-h", "--help"}, "Print usage");
}

bool AnalysisBase::ReadCLI(int argc, char **argv) {
    _clParser.parseCmdLine(argc, argv);

    if (_clParser.isOptionTriggered("help")) {
        std::cout << _clParser.getConfigSummary();
        exit(1);
    }

    setInputFile(_clParser.getOptionVal<TString>("input_file", "", 0));
    setOutputFile(_clParser.getOptionVal<TString>("output_file", "", 0));
    _readerThreads = _clParser.getOptionVal<uint>("nthreads", 1, 0);

    setParamFile(_clParser.getOptionVal<TString>("param_file", "", 0));
    setVerbosity(_clParser.getOptionVal<int>("verbosity", _verbose, 0));

    setStartID(_clParser.getOptionVal<int>("start_id", _start_ID, 0));
    setEndID(_clParser.getOptionVal<int>("end_id", _end_ID, 0));

    setBatchMode(_clParser.isOptionTriggered("batch"));
    setDebugMode(_clParser.isOptionTriggered("debug"));
    setOverwrite(_clParser.isOptionTriggered("overwrite"));

    if (!_batch)
        _app = new TApplication("app", &argc, argv);

    return true;
}

//******************************************************************************
bool AnalysisBase::Initialize() {
//******************************************************************************
    // Read parameter file
    if (!ReadParamFile()) {
        std::cerr << "ERROR! " << __func__ << "(). Parameter file is not read" << std::endl;
        exit(1);
    }

    if (_file_out_name == "") {
        std::cerr << "ERROR. " << __func__ << "() No output file specified" << std::endl;
        exit(1);
    }

    TString filename = Interface::getRootFile(_file_in_name);
    _interface.reserve(_readerThreads - 1);
    for (auto i = 0; i < std::max(1U, _readerThreads - 1); ++i) {
        auto interfaceType = Interface::getFileType(filename);
        _interface.emplace_back(interfaceFactory::get(filename, interfaceType));
        _interface.back()->Initialize();
    }

    std::cout << "Initializing analysis base...............";

    // setup the T2K style
    Int_t T2KstyleIndex = 2;
    // Official T2K style as described in http://www.t2k.org/comm/pubboard/style/index_html
    TString localStyleName = "T2K";
    // -- WhichStyle --
    // 1 = presentation large fonts
    // 2 = presentation small fonts
    // 3 = publication/paper
    _t2kstyle = std::unique_ptr<TStyle>(T2K().SetT2KStyle(T2KstyleIndex, localStyleName));

    gROOT->SetStyle(_t2kstyle->GetName());
    gROOT->ForceStyle();

    Long64_t N_events = _interface[0]->getEntries();
    _eventList.reserve(N_events);
    for (auto i = 0; i < N_events; ++i)
        _eventList.push_back(i);

    // Open the output file
    if (_overwrite)
        _file_out = std::unique_ptr<TFile>(TFile::Open(_file_out_name.Data(), "RECREATE"));
    else
        _file_out = std::unique_ptr<TFile>(TFile::Open(_file_out_name.Data(), "NEW"));

    if (!_file_out->IsOpen()) {
        std::cerr << "ERROR. AnalysisBase::Initialize()" << std::endl;
        std::cerr << "File already exists or directory is not writable" << std::endl;
        std::cerr << "To prevent overwriting of the previous result the program will exit" << std::endl;
        exit(1);
    }

    if (_file_out)
        _file_out->cd();

    // Initialize histoes
    // * do it in your analysis *

    // Initialize selection
    // * do it in your analysis *

    std::cout << "done" << std::endl;

    return true;
}

//******************************************************************************
bool AnalysisBase::Loop() {
//******************************************************************************
    auto N_events = num::cast<int>(_eventList.size());
    if (_test_mode)
        N_events = std::min(num::cast<int>(_eventList.size()), 100);

    if (_start_ID < 0)
        _start_ID = 0;
    if (_end_ID > 0) {
        _end_ID = std::min(_end_ID, N_events);
    } else
        _end_ID = N_events;

    if (_verbose >= static_cast<int>(verbosity_base::v_progress)) {
        std::cout << "Input file............................... " << _file_in_name << std::endl;
        std::cout << "Output file.............................. " << _file_out_name << std::endl;
        std::cout << "Processing" << std::endl;
        std::cout << "[                              ]   Nevents = " << _end_ID - _start_ID << "\r[";
    }

    int denominator = 100;
    if (N_events < 100)
        denominator = N_events;

    std::atomic<int> read{_start_ID};
    std::atomic<int> readAhead{0};
    std::vector<std::thread> inputReader;
    inputReader.reserve(_readerThreads - 1);
    for (auto i = 0; i < _readerThreads - 1; ++i) {
        inputReader.emplace_back([interface = std::move(_interface[i]), &read, &readAhead, this]() {
            while (read < _end_ID) {
                // do not read all events in a file to save RAM
                if (readAhead > 50)
                    continue;
                int toRead = read++;
                auto event = interface->getEvent(_eventList[toRead]);
                std::lock_guard<std::mutex> lock(_mu);
                _TEventList.emplace_back(std::move(event));
                ++readAhead;
            }
        });
    }

    // Event loop
    int eventID = _start_ID;
    int prevDump = eventID - 1;
    while (eventID < _end_ID) {
        if (_verbose >= static_cast<int>(verbosity_base::v_event_number)) {
            std::cout << "*************************************" << std::endl;
            std::cout << "Event " << eventID << std::endl;
            std::cout << "*************************************" << std::endl;
        }

        // Dump progress in command line
        if ((eventID % (N_events / denominator)) == 0 && _verbose == static_cast<int>(verbosity_base::v_progress)) {
            if (prevDump != eventID) {
                _clDump.CL_progress_dump(eventID - _start_ID,
                                         _end_ID - _start_ID,
                                         _selected);
                prevDump = eventID;
            }
        }

        // start the timer
        GenericToolbox::getElapsedTimeSinceLastCallInMicroSeconds("read");

        // in case of one thread, read the interface
        if (_readerThreads == 1) {
            _TEventList.emplace_back(_interface[0]->getEvent(_eventList[eventID]));
        }

        // wait for interface to read an event
        if (_TEventList.empty())
            continue;

        auto tEvent = _TEventList.front();
        ++eventID;
        --readAhead;
        _TEventList.erase(_TEventList.begin());

        _read_time += GenericToolbox::getElapsedTimeSinceLastCallInMicroSeconds("read");

        _store_event = false;

        GenericToolbox::getElapsedTimeSinceLastCallInMicroSeconds("reco");

        // copy event to a child class to be filled with reconstruction
        bool sel = _reconstruction->ReconstructEvent(tEvent);
        _reco_time += GenericToolbox::getElapsedTimeSinceLastCallInMicroSeconds("reco");
        if (!sel) {
            continue;
        }

        ++_reconstructed;

        // do basic plotting
        auto c = std::make_unique<TCanvas>();
        if (!_batch) {
            c = DrawSelection(tEvent);
            c->SetTitle(Form("Event %i", tEvent->GetID()));
            c->Draw();
            c->WaitPrimitive();
        }

        GenericToolbox::getElapsedTimeSinceLastCallInMicroSeconds("ana");
        ProcessEvent(tEvent);
        _ana_time += GenericToolbox::getElapsedTimeSinceLastCallInMicroSeconds("ana");

        if (_store_event)
            ++_selected;
    } // end of event loop
    for (auto i = 0; i < _readerThreads - 1; ++i) {
        inputReader[i].join();
    }
    // if progress bar is active --> go to the next line
    if (_verbose == static_cast<int>(verbosity_base::v_progress))
        std::cout << std::endl;

    return true;
}

//******************************************************************************
bool AnalysisBase::ProcessEvent(const std::shared_ptr<TEvent> &event) {
//******************************************************************************
    (void) event;
    throw std::logic_error("Event processing should be defined in your analysis");
}

//******************************************************************************
bool AnalysisBase::WriteOutput() {
//******************************************************************************
    //if(_test_mode) return true;
    if (!_file_out || !_file_out->IsOpen()) {
        std::cout << "AnalysisBase::WriteOutput   _file_out is not Open!" << std::endl;
        return false;
    }

    std::cout << "Writing standard output..................";

    _file_out->cd();

    auto size = _output_vector.size();
    for (auto i = 0; i < size; ++i) {
        if (!_output_vector[i])
            std::cerr << "ERROR! " << __func__ << "()  output object pointer is nullptr" << std::endl;
        _output_vector[i]->Write();
    }

    _file_out->Close();

    std::cout << "done     " << "Write  " << size << " objects" << std::endl;

    return true;
}

//******************************************************************************
std::unique_ptr<TCanvas> AnalysisBase::DrawSelection(
    const std::shared_ptr<TEvent> &reco_event
) {
//******************************************************************************
    gStyle->SetCanvasColor(0);
    gStyle->SetMarkerStyle(21);
    gStyle->SetMarkerSize(1.05);
    TH2F MM("MM", "", Geom::nPadx, 0, Geom::nPadx, Geom::nPady, 0, Geom::nPady);
    TH2F MMsel("MMsel", "", Geom::nPadx, 0, Geom::nPadx, Geom::nPady, 0, Geom::nPady);
    TNtuple event3D("event3D", "event3D", "x:y:z:c");

    // all hits
    for (const auto &module : reco_event->GetAllHits()) {
        for (const auto& h : module.second) {
            MM.Fill(h->GetCol(), h->GetRow(), h->GetQMax());
        }
    }

    // TrackSel hits
    // FIXME loop over patterns
    for (const auto &pattern : reco_event->GetAllHits()) {
        for (const auto &h : pattern.second) {
            if (!h->GetQMax()) continue;
            event3D.Fill((Float_t) h->GetTime(), (Float_t) h->GetRow(), (Float_t) h->GetCol(), (Float_t) h->GetQMax());
            MMsel.Fill(h->GetCol(), h->GetRow(), h->GetQMax());
        }
    }

    auto canv = std::make_unique<TCanvas>("canv", "canv", 0., 0., 1400., 600.);
    canv->Divide(3, 1);
    canv->cd(1);
    MM.Draw("COLZ");
    canv->cd(2);
    MMsel.Draw("COLZ");

    canv->cd(3);
    event3D.Draw("x:y:z:c", "", "box2");
    auto htemp = (TH3F *) gPad->GetPrimitive("htemp");
    htemp->GetXaxis()->SetLimits(0, Geom::nPadx);
    htemp->GetYaxis()->SetLimits(0, Geom::nPady);
    htemp->GetZaxis()->SetLimits(0, 500);
    htemp->SetTitle("");
    canv->Update();
    return canv;
}

//******************************************************************************
//******************************************************************************
// ***** Functions below are utils: params reading, progress dump **************
//******************************************************************************
//******************************************************************************

//******************************************************************************
bool AnalysisBase::ReadParamFile() {
//******************************************************************************
    if (_param_file_name == "") {
        auto source = std::string(__FILE__);
        auto found = source.find_last_of('/');
        _param_file_name = source.substr(0, found) + "/../../params/default.ini";
    }
    std::cout << "*****************************************" << std::endl;
    std::cout << "Read parameters from " << _param_file_name << std::endl;
    std::ifstream cFile(_param_file_name);
    if (cFile.is_open()) {
        std::string line;
        while (getline(cFile, line)) {
            line.erase(std::remove_if(line.begin(),
                                      line.end(),
                                      isspace
                       ),
                       line.end()
            );

            if (line[0] == '#' || line.empty())
                continue;
            auto delimiterPos = line.find('=');
            auto name = line.substr(0, delimiterPos);
            auto value = line.substr(delimiterPos + 1);

            if (name == "invert") {
                if (value == "1") {
                    _invert = true;
                    _clustering->setInvert(_invert);
                    std::cout << "Inverted Geometry used" << std::endl;
                }
            } else if (name == "cluster") {
                if (value == "column") {
                    std::cout << "Column cluster is used" << std::endl;
                } else if (value == "diag") {
                    _clustering = std::make_unique<Clustering>(clusterType::kDiagonal, _invert);
                    std::cout << "Diagonal cluster is used" << std::endl;
                } else if (value == "2by1") {
                    _clustering = std::make_unique<Clustering>(clusterType::k2by1, _invert);
                    std::cout << "2by1 cluster is used" << std::endl;
                } else if (value == "3by1") {
                    _clustering = std::make_unique<Clustering>(clusterType::k3by1, _invert);
                    std::cout << "3by1 cluster is used" << std::endl;
                } else {
                    std::cerr << "ERROR. Unknown clustering " << value << std::endl;
                    return false;
                }
            } else  if (name == "prf_shape") {
                if (value == "gaus_lorentz") {
                    _gaus_lorentz_PRF = true;
                    std::cout << "PRF is fit with Gaussian-Lorentzian" << std::endl;
                } else if (value == "pol4") {
                    std::cout << "PRF is fit with 4th degree polinom" << std::endl;
                } else {
                    std::cerr << "ERROR. Unknown PRF function " << value << std::endl;
                    return false;
                }
            } else if (name == "individual_prf") {
                if (value == "1") {
                    _individual_column_PRF = true;
                    std::cout << "Individual PRF for each column is used" << std::endl;
                }
            } else if (name == "prf_centre_freedom") {
                if (value == "1") {
                    _prf_free_centre = true;
                    std::cout << "PRF centre position is a free parameter of the fit" << std::endl;
                }
            } else if (name == "track_shape") {
                if (value == "parabola") {
                    _do_para_fit = true;
                    std::cout << "Parabola track fit is used" << std::endl;
                } else if (value == "linear") {
                    _do_linear_fit = true;
                    std::cout << "Linear track fit is used" << std::endl;
                } else if (value == "arc") {
                    std::cout << "Arc track fit is used" << std::endl;
                } else {
                    std::cerr << "ERROR. Unknown track shape " << value << std::endl;
                    return false;
                }
            } else if (name == "max_mult") {
                _max_mult = TString(value).Atoi();
            } else if (name == "max_mean_mult") {
                _max_mean_mult = (Float_t) TString(value).Atof();
            } else if (name == "cut_gap") {
                if (value == "0")
                    _cut_gap = false;
            } else if (name == "cluster_min") {
                _min_clusters = TString(value).Atoi();
            } else if (name == "max_phi") {
                _max_phi = (Float_t) TString(value).Atof();
            } else if (name == "max_theta") {
                _max_theta = (Float_t) TString(value).Atof();
                //switch to WF storage
            } else if (name == "to_store_wf") {
                if (value == "0") {
                    _to_store_wf = false;
                    std::cout << "WFs will NOT be stored" << std::endl;
                } else {
                    std::cout << "WFs will be stored. Analysis will be slowed down" << std::endl;
                }
            } else if (name == "cross_talk") {
                if (value == "suppress") {
                    _cross_talk_treat = cross_talk::suppress;
                    std::cout << "Cross-talk will be suppressed" << std::endl;
                } else if (value == "cherry_pick") {
                    _cross_talk_treat = cross_talk::cherry_pick;
                    std::cout << "Cross-talk will be cherry-picked" << std::endl;
                } else {
                    _cross_talk_treat = cross_talk::defaultCt;
                    std::cout << "Cross-talk will not be treated" << std::endl;
                }
            } else if (name == "dead") {
                auto dead_pads = GenericToolbox::splitString(value, ";");
                if (!dead_pads.empty())
                    std::cout << "Dead pads: ";
                for (const auto &pad : dead_pads) {
                    auto coordinates = GenericToolbox::splitString(pad, ",");
                    if (coordinates.size() != 2) {
                        continue;
                    }
                    _broken_pads.push_back({TString(coordinates[0]).Atoi(),
                                            TString(coordinates[1]).Atoi(),
                                            TString(coordinates[2]).Atoi()}
                    );
                    std::cout << "module " << _broken_pads.back()[0]
                              << "row: " << _broken_pads.back()[1]
                              << "column " << _broken_pads.back()[2] << "; ";
                }
                std::cout << std::endl;
            } else if (name == "time_min") {
                _time_min = (Int_t) TString(value).Atoi();
                std::cout << "Minimum time " << _time_min << std::endl;
            } else if (name == "time_max") {
                _time_max = (Int_t) TString(value).Atoi();
                std::cout << "Maximum time " << _time_max << std::endl;
            }
        }
    } else {
        return false;
    }
    std::cout << "*****************************************" << std::endl;
    return true;
}
