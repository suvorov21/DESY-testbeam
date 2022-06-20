//
// Created by SERGEY SUVOROV on 06/04/2022.
//

#include <iostream>
#include <fstream>

#include "TFile.h"
#include "TTree.h"

#include "Interface.hxx"

TString Interface::getRootFile(const TString& filename) {
    if (filename == "") {
        std::cerr << "ERROR. " << __func__ << "() No input file specified" << std::endl;
        exit(1);
    }

    if (!filename.Contains(".root")) {
        std::ifstream fList(filename.Data());
        if (!fList.good()) {
            std::cerr << "Can not read input " << filename << std::endl;
            exit(1);
        }
        std::string temp;
        getline(fList, temp);
        return static_cast<std::string>(temp);
    }

    return filename;
}

interfaceType Interface::getFileType(const TString & filename) {
    TFile file(filename.Data(), "READ");

    if (file.Get<TTree>("EventTree")) {
        return interfaceType::kRawEvent;
    }

    if (auto tree = file.Get<TTree>("tree")) {
        TString branch_name = tree->GetBranch("PadAmpl")->GetTitle();
        if (branch_name.Contains("[510]")) {
            return interfaceType::kArray510;
        } else if (branch_name.Contains("[511]")) {
            return interfaceType::kArray511;
        }
    }

    throw std::logic_error("Unknown file type");
}

bool Interface::chainInputFiles(const TString& treeName) {
    _chain = std::make_unique<TChain>(treeName);
    if (_file_in_name.Contains(".root")) {
        _chain->AddFile(_file_in_name);
    } else {
        std::ifstream fList(_file_in_name.Data());
        if (!fList.good()) {
            std::cerr << "Can not read input " << _file_in_name << std::endl;
            exit(1);
        }
        while (fList.good()) {
            std::string temp_filename;
            getline(fList, temp_filename);
            if (fList.eof()) break;
            _chain->AddFile(temp_filename.c_str());
        }
    }

    return true;
}

Long64_t Interface::getEntries() {
    if (_chain) {
        return _chain->GetEntries();
    }
    return 0;
}

template<short timeSize>
void interfaceRoot<timeSize>::Initialize() {
    Interface::chainInputFiles("tree");
    _chain->SetBranchAddress("PadAmpl", _padAmpl);
}

template<short timeSize>
std::shared_ptr<TEvent> interfaceRoot<timeSize>::getEvent(Int_t i) {
    // create TRawEvent from 3D array
    std::shared_ptr<TEvent> event = std::make_shared<TEvent>(i);
    _chain->GetEntry(i);
    // Subtract the pedestal
    for (short x = 0; x < Geom::nPadx; ++x) {
        for (short y = 0; y < Geom::nPady; ++y) {
            auto elec = Geom::get().GetPadToEle(x, y);
            short Qmax = -1;
            short Tmax = -1;
            // array is created (and destructed) per WF
            // if the pad has a meaningful signal the array will be cast to vector later
            // creating/destructing a vector at this step leads to significant
            // performance loss
            std::array<short, timeSize> wf{};
            for (short t = 0; t < timeSize; ++t) {
                short q;
                try {
                    q = num::cast<short>(_padAmpl[x][y][t] - 250);
                } catch (const std::bad_cast& e) {}

                if (q < -249)
                    continue;

                wf[t] = q;
                if (q > Qmax) {
                    Qmax = q;
                    Tmax = t;
                }
            } // time

            if (Qmax < 0)
                continue;

            // compute FWHM
            short fwhm = 0;
            short width = 0;
            for (auto t = 0; t < timeSize; ++t) {
                if (wf[t] > Qmax / 2)
                    fwhm += 1;
                if (wf[t] > 0)
                    width += 1;
            }
            auto hit = std::make_shared<THit>(x, y,
                                              0, elec.first, elec.second);
            hit->SetFWHM(fwhm);
            hit->SetADCvector(std::vector<short>(wf.begin(), wf.end()));
            hit->SetWidth(width);
            hit->SetQMax(Qmax);
            hit->SetTimeMax(Tmax);
            hit->ShrinkWF();
            event->AddHitPtr(hit);
        } // over Y
    } // over X
    return event;
}

template class interfaceRoot<510>;
template class interfaceRoot<511>;

void interfaceRawEvent::Initialize() {
    Interface::chainInputFiles("EventTree");
    _event = new TRawEvent();
    _chain->SetBranchAddress("TRawEvent", &_event);
}

std::shared_ptr<TEvent> interfaceRawEvent::getEvent(Int_t i) {
    _chain->GetEntry(i);
    return std::make_shared<TEvent>(*_event);
}

