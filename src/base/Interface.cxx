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

    if (file.Get<TTree>("event_tree")) {
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

template<int timeSize>
void interfaceRoot<timeSize>::Initialize() {
    Interface::chainInputFiles("tree");
    _chain->SetBranchAddress("PadAmpl", _padAmpl);
}

template<int timeSize>
std::shared_ptr<TRawEvent> interfaceRoot<timeSize>::getEvent(Int_t i) {
    // create TRawEvent from 3D array
    std::shared_ptr<TRawEvent> event = std::make_shared<TRawEvent>(i);
    _chain->GetEntry(i);
    // Subtract the pedestal
    for (auto x = 0; x < Geom::nPadx; ++x) {
        for (auto y = 0; y < Geom::nPady; ++y) {
            auto elec = Geom::get().GetPadToEle(x, y);
            auto hit = new TRawHit(0, elec.first, elec.second);
            hit->ResetWF();
            int max = -260;
            for (auto t = 0; t < timeSize; ++t) {
                hit->SetADCunit(t, _padAmpl[x][y][t]);
                if (_padAmpl[x][y][t] > max)
                    max = _padAmpl[x][y][t];
            } // time
            hit->ShrinkWF();
            if (max > 0)
                event->AddHit(hit);
            else
                delete hit;
        } // over Y
    } // over X
    return event;
}

template class interfaceRoot<510>;
template class interfaceRoot<511>;

void interfaceRawEvent::Initialize() {
    Interface::chainInputFiles("event");
    _event = new TRawEvent();
    _chain->SetBranchAddress("Event", &_event);
}

std::shared_ptr<TRawEvent> interfaceRawEvent::getEvent(Int_t i) {
    _chain->GetEntry(i);
    return std::make_shared<TRawEvent>(_event);
}

