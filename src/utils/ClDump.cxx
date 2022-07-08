//
// Created by SERGEY SUVOROV on 11/04/2022.
//

#include <GenericToolbox.h>

#include "ClDump.hxx"
#include <Geom.hxx>

//******************************************************************************
void ClDump::CL_progress_dump(int eventID,
                              int N_events,
                              int selected) {
//******************************************************************************
    auto mem = GenericToolbox::getProcessMemoryUsage();
    if (eventID)
        _loop_start_ts += GenericToolbox::getElapsedTimeSinceLastCallInMicroSeconds("loop");
    else
        GenericToolbox::getElapsedTimeSinceLastCallInMicroSeconds("loop");

    int m, s;
    if (eventID) {
        long long EET = num::cast<long long>(((N_events - eventID) * _loop_start_ts / eventID));
        EET /= 1e6;
        m = num::cast<int>(EET / 60);
        s = num::cast<int>(EET % 60);
    }

    for (auto i = 0; i < 30; ++i)
        if (i < 30. * eventID / N_events) std::cout << "#";
        else std::cout << " ";
    std::cout << "]   Nevents = " << N_events << "\t" << round(1. * eventID / N_events * 100) << "%";
    std::cout << "\t Selected  " << selected << " (" << round(1. * selected / eventID * 100) << "%)";
    std::cout << "\t Memory  " << std::setw(4) << mem / 1048576 << " " << "MB";
    if (eventID) {
        std::cout.precision(4);
        std::cout << "\t Av speed " << std::setw(5) << num::cast<double>(_loop_start_ts) / 1e3 / eventID << " ms/event";
        std::cout << "\t EET " << std::setw(2) << m << ":";
        std::cout << std::setw(2) << std::setfill('0');
        std::cout << s;
        std::cout << std::setfill(' ');
    }
    std::cout << "      \r[" << std::flush;
}
