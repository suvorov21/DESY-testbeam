//
// Created by SERGEY SUVOROV on 11/04/2022.
//

#ifndef DESY_TESTBEAM_SRC_UTILS_CL_DUMP_HXX_
#define DESY_TESTBEAM_SRC_UTILS_CL_DUMP_HXX_

class ClDump {
    long long _loop_start_ts;
 public:
    void CL_progress_dump(int eventID,
                          int N_events,
                          int selected);
};

#endif //DESY_TESTBEAM_SRC_UTILS_CL_DUMP_HXX_
