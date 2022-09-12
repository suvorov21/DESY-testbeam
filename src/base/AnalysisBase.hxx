#ifndef SRC_BASE_ANALYSISBASE_HXX_
#define SRC_BASE_ANALYSISBASE_HXX_

/** @cond */
#include <numeric>
#include <utility>
#include <list>

#include "TString.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TChain.h"
#include "TApplication.h"
#include "TCanvas.h"
#include <TNtuple.h>
#include <TF1.h>
#include <TH2F.h>
#include <TH3F.h>
#include "TStopwatch.h"
#include "TGraphErrors.h"

#include "SetT2KStyle.hxx"
#include "CmdLineParser.h"
/** @endcond */

#include "ReconstructionBase.hxx"
#include "Interface.hxx"
#include "ClDump.hxx"

class TCluster;

/// Main analysis template
/**
* The basic class that is used for the track analysis. It's responsible for
* opening the data file, loop over events, call the reconstruction,
* recording the output.
*/
class AnalysisBase {
 public:
    AnalysisBase();

    /// Read CL arguments
    virtual bool ReadCLI(int argc, char **argv);

    /// Initialise histoes, input files, selections
    virtual bool Initialize();

    /// Loop over TChain entries. Can use pre-defined event list
    /** Normally for the first iteration loop over all events is performed.
    * For the following iterations only events that passed reconstruction
    * and selection are used.
    */
    virtual bool Loop();
    /// Process the reconstruction output
    /**
    * Function should be defined in the derived analysis. It is responsible
    * for applying the selection for the successfully reconstructed tracks,
    * performing the analysis itself, fill Histograms and TTree branches.
    * @param event The event that is returned by pattern recognition
    */
    virtual bool ProcessEvent(const std::shared_ptr<TEvent> &event);

    /// Write the output file (histos, trees)
    virtual bool WriteOutput();

    /// Read parameter file
    bool ReadParamFile();

    // setters
    void setInputFile(const TString &var) { _file_in_name = var; }
    void setOutputFile(const TString &var) { _file_out_name = var; }
    void setParamFile(const TString &var) { _param_file_name = var; }
    void setVerbosity(const int &var) { _verbose = var; }
    void setStartID(const int &var) { _start_ID = var; }
    void setEndID(const int &var) { _end_ID = var; }
    void setBatchMode(const bool var) { _batch = var; }
    void setDebugMode(const bool var) { _test_mode = var; }
    void setOverwrite(const bool var) { _overwrite = var; }

    /************************** Utilities functions *****************************/

    /// Set a vector of events that will be processed.
    /** Used in the analysis with few iterations. After the first iteration
    * the list of files that passed the selection is established and
    * there is no need to go through all events again, but only through
    * those who passed the reconstruction and selection.
    */
    void SetEventList(const std::vector<Int_t> &var) {
        _eventList.clear();
        _eventList = var;
    }

    /// Draw the selected event
    static std::unique_ptr<TCanvas> DrawSelection(
        const std::shared_ptr<TEvent> &reco_event
    );

    /// verbosity levels
    enum class verbosity_base {
        v_progress = 1,
        v_event_number
    };
 protected:
    /// output file name
    TString _file_out_name{""};

    /// output file
    std::unique_ptr<TFile> _file_out{nullptr};

    /// output vector to put in the file
    std::vector<TObject *> _output_vector{};

    /// Reconstruction used in the analysis.
    /**  You can use plenty in the analysis. At least one should be defines */
    std::unique_ptr<ReconstructionBase> _reconstruction{nullptr};

    /// CLI parser
    CmdLineParser _clParser;

    /// DEBUG vars
    Int_t _verbose{0};
    bool _batch{false};
    bool _test_mode{false};

    /// Whether to make PRF center position a free parameter
    bool _prf_free_centre{false};
    /// Whether to use Gaussian lorentzian PRf fit over polynomial
    bool _gaus_lorentz_PRF{false};

    /// Selection parameters
    /// The maximum multiplicity of the track
    Int_t _max_mult{0};

    /// The maximum mean multiplicity for track
    Float_t _max_mean_mult{0};

    /// Whether to cut tracks with gap in the cluster
    /** E.g. one missed pad in the middle of the column or in the
    * middle of the diagonal
    */
    bool _cut_gap{false};

    /// Vector of broken pads to be excluded from the analysis
    broken_pads_t _broken_pads{};

    /// Minimum number of clusters in the track
    Int_t _min_clusters{0};

    /// Maximum angle (abs(tan)) in the MM plane
    Float_t _max_phi{0};

    /// Maximum angle (abs(tan)) w.r.t. MM plane
    Float_t _max_theta{0};

    /// Time boundaries of the selection
    Int_t _time_min{-1};
    Int_t _time_max{-1};

    /// Whether to invert track analysis logic
    /** E.g. analyse cosmic tracks. Rows and columns will be replaced */
    bool _invert{false};

    /// Whether to use individual PRF
    bool _individual_column_PRF{false};

    /// Whether to use arc function for track fitting
    bool _do_linear_fit{false};
    bool _do_para_fit{false};

    /// Whether to store the WFs
    bool _to_store_wf{false};

    /// An actual clustering procedure
    std::unique_ptr<Clustering> _clustering;

    /// Whether to store particular event
    bool _store_event{false};

    /// vector of event IDs that will be analysed
    std::vector<Int_t> _eventList{};

    /// Hot to treat cross-talk
    enum class cross_talk {
        defaultCt = 0,
        suppress,
        cherry_pick
    };

    cross_talk _cross_talk_treat{cross_talk::defaultCt};

    /// number of selected events
    int _selected{-999};

    /// number of reconstructed events
    int _reconstructed{-999};

    /// time controller
    long long _read_time{0};
    long long _reco_time{0};
    long long _ana_time{0};
    long long _column_time{0};
    long long _fitters_time{0};
    long long _filling_time{0};
    long long _sel_time{0};

 private:
    /// input file name
    TString _file_in_name{""};

    /// Number of data readers thread
    uint _readerThreads{1};
    mutable std::mutex _mu;
    /// vector of interfaces for the data reading
    std::vector<std::unique_ptr<Interface>> _interface;

    /// list of raw events to be filled in parallel
    std::list<std::shared_ptr<TEvent>> _TEventList;

    /// name of the parameter file
    TString _param_file_name{""};

    /// Event to start the loop
    int _start_ID{-999};
    /// Last event to analyse
    int _end_ID{-999};

    /// T2K plotting style
    std::unique_ptr<TStyle> _t2kstyle{nullptr};

    bool _overwrite{false};

    /// Dump the progress into CL
    ClDump _clDump{};

    /// Time control system
    TApplication* _app{nullptr};
};

#endif  // SRC_BASE_ANALYSISBASE_HXX_
