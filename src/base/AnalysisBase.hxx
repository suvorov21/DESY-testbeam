#ifndef SRC_BASE_ANALYSISBASE_HXX_
#define SRC_BASE_ANALYSISBASE_HXX_

/** @cond */
#include <numeric>
#include <utility>

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

class TCluster;

/// Class that keeps rules for track clusterisation
class Clustering {
 public:
    Clustering(Float_t a,
               int n_pads
    ) : angle(a), n_pads(n_pads) { coeff = (n_pads == 0) ? 0.001 : 1. / n_pads; }
    /// angle of a reference frame rotation
    Double_t angle;
    /// Number of pads in a row
    int n_pads;
    /// Slope coefficient. 0 corresponds to columns/rows. 1 to diagonals and so on
    Double_t coeff;
    /// Function of row and column that is constant for a given cluster
    [[nodiscard]] int GetConstant(int row, int col) const {
        if (n_pads == 0)
            return col;
        else {
            return floor(coeff * col - row);
        }
    }
};

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
    void setInvert(const bool var) { _invert = var; }

    void setClusterisation(std::shared_ptr<Clustering> var) { _clustering = std::move(var); }

    std::shared_ptr<Clustering> CL_col{nullptr};
    std::shared_ptr<Clustering> CL_diag{nullptr};
    std::shared_ptr<Clustering> CL_2by1{nullptr};
    std::shared_ptr<Clustering> CL_3by1{nullptr};
//  Clustering* CL_3by2;

    /// An actual clustering procedure
    std::shared_ptr<Clustering> _clustering;

    /// Split track into clusters
    /** Extract the vector of clusters from the whole track.
    * the logic of clusterisation is given with the function of the Clustering object
    * The function takes (row, column) and return a value
    * that is constant for a given cluster.
    * For example for clustering with columns the rule column == const is constant.
    * For diagonals column - row = const and so on.
    */
    [[nodiscard]] TClusterPtrVec ClusterTrack(const THitPtrVec &tr) const;

    /************************** Utilities functions *****************************/

    /// Dump progress in the command line
    virtual void CL_progress_dump(int eventID, int Nevents);

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
    [[nodiscard]] std::vector<Int_t> GetEventList() const { return _eventList; }

    /// Draw the selected event
    std::unique_ptr<TCanvas> DrawSelection(
        const std::shared_ptr<TRawEvent> &raw_event,
        const std::shared_ptr<TEvent> &reco_event
    );

    /// verbosity levels
    enum class verbosity_base {
        v_progress = 1,
        v_event_number
    };

 protected:
    /// input file name
    TString _file_in_name{""};
    /// output file name
    TString _file_out_name{""};

    std::unique_ptr<Interface> _interface;

    /// name of the parameter file
    TString _param_file_name{""};

    /// CLI parser
    CmdLineParser _clParser;

    /// vector of event IDs that will be analysed
    std::vector<Int_t> _eventList{};
    /// Whether to store particular event
    bool _store_event{false};
    /// Event to start the loop
    int _start_ID{-999};
    /// Last event to analyse
    int _end_ID{-999};
    /// number of selected events
    int _selected{-999};

    /// number of reconstructed events
    int _reconstructed{-999};

    /// output file
    TFile *_file_out{nullptr};

    /// output vector to put in the file
    std::vector<TObject *> _output_vector{};

    /// Reconstruction used in the analysis.
    /**  You can use plenty in the analysis. At least one should be defines */
    std::unique_ptr<ReconstructionBase> _reconstruction{nullptr};

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
    std::vector<std::pair<int, int>> _broken_pads{};

    /// Minimum number of clusters in the track
    Int_t _min_clusters{0};

    /// Maximum angle (abs(tan)) in the MM plane
    Float_t _max_phi{0};

    /// Maximum angle (abs(tan)) w.r.t. MM plane
    Float_t _max_theta{0};

    /// Hot to treat cross-talk
    enum cross_talk {
        def = 0,
        suppress,
        cherry_pick
    };

    cross_talk _cross_talk_treat{def};

    /// T2K plotting style
    TStyle *_t2kstyle{nullptr};

    /// DEBUG vars
    Int_t _verbose{0};
    bool _batch{false};
    bool _test_mode{false};
    bool _overwrite{false};

    /// Whether to invert track analysis logic
    /** E.g. analyse cosmic tracks. Rows and columns will be replaced */
    bool _invert{false};

    /// Whether to use Gaussian lorentzian PRf fit over polynomial
    bool _gaus_lorentz_PRF{false};

    /// Whether to use individual PRF
    bool _individual_column_PRF{false};

    /// Whether to make PRF center position a free parameter
    bool _prf_free_centre{false};

    /// Whether to use arc function for track fitting
    bool _do_linear_fit{false};
    bool _do_para_fit{false};

    /// Whether to store the WFs
    bool _to_store_wf{false};

    /// Time control system
    TApplication *_app{nullptr};
    /// time controller
    long long _loop_start_ts{0};
    long long _read_time{0};
    long long _reco_time{0};
    long long _ana_time{0};
    long long _column_time{0};
    long long _fitters_time{0};
    long long _filling_time{0};
    long long _sel_time{0};
};

#endif  // SRC_BASE_ANALYSISBASE_HXX_
