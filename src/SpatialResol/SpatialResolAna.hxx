#ifndef SRC_SPATIALRESOL_SPATIALRESOLANA_HXX_
#define SRC_SPATIALRESOL_SPATIALRESOLANA_HXX_

#include "AnalysisBase.hxx"
//#include "CrossingReconstruction.hxx"
#include "DBSCANReconstruction.hxx"
#include "TrackSelection.hxx"
#include "PadSelection.hxx"
#include "TrackFitter.hxx"

#include <set>

constexpr int Nclusters = 200;

/// Spatial resolution analysis
class SpatialResolAna : public AnalysisBase {
 public:
    SpatialResolAna();

    /// Initialise histoes, input files, selections
    bool Initialize() override;

    /// Define the output tree
    void InitializeTree();

    /// Define histoes to be filed
    void InitializeHisto();

    /// Read the file from the previous iteration
    void ReadPrevIter();

    /// Read CLI
    bool ReadCLI(int argc, char **argv) override;

    /// Process the selection output called Event
    bool ProcessEvent(const std::shared_ptr<TEvent> &event) override;

    /// Fit and fill cluster information
    void ProcessCluster(const TClusterPtr &cluster, uint id);

    /// Fit the global track
    std::shared_ptr<TF1> ProcessTrack(const TClusterPtrVec &track);

    /// Treat cross-talk. Either suppress or cherry-pick
    void TreatCrossTalk(const THitPtr &pad, THitPtrVec &robust_pads, int &pad_id);

    /// Fill output tree vars for a pad
    void FillPadOutput(const THitPtr &pad, const uint &clusterId, const int &padId);

    /// Fill SR info
    void FillSR(const TClusterPtr &cluster,
                const uint &clusterId,
                const std::shared_ptr<TF1> &fit,
                const std::array<std::shared_ptr<TF1>, Nclusters> &fit1);
    /// Fill PRF Info
    void FillPRF(const THitPtr &pad,
                 int &padId,
                 const uint &clusterId,
                 const std::shared_ptr<TF1> &fit);

    /// Reset stored vars
    void Reset(int id);

    /// Compute dE/dx
    Double_t ComputedEdx();
    Double_t ComputedEdx_sumWF();
    //for length per pad calculations
    std::pair<double, double> XYintersect(double x1,double y1,double x2,double y2,double intr,double slp);
    bool Select_Length(double x1,double y1,double x2,double y2,double intr,double slp);

    static void DeletePadFromCluster(THitPtrVec &robust_pads, int &pad_id);

    /// Write output files (histos, trees)
    /** Specify only for the values that are not included in the vector */
    bool WriteOutput() override;

    /// Profile PRF with peak and RMS
    static bool ProfilePRF(const TH2F *PRF_h, TGraphErrors *gr);
    /// Profile PRF along X axis
    static bool ProfilePRF_X(const TH2F *PRF_h, TGraphErrors *gr, TH1F *errors);

    /**
     * Initialise PRF with expected params
     * @param name name of the histo
     * @param shift Whether the PRF center position is a free parameter
     * @return initialised PRF histo object
     */
    static TF1 *InitializePRF(const TString &name,
                              bool shift = false,
                              bool gaus_lorentz_PRF = false);

    /// Draw the histograms of interest
    bool Draw();

    /// verbosity levels
    enum class verbosity_SR {
        v_analysis_steps = static_cast<int>(verbosity_base::v_event_number) + 1,
        v_fit_details,
        v_residuals,
        v_prf
    };

 private:
    /// Previous iteration output to extract PRF
    TFile *_prev_iter_file{nullptr};

    /// Name of the file from previous iteration
    TString _prev_iter_name{""};
    /// iteration number. Starting from 0
    Int_t _iteration{0};

    /// Whether to apply correction of spatial resolution (take geometrical mean)
    bool _correction{false};

    /// selection that will be used
    std::unique_ptr<TrackSel> _selection;

    /// output tree
    TTree *_tree{nullptr};
    /// Output tree vars
    // Event vars
    /// event number
    Int_t _ev{-999};
    /// track number in the event
    Int_t _track;
    /// angle in MM plane
    Double_t _angle_xy{-999};
    /// angle w.r.t. MM
    Double_t _angle_yz{-999};

    /// Number of robust clusters in event
    Int_t _rob_clusters{0};

    /// track fit quality Chi2/NDF
    Double_t _quality{-999};
    /// momentum
    Double_t _mom{-999};
    /// angle
    Double_t _sin_alpha{-999};
    /// offset
    Double_t _offset{-999};

    /// maximum of the multiplicity
    Int_t _m_max{-999};

    /// mean multiplicity
    Double_t _m_mean{-999};

    /// Cluster vars
    /// Position of the cluster
    Double_t _clust_pos[Nclusters]{0};
    /// Error from the PRF method fit
    Double_t _clust_pos_error[Nclusters]{0};
    /// X position of the cluster
    Double_t _x[Nclusters]{0};
    /// X position of the avareged cluster
    Double_t _x_av[Nclusters]{0};
    /// Position of the track
    Double_t _track_pos[Nclusters]{0};
    /// Residuals (X_track-X_cluster)
    Double_t _residual[Nclusters]{0};
    /// Residuals (X_track-X_cluster) w/o the given cluster in the fit
    Double_t _residual_corr[Nclusters]{0};
    /// charge in the cluster
    Int_t _charge[Nclusters]{0};
    /// multiplicity of the cluster
    Int_t _multiplicity[Nclusters]{0};
    /// MM module number
    Int_t _module[Nclusters]{0};

    /// Track analytical fit function for plotting
    std::shared_ptr<TF1> _track_fit_func{nullptr};

    /// Pad vars
    /// X_track - X_pad --> X axis of the PRF
    Double_t _dx[Nclusters][10]{0};
    /// Fraction of charge Q_pad / Q_cluster --> Y axis of PRF
    Double_t _qfrac[Nclusters][10]{0};
    /// time of the pad
    Int_t _time[Nclusters][10]{0};

    /** variables for the dEdx analysis */
    /// dE/dx
    Double_t _dEdx{0};
    Double_t _dEdx_sumWF{0};
    Int_t _pad_charge[Nclusters][10]{0};
    Int_t _pad_time[Nclusters][10]{0};

    /// pad X position == column
    Int_t _pad_x[Nclusters][10]{0};
    /// pad Y position === row
    Int_t _pad_y[Nclusters][10]{0};

    /// Waveform width
    Int_t _wf_width[Nclusters][10]{0};
    /// Waveform width at half maximum
    Int_t _wf_fwhm[Nclusters][10]{0};

    Int_t _pad_wf_q[Nclusters][10][520]{0};


    //Track length per pad
    Float_t _pad_lenTr[Nclusters][10];
    //Track length per cluster
    Float_t _cluster_lenTr[Nclusters];
    // To store WFs max of the sum per cluster
    Int_t _cluster_WF_q[Nclusters];
    // To store time variable of WFs max of the sum per cluster
    Int_t _cluster_WF_t[Nclusters];

    /** Histograms **/
    /** Pad response function block **/
    /// PRF function from the previous step. Used for Chi2 fit
    TF1 *_prf_function{nullptr};
    TF1 **_prf_function_arr{nullptr};
    /// PRF histoes
    TH2F *_prf_histo{nullptr};
    // PRF profiling graphs
    TGraphErrors *_prf_graph{nullptr};

    /// Pad response function in time
    TH2F *_prf_time{nullptr};

    /// analytical PRF time function
    TF1 *_prf_time_func{nullptr};

    /// uncertainties of the time profile
    TGraphErrors *_prf_time_error{nullptr};

    TH1F *_prf_time_e{nullptr};

    /// Fitter class for the track and cluster fitting
    std::unique_ptr<TrackFitCern> _fitter{nullptr};

    /** Switchers **/
    /// Whether fit all the pads separately
    bool _do_separate_pad_fit{false};

    /// Whether to assign uncertainty to charge
    bool _charge_uncertainty{false};

    /// Whether to follow the selected events from previous iteration
    bool _processAll{false};

    /// errors vs the PRF value
    // TODO consider removing
    static const int prf_error_bins = 10;
    Double_t prf_error_bins_arr[prf_error_bins] = {0., 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 1.01};
    TAxis *_prf_error_axis = new TAxis(prf_error_bins - 1, prf_error_bins_arr);
    TGraphErrors *_uncertainty_vs_prf_gr{nullptr};
    TGraphErrors *_uncertainty_vs_prf_gr_prev{nullptr};
    TH1F *_uncertainty_vs_prf_histo{nullptr};

    /// Average uncertainty from the previous iteration
    Double_t _uncertainty{-999};

    /// vector of events IDs that passed the Reco and selection
    std::set<UInt_t> _passed_events{};

    // [units are meters]
    static constexpr float prf_min = -0.027;
    static constexpr float prf_max = 0.027;
    static constexpr int prf_bin = 180;

    static constexpr float resol_min = -0.004;
    static constexpr float resol_max = 0.004;
    static constexpr int resol_bin = 200.;

    /// space limit for the PRF usage
    /// pads that are far away are supposed to be unreliable
    const float fit_bound_left = -0.025;
    const float fit_bound_right = 0.025;
};

#endif  // SRC_SPATIALRESOL_SPATIALRESOLANA_HXX_
