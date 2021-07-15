#ifndef SRC_BASE_ANALYSISBASE_HXX_
#define SRC_BASE_ANALYSISBASE_HXX_

/** @cond */
#include <numeric>

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
/** @endcond */

#include "ReconstructionBase.hxx"


class TCluster;

/// Class that keeps rules for track clusterisation
class Clustering {
public:
  Clustering(Float_t a,
             int n_pads
             ):angle(a), n_pads(n_pads) {coeff = (n_pads == 0)?0.001:1./n_pads;}
  virtual ~Clustering(){;};
  /// angle of a reference frame rotation
  Float_t angle;
  /// Number of pads in a row
  int n_pads;
  /// Slope coefficient. 0 corresponds to columns/rows. 1 to diagonals and so on
  Float_t coeff;
  /// Function of row and column that is constant for a given cluster
  int GetConstant(int row, int col) {
    if (n_pads == 0)
      return col;
    else {
      return (floor(coeff * col - row));
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
  AnalysisBase(int argc, char** argv);
  virtual ~AnalysisBase() {;}

  /// Initialise histoes, input files, selections
  virtual bool Initialize();
  /// Loop over TChain entries. Can use pre-defined event list
  /** Normally for the first iteration loop over all events is performed.
  * For the following iterations only events that passed reconstruction
  * and selection are used.
  */
  virtual bool Loop(std::vector<Int_t> EventList);
  /// Process the reconstruction output
  /**
  * Function should be defined in the derived analysis. It is responsible
  * for applying the selection for the successfully reconstructed tracks,
  * performing the analysis itself, fill jistograms and TTree branches.
  */
  virtual bool ProcessEvent(const TEvent* event);
  /// Write the output file (histos, trees)
  virtual bool WriteOutput();

  /// Read parameter file
  bool ReadParamFile();

  /// Process a cluster and return only pads that are suggested to be robust
  /** E.g. function can return only 2 pads in a column.
   * Another use case is to ommit pads with wrong timestamps.
   * Any user defined selection may be applied.
   */
  std::vector<THit*> GetRobustPadsInCluster(std::vector<THit*> col);
  /// Return only robust clusters
  /** E.g. apply a trunccation - ommit clusters with relatively large charge
   * Or put a strong upper limit on cluster charge.
   * Any condition can be specified.
   */
  std::vector<TCluster*> GetRobustClusters(std::vector<TCluster*> tr);

  // TODO consider a better implementation. No need to create all instances
  Clustering* CL_col;
  Clustering* CL_diag;
  Clustering* CL_2by1;
  Clustering* CL_3by1;
  Clustering* CL_3by2;

  /// An actual clustering precedure
  Clustering* _clustering;

  /// Split track into clusters
  /** Extract the vector of clusters from the whole track.
  * the logic of clusterisation is given with the function of the Clustering object
  * The function takes (row, column) and return a value
  * that is constant for a given cluster.
  * For example for clustering with columns the rule column == const is constant.
  * For diagonals column - row = connst and so on.
  */
  std::vector<TCluster*> ClusterTrack(const std::vector<THit*> &tr);


  /************************** Utilities functions *****************************/

  /// Print usage
  void help(const std::string& name);

  /// Dump progress in the command line
  virtual void CL_progress_dump(int eventID, int Nevents);

  /// Dump RAM usage in CL
  void process_mem_usage(double& vm_usage, double& resident_set);

  /// Set a vector of events that will be processed.
  /** Used in the analysis with few iterations. After the first iteration
  * the list of files that passde the selection is established and
  * there is no need to go through all events again, but only through
  * those who passed the resonstruction and selection.
  */
  void SetEventList(const std::vector<Int_t>& var) {_EventList.clear(); _EventList = var;}
  std::vector<Int_t> GetEventList() const {return _EventList;}

  /// Draw the selected event
  virtual void DrawSelection(const TEvent *event);

  AnalysisBase(const AnalysisBase& ana){(void)ana;
    std::cerr << "Copy constructor is depricated" << std::endl; exit(1);}
  bool operator==(const AnalysisBase* ana){(void)ana;
    std::cerr << "Comparison is depricated" << std::endl; exit(1);}

  bool ChainInputFiles(TString tree_name);

  /// verbosity levels
  enum verbosity_base {
    v_progress = 1,
    v_event_number,
    v_base_last
  };

 protected:
  /// input file name
  TString _file_in_name;
  /// output file name
  TString _file_out_name;

  /// name of the parameter file
  TString _param_file_name;

  /// vector of event IDs that will be analysed
  std::vector<Int_t> _EventList;
  /// Whether to store particular event
  bool    _store_event;
  /// Event to start the loop
  int     _start_ID;
  /// Last event to analyse
  int     _end_ID;
  /// number of selected events
  int     _selected;

  /// The current processing event
  TRawEvent* _event;
  bool    _work_with_event_file;

  /// input file
  /* May be a single ROOT file or a list of files*/
  TFile* _file_in;
  /// output file
  TFile* _file_out;

  /// chain with inpout files
  TChain* _chain;

  /// Name of the file from previous iteration
  TString _Prev_iter_name;
  /// iteration number. Starting from 0
  Int_t   _iteration;

  /// Whether to apply correction of spatial resolution (take geometrical mean)
  bool _correction;

  /// choose the input array size, whether to use 510 or 511 time bins
  bool _saclay_cosmics;
  Int_t _padAmpl[geom::nPadx][geom::nPady][geom::Nsamples];
  Int_t _padAmpl_saclay[geom::nPadx][geom::nPady][geom::Nsamples_saclay];

  /// output vector to put in the file
  std::vector<TObject*> _output_vector;

  /// Reconstruction used in the analysis.
  /**  You can use plenty in the analysis. At least one should be defines */
  ReconstructionBase* _reconstruction;

  /// Selection parameters
  /// The maximum multiplicity of the track
  Int_t _max_mult;

  /// The maximum mean multiplicity for track
  Float_t _max_mean_mult;

  /// Whether to cut tracks with gap in the cluster
  /** E.g. one missed pad in the middle of the column or in the
  * middle of the diagonal
  */
  bool _cut_gap;

  /// Minimum number of clusters in the track
  Int_t _min_clusters;

  /// Maximum angle (abs(tan)) in the MM plane
  Float_t _max_phi;

  /// Maximum angle (abs(tan)) w.r.t. MM plane
  Float_t _max_theta;

  /// Hot to treat cross-talk
  enum cross_talk {
    def = 0,
    suppress,
    cherry_pick
  };

  cross_talk _cross_talk_treat;

  /// T2K plotting style
  TStyle* _t2kstyle;

  /// DEBUG vars
  Int_t _verbose;
  bool _batch;
  bool _test_mode;
  bool _overwrite;

  /// Whether to invert track analysis logic
  /** E.g. analyse cosmic tracks. Rows and columns will be replaced */
  bool _invert;

  /// Whether to use Gaussian lorentzian PRf fit over polynomial
  bool _gaus_lorentz_PRF;

  /// Wheather to use individual PRF
  bool _individual_column_PRF;

  /// Wheather to make PRF center position a free parameter
  bool _PRF_free_centre;

  /// Whether to use arc function for track fitting
  bool _do_linear_fit;
  bool _do_para_fit;

  /// Wether to store the WFs
  bool _to_store_wf;

  /// Time control system
  TApplication* _app;
  TStopwatch* _sw_event;

  TStopwatch* _sw_partial[6];
};


#endif  // SRC_BASE_ANALYSISBASE_HXX_
