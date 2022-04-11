#include "TCluster.hxx"
#include "Geom.hxx"

//// TCLUSTER
TCluster::TCluster() : _x(-999), _y(-999), _y_error(-999), _charge(0) {}

TCluster::TCluster(const std::shared_ptr<THit> &pad) {
    _hits.push_back(pad);
    _x = (float_t) geom::GetXposPad(pad);
    _y = (float_t) geom::GetYposPad(pad);
    _y_error = 0;
    _charge = 0;
}

typename std::vector<std::shared_ptr<THit>>::iterator TCluster::begin() {
    return _hits.begin();
}

typename std::vector<std::shared_ptr<THit>>::iterator TCluster::end() {
    return _hits.end();
}

std::shared_ptr<THit> TCluster::operator[](size_t index) {
    return _hits[index];
}

Clustering::Clustering(clusterType type, bool invert) :
_invert(invert), _type(type) {
    switch (type) {
        case clusterType::kRowColumn:
            _angle = 0.;
            _n_pads = 0.;
            break;
        case clusterType::kDiagonal:
            _angle = units::a45;
            _n_pads = 1;
            break;
        case clusterType::k2by1:
            _angle = invert ? units::a2_inv : units::a2;
            _n_pads = 2;
            break;
        case clusterType::k3by1:
            _angle = invert ? units::a3_inv : units::a3;
            _n_pads = 3;
            break;
        default: std::cerr << "Clustering is not defined!" << std::endl;
    }
    _coeff = (_n_pads == 0) ? 0.001 : 1. / _n_pads;

}

//******************************************************************************
TClusterPtrVec Clustering::ClusterTrack(const THitPtrVec &tr) const {
//******************************************************************************
    TClusterPtrVec cluster_v;
    for (const auto &pad : tr) {
        auto col_id = pad->GetCol(_invert);
        auto row_id = pad->GetRow(_invert);

        // skip first and last row/column
        if (row_id == 0 || row_id == geom::GetNRow(_invert) - 1 ||
            col_id == 0 || col_id == geom::GetNColumn(_invert) - 1)
            continue;

        auto cons = GetConstant(row_id, col_id);

        // search if the cluster is already considered
        TClusterPtrVec::iterator it;
        for (it = cluster_v.begin(); it < cluster_v.end(); ++it) {
            if (!(**it)[0]) {
                continue;
            }

            auto cluster_col = (**it)[0]->GetCol(_invert);
            auto cluster_row = (**it)[0]->GetRow(_invert);
            if (GetConstant(cluster_row, cluster_col) == cons) {
                (*it)->AddHit(pad);
                (*it)->AddCharge(pad->GetQ());
                /** update X position */
                auto x_pad = geom::GetXposPad(pad, _invert, _angle);
                auto mult = (*it)->GetSize();
                auto x_new = ((*it)->GetX() * ((Float_t) mult - 1) + x_pad) / num::cast<double>(mult);
                (*it)->SetX((float_t) x_new);
                /** */

                break;
            }
        } // loop over track clusters
        // add new cluster
        if (it == cluster_v.end()) {
            auto first_cluster = std::make_unique<TCluster>(pad);
            first_cluster->SetX((float_t) geom::GetXposPad(pad, _invert, _angle));
            first_cluster->SetCharge(pad->GetQ());
            cluster_v.push_back(std::move(first_cluster));
        }
    } // over pads

    return cluster_v;
}

void Clustering::setInvert(bool invert) {
    _invert = invert;
    if (_type == clusterType::k3by1) {
        _angle = invert ? units::a3_inv : units::a3;
    }
    if (_type == clusterType::k2by1) {
        _angle = invert ? units::a2_inv : units::a2;
    }
}