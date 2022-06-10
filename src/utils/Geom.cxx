#include "Geom.hxx"
#include <fstream>

Geom::Geom() {
    auto source = std::string(__FILE__);
    auto found = source.find_last_of('/');
    std::string geom_file_name = source.substr(0, found) + "/geom.dat";

    std::ifstream file(geom_file_name);
    if (file.is_open()) {
        std::string line;
        getline(file, line);
        int x, y, chip, channel;
        while (!file.eof()) {
            file >> x >> y >> chip >> channel;
            padToEle[x * nPady + y] = std::make_pair(chip, channel);
            eleToPad[chip * nChannel + channel] = std::make_pair(x, y);
        }
    }
}

Geom &Geom::get() {
    static Geom g;
    return g;
}

std::pair<short, short> Geom::GetPadToEle(const int x, const int y) {
    return padToEle[x * nPady + y];
}

std::pair<short, short> Geom::GetEleToPad(const int chip, const int channel) {
    return eleToPad[chip * nChannel + channel];
}

Double_t Geom::GetXposPad(const THitPtr &h, bool invert, Double_t angle) {
    Double_t x_offset{0};
    Double_t y_offset{0};
    if (invert) {
        x_offset = -1 * num::cast<Double_t>(h->GetCard() / 4 - 1) * yModuleOffset;
        y_offset = num::cast<Double_t>(h->GetCard() % 4) * xModuleOffset;
    } else {
        x_offset = num::cast<Double_t>(h->GetCard() % 4) * xModuleOffset;
        y_offset = -1 * num::cast<Double_t>(h->GetCard() / 4 - 1) * yModuleOffset;
    }

    Double_t x_flat = Geom::GetXpos(h, invert) + x_offset;
    Double_t y_flat = Geom::GetYpos(h, invert) + y_offset;
    return x_flat * TMath::Cos(angle) -
        y_flat * TMath::Sin(angle);
}

Double_t Geom::GetYposPad(const THitPtr &h, bool invert, Double_t angle) {
    Double_t x_offset{0};
    Double_t y_offset{0};
    if (invert) {
        x_offset = -1 * num::cast<Double_t>(h->GetCard() / 4 - 1) * yModuleOffset;
        y_offset = num::cast<Double_t>(h->GetCard() % 4) * xModuleOffset;
    } else {
        x_offset = num::cast<Double_t>(h->GetCard() % 4) * xModuleOffset;
        y_offset = -1 * num::cast<Double_t>(h->GetCard() / 4 - 1) * yModuleOffset;
    }
    Double_t x_flat = Geom::GetXpos(h, invert) + x_offset;
    Double_t y_flat = Geom::GetYpos(h, invert) + y_offset;
    return x_flat * TMath::Sin(angle) +
        y_flat * TMath::Cos(angle);
}

Double_t Geom::GetXpos(const THitPtr &h, bool invert) {
    return Geom::GetXpos(h->GetCol(invert), invert);
}

Double_t Geom::GetYpos(const THitPtr &h, bool invert) {
    return Geom::GetYpos(h->GetRow(invert), invert);
}

Double_t Geom::GetXpos(int it_x, bool invert) {
    if ((!invert && it_x >= Geom::nPadx) ||
        (invert && it_x >= Geom::nPady) || it_x < 0) {
        TString msg = __func__;
        msg += "Wrong Index " + std::to_string(it_x) + "\t" + std::to_string(invert);
        throw std::logic_error(msg);
    }
    if (!invert)
        return Geom::GetXraw(it_x);
    else
        return Geom::GetYraw(it_x);
}

Double_t Geom::GetYpos(int it_y, bool invert) {
    if ((!invert && it_y >= Geom::nPady) ||
        (invert && it_y >= Geom::nPadx) || it_y < 0) {
        TString msg = __func__;
        msg += "Wrong Index " + std::to_string(it_y) + "\t" + std::to_string(invert);
        throw std::logic_error(msg);
    }
    if (!invert)
        return Geom::GetYraw(it_y);
    else
        return Geom::GetXraw(it_y);
}

int Geom::GetNColumn(bool invert) {
    if (!invert)
        return nPadx;
    else
        return nPady;
}

int Geom::GetNRow(bool invert) {
    return Geom::GetNColumn(!invert);
}

constexpr Double_t Geom::GetXraw(int padXid) {
    return dx * (padXid - 1. * (nPadx - 1) / 2);
}

constexpr Double_t Geom::GetYraw(int padYid) {
    return dy * (padYid - 1. * (nPady - 1) / 2);
}
