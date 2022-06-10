/** @cond */
#include <utility>
#include <memory>

#include "TString.h"
#include "TChain.h"
/** @endcond */

#include "TEvent.hxx"
#include "Geom.hxx"

//
// Created by SERGEY SUVOROV on 06/04/2022.
//

#ifndef DESY_TESTBEAM_INTERFACE_HXX
#define DESY_TESTBEAM_INTERFACE_HXX


enum class interfaceType {
    kArray510,
    kArray511,
    kRawEvent
};

class Interface {
 protected:
    /// input file name
    TString _file_in_name{""};

    /// chain with input files
    std::unique_ptr<TChain> _chain{nullptr};
 public:
    explicit Interface(TString var) : _file_in_name(std::move(var)) {}
    virtual ~Interface() = default;

    /// extract the name of the input ROOT file
    /// in case of list input take the first ROOT file
    static TString getRootFile(const TString&);
    static interfaceType getFileType(const TString&);

    /// chain input files
    bool chainInputFiles(const TString& treeName);

    Long64_t getEntries();

    virtual void Initialize() = 0;
    virtual std::shared_ptr<TEvent> getEvent(Int_t) = 0;

};

template<short timeSize = 511>
class interfaceRoot : public Interface {
    /// array of WFs
    Int_t _padAmpl[Geom::nPadx][Geom::nPady][timeSize]{-260};
 public:
    explicit interfaceRoot(TString var) : Interface(std::move(var)) {}
    ~interfaceRoot() override = default;
    void Initialize() override;
    std::shared_ptr<TEvent> getEvent(Int_t) override;
};

class interfaceRawEvent : public Interface {
 private:
    TRawEvent* _event{nullptr};
 public:
    explicit interfaceRawEvent(TString var) : Interface(std::move(var)) {}
    void Initialize() override;
    std::shared_ptr<TEvent> getEvent(Int_t) override;
};

/// Factory function that builds an interface
 class interfaceFactory {public:
    static std::unique_ptr<Interface> get(const TString &filename,
                                                interfaceType type) {
        switch (type) {
            case interfaceType::kArray510 :
                return std::make_unique<interfaceRoot<510>>(filename);
            case interfaceType::kArray511:
                return std::make_unique<interfaceRoot<511>>(filename);
            case interfaceType::kRawEvent:
                return std::make_unique<interfaceRawEvent>(filename);
            default: return nullptr;
        }
    }
};


#endif // DESY_TESTBEAM_INTERFACE_HXX
