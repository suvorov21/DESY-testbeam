#include "ReconstructionBase.hxx"

bool ReconstructionBase::Initialize(int verbose) {
    (void) verbose;
    std::cout << "WARNING. The default Reconstruction is initialised. The result is always true" << std::endl;
    return true;
}

bool ReconstructionBase::ReconstructEvent(const std::shared_ptr<TEvent> &event) {
    (void) event;

    return true;
}
