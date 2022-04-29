#include "SpatialResolAna.hxx"

int main(int argc, char** argv) {
  auto ana = std::make_unique<SpatialResolAna>();
  if (!ana->ReadCLI(argc, argv))               return -1;
  if (!ana->Initialize())               return -1;
  if (!ana->Loop())  return -1;
  if (!ana->WriteOutput())              return -1;

  return 0;
}