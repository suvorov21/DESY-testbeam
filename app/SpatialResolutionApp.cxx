#include "SpatialResolAna.hxx"

int main(int argc, char** argv) {
  auto ana = new SpatialResolAna();
  if (!ana->Initialize(argc, argv))               return -1;
  if (!ana->Loop(ana->GetEventList()))  return -1;
  if (!ana->WriteOutput())              return -1;

  return 0;
}