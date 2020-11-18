#include "InclinedTracks.hxx"

int main(int argc, char** argv) {
  auto ana = new InclinedTracks(argc, argv);
  if (!ana->Initialize())               return -1;
  if (!ana->Loop(ana->GetEventList()))  return -1;
  if (!ana->WriteOutput())              return -1;

  return 0;
}