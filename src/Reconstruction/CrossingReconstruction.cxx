#include "CrossingReconstruction.hxx"
#include "line.hxx"

CrossingReconstruction::CrossingReconstruction(): ReconstructionBase() {
  ;
}

bool CrossingReconstruction::Initialize() {
  std::cout << "Initialize crossing Reconstruction............";

  std::cout << "done" << std::endl;
  return true;
}

bool CrossingReconstruction::SelectEvent(const Int_t padAmpl[geom::nPadx][geom::nPady][geom::Nsamples],
                                TEvent* event) {

  std::vector<float> front_cluster_t;
  std::vector<float> front_cluster_y;
  std::vector<float> front_cluster_x;
  std::vector<int>   front_cluster_N;
  std::vector<int>   front_cluster_charge;

  // fill max per WF ampl
  int PadDisplay[geom::nPadx][geom::nPady]  = {0};
  int PadTime[geom::nPadx][geom::nPady]     = {0};
  for (auto it_x = 0; it_x < geom::nPadx; ++it_x)
    for (auto it_y = 0; it_y < geom::nPady; ++it_y)
      for (auto it_t = 0; it_t < geom::Nsamples; ++it_t)
        if (padAmpl[it_x][it_y][it_t] > PadDisplay[it_x][it_y]) {
          PadDisplay[it_x][it_y]  = padAmpl[it_x][it_y][it_t];
          PadTime[it_x][it_y]     = it_t;
        }

  // looping over all side pads to find time & space cluster
  for (int it_y = 0; it_y < geom::nPady; ++it_y) {
    for (int it_x = 0; it_x < cluster_offset; ++it_x) {
      if (PadDisplay[it_x][it_y] <= cluster_threshold)
        continue;

      int index = -1;
      for (uint it = 0; it < front_cluster_t.size(); ++it) {
        if (abs(PadTime[it_x][it_y] - front_cluster_t[it])  < cluster_time_gate
         && abs(it_y             - front_cluster_y[it]) < cluster_space_gate ) {
          index = it;
          break;
        }
      }

      if (index == -1) {
        front_cluster_t.push_back(1. * PadTime[it_x][it_y]);
        front_cluster_y.push_back(1. * it_y);
        front_cluster_x.push_back((cluster_offset - 1) / 2);
        front_cluster_N.push_back(1);
        front_cluster_charge.push_back(PadDisplay[it_x][it_y]);
      } else {
        // assign to existing cluster
        // average the space anf time
        front_cluster_t[index] = 1. * front_cluster_t[index] * front_cluster_N[index] + PadTime[it_x][it_y];
        front_cluster_y[index] = 1. * front_cluster_y[index] * front_cluster_N[index] + it_y;
        ++front_cluster_N[index];
        front_cluster_t[index] /= front_cluster_N[index];
        front_cluster_y[index] /= front_cluster_N[index];
        front_cluster_charge[index] += PadDisplay[it_x][it_y];
      }
    } // over it_x within cluster_offset
  } // over it_y

  if (front_cluster_t.size() > cluster_lim || !front_cluster_t.size())
    return false;

  std::vector<float> back_cluster_t;
  std::vector<float> back_cluster_y;
  std::vector<float> back_cluster_x;
  std::vector<int>   back_cluster_N;
  std::vector<int>   back_cluster_charge;

  for (int it_y = 0; it_y < geom::nPady; ++it_y) {
    for (int it_x = geom::nPadx - cluster_offset; it_x < geom::nPadx; ++it_x) {
      if (PadDisplay[it_x][it_y] <= cluster_threshold)
        continue;

      int index = -1;
      for (uint it = 0; it < back_cluster_t.size(); ++it) {
        if (abs(PadTime[it_x][it_y] - back_cluster_t[it])  < cluster_time_gate
         && abs(it_y             - back_cluster_y[it]) < cluster_space_gate ) {
          index = it;
          break;
        }
      }

      if (index == -1) {
        // create new cluster
        back_cluster_t.push_back(1. * PadTime[it_x][it_y]);
        back_cluster_y.push_back(1. * it_y);
        back_cluster_x.push_back(geom::nPadx - (cluster_offset - 1) / 2);
        back_cluster_N.push_back(1);
        back_cluster_charge.push_back(PadDisplay[it_x][it_y]);
      } else {
        // assign to existing cluster
        // average the space anf time
        back_cluster_t[index] = 1. * back_cluster_t[index] * back_cluster_N[index] + PadTime[it_x][it_y];
        back_cluster_y[index] = 1. * back_cluster_y[index] * back_cluster_N[index] + it_y;
        ++back_cluster_N[index];
        back_cluster_t[index] /= back_cluster_N[index];
        back_cluster_y[index] /= back_cluster_N[index];
        back_cluster_charge[index] += PadDisplay[it_x][it_y];
      }
    } // over it_x within cluster_offset
  } // over it_y

  if (back_cluster_t.size() > cluster_lim || !back_cluster_t.size())
    return false;

   // store the tracks info
  std::vector<std::vector<std::vector<short> > > track_container;
  std::vector<std::vector<std::vector<short> > > track_container_time;
  std::vector<int> track_col;

  bool pileup = false;

  // loop over all possible combinations
  for (uint trackStart = 0; trackStart < front_cluster_t.size(); ++trackStart) {
    for (uint trackEnd = 0; trackEnd < back_cluster_t.size(); ++trackEnd) {

      // create a draft for tarck
      std::vector<std::vector<short> > TrackDraft;
      std::vector<std::vector<short> > TrackDraft_time;
      TrackDraft.resize(geom::nPadx);
      TrackDraft_time.resize(geom::nPadx);
      for (int z=0; z < geom::nPadx; ++z) {
        TrackDraft[z].resize(geom::nPady, 0);
        TrackDraft_time[z].resize(geom::nPady, 0);
      }

      TVector3 vec1 = TVector3(front_cluster_x[trackStart], front_cluster_y[trackStart], front_cluster_t[trackStart]);
      TVector3 vec2 = TVector3(back_cluster_x[trackEnd], back_cluster_y[trackEnd], back_cluster_t[trackEnd]);
      TLine_att line(vec1, vec2);

      // scan the pad with the threshold and within a box
      // count number of columns
      int columns = 0;

      int breaking = 0;
      bool gap = false;

      // pads out of time but above threshold
      int pad_OOT = 0;

      for (int it_x = 0; it_x < geom::nPadx; ++it_x) {
        bool column = false;
        int track_max   = -1;
        int track_max_y = -1;
        int track_max_t = -1;
        bool columan_max_adc = false;
        (void)columan_max_adc;

        for (int it_y = 0; it_y < geom::nPady; ++it_y) {

          if (!PadDisplay[it_x][it_y])
            continue;

          TVector3 track_pos = line.EvalX(it_x);

          int found_max       = -1;
          int found_max_time  = -1;

          if (abs(track_pos.Z() - PadTime[it_x][it_y]) <= box_width_t1) {
            found_max       = PadDisplay[it_x][it_y];
            found_max_time  =  PadTime[it_x][it_y];
            columan_max_adc = true;
          } else {
            int first_time    = std::max((int)track_pos.Z() - scan_time_down, 0);
            int last_time     = std::min((int)track_pos.Z() + scan_time_up, geom::Nsamples);

            int it_t = first_time;

            while (it_t < geom::Nsamples) {
              if (found_max && !padAmpl[it_x][it_y][it_t])
                break;

              if (found_max == -1 && it_t > last_time)
                break;

              if (padAmpl[it_x][it_y][it_t] > found_max) {
                found_max       = padAmpl[it_x][it_y][it_t];
                found_max_time  = it_t;
              }

              ++it_t;
            } // loop over time
          }

          if (found_max == -1)
            continue;

          if (abs(track_pos.Y() - it_y) > box_width_j) {
            ++pad_OOT;
            continue;
          }

          if (found_max > track_max) {
            track_max   = found_max;
            track_max_y = it_y;
            track_max_t = found_max_time;
            column = true;
          }
        } // over it_y

        //if (!columan_max_adc)
         // pileup = true;

        if (column) {
          ++columns;
          TrackDraft[it_x][track_max_y]       = track_max;
          TrackDraft_time[it_x][track_max_y]  = track_max_t;
          if (gap) {
            ++breaking;
            gap = false;
          }
        } else if (columns)
          gap = true;
      } // over it_x

      // check number of columns and number of gaps
      if (breaking > breaking_thr || columns < column_thr || pad_OOT > OOT_thr)
        continue;

      if (pileup && pile_up_cut)
        continue;

      track_container.push_back(TrackDraft);
      track_container_time.push_back(TrackDraft_time);
      track_col.push_back(columns);
    } // track end
  } // track start

  if (!track_container.size())
    return false;

  std::vector<TTrack*> track_v;
  for (uint trackId = 0; trackId < track_container.size(); ++trackId) {
    for (int it_x = 0; it_x < geom::nPadx; ++it_x) {
      int it_y = 0;
      while (!track_container[trackId][it_x][it_y] && it_y < geom::nPady-1)
        ++it_y;

      if (it_y == geom::nPady-1)
        continue;

      int it_y_track = it_y;

      int first_time    = std::max(PadTime[it_x][it_y_track] - scan_time_down, 0);
      int last_time     = std::min(PadTime[it_x][it_y_track] + scan_time_up, geom::Nsamples);

      while (PadDisplay[it_x][it_y] && it_y < geom::nPady-1
        && abs(it_y - it_y_track) < track_collection_dist) {

        // if the maximum is inside the time box --> take it
        if(PadTime[it_x][it_y] - PadTime[it_x][it_y_track] < box_width_t2 && PadTime[it_x][it_y_track] - PadTime[it_x][it_y] < box_width_t1) {
          track_container[trackId][it_x][it_y]      = PadDisplay[it_x][it_y];
          track_container_time[trackId][it_x][it_y] = PadTime[it_x][it_y];
        // if not --> scan ADC vs. time for this pad
        } else {
          int adc_max_inbox = -1;
          int adc_max_t     = -1;

          int it_t = first_time;
          while (it_t < geom::Nsamples) {
            if (adc_max_inbox && !padAmpl[it_x][it_y][it_t])
              break;

            if (adc_max_inbox == -1 && it_t > last_time)
               break;

            if (padAmpl[it_x][it_y][it_t] > adc_max_inbox) {
               adc_max_inbox       = padAmpl[it_x][it_y][it_t];
               adc_max_t  = it_t;
             }
             ++it_t;
          } // loop over time

          // if the hit in the box is befor the max hit --> OK
          // if not --> check that we are taking the peak not the tail
          if (adc_max_t > PadTime[it_x][it_y] && adc_max_t - PadTime[it_x][it_y] < scan_time_up && pile_up_cut) {
            pileup = true;
            break;
          }
          if (adc_max_inbox != -1) {
            track_container[trackId][it_x][it_y]      = adc_max_inbox;
            track_container_time[trackId][it_x][it_y] = adc_max_t;
          } else {
            if (PadDisplay[it_x][it_y] > cluster_threshold && pile_up_cut) {
              //pileup = true;
            }
            break;
          }
        } // look for the hit in time window
        ++it_y;
      } // go up

      it_y = it_y_track;
      while (PadDisplay[it_x][it_y] && it_y > 0
        && abs(it_y - it_y_track) < track_collection_dist) {
        // if the maximum is inside the time box --> take it
        if(PadTime[it_x][it_y] - PadTime[it_x][it_y_track] < box_width_t2 && PadTime[it_x][it_y_track] - PadTime[it_x][it_y] < box_width_t1) {
          track_container[trackId][it_x][it_y] = PadDisplay[it_x][it_y];
          track_container_time[trackId][it_x][it_y] = PadTime[it_x][it_y];
        // if not --> scan ADC vs. time for this pad
        } else {
          int adc_max_inbox = -1;
          int adc_max_t     = -1;


          int it_t = first_time;
          while (it_t < geom::Nsamples) {
            if (adc_max_inbox && !padAmpl[it_x][it_y][it_t])
              break;

            if (adc_max_inbox == -1 && it_t > last_time)
               break;

            if (padAmpl[it_x][it_y][it_t] > adc_max_inbox) {
               adc_max_inbox       = padAmpl[it_x][it_y][it_t];
               adc_max_t  = it_t;
             }
             ++it_t;
           } // loop over time

          // if the hit in the box is befor the max hit --> OK
          // if not --> check that we are taking the peak not the tail
          if (adc_max_t > PadTime[it_x][it_y] && adc_max_t - PadTime[it_x][it_y] < scan_time_up && pile_up_cut) {
            pileup = true;
            break;
          }

          if (adc_max_inbox != -1) {
            track_container[trackId][it_x][it_y] = adc_max_inbox;
            track_container_time[trackId][it_x][it_y] = adc_max_t;
          } else {
            if (PadDisplay[it_x][it_y] > cluster_threshold && pile_up_cut) {
              //pileup = true;
            }

            break;
          }
        } // look for the hit in time window
        --it_y;
      } // go down (it_y)
    if (pileup && pile_up_cut)
      break;
    } // it_x
    if (pileup && pile_up_cut)
      continue;

    TTrack* track = new TTrack();
    std::vector<THit*> hits_v;
    for (auto x = 0; x < geom::nPadx; ++x) {
      for (auto y = 0; y < geom::nPady; ++y) {
        if (track_container[trackId][x][y]) {
          THit* hit = new THit();
          hit->SetQ(track_container[trackId][x][y]);
          hit->SetCol(x);
          hit->SetRow(y);
          hit->SetTime(track_container_time[trackId][x][y]);

          track->AddColHit(hit);
          track->AddRowHit(hit);
          hits_v.push_back(hit);
        } // not empty
      } // y
    }// x
    track->SetHits(hits_v);
    track_v.push_back(track);

  } // loop over tracks

  event->SetTracks(track_v);

  return true;
}

