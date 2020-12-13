#!/usr/bin/env python3

import os
from subprocess import run, Popen, PIPE

def main():
    # input definition
    com_path = 'build'
    command = 'SpatialResol.exe'
    flags = ['-b', '-r']

    multithread = False

    TEvent = True
    Niter = 10

    split = 10

    input_path  = "/eos/experiment/neutplatform/t2knd280/DESY_TPC/ROOT/v1/2019_06_14/"
    # input_path  = "/eos/user/s/ssuvorov/DESY_testbeam/tree/"
    # input_path  = '~/DATA/'
    input_name  = "/R2019_06_14-16_42_26-000.root"

    out_path    = "/eos/user/s/ssuvorov/DESY_testbeam/tree/"
    # out_path    = '~/DATA/'
    out_name    = "g_360_200_95p_nbr"

    ####################################################################
    ########### end of input definition ################################
    ####################################################################

    project_path = os.path.abspath(os.path.dirname(
                                                   os.path.abspath(__file__)
                                                   ) + '/../' + com_path)
    command = project_path + '/' + command
    comm = [command]
    for flag in flags:
        comm.append(flag)

    if multithread:
        import ROOT
        input_file = input_path+input_name
        f = ROOT.TFile(input_file, 'READ')
        try:
            Nentries = f.tree.GetEntries()
        except AttributeError:
            try:
                Nentries = f.event_tree.GetEntries()
            except:
                print('Unknown tree')
                return

    for iterId in range (0, Niter):
        exe = comm.copy()
        exe.append('-t'+str(iterId))
        # generate TEvent for the first time
        input_file = input_path+input_name
        if TEvent:
            if iterId == 0:
                exe.append('-s')
            else:
                input_file = out_path+input_name

        exe.append('-i ' + input_file)
        output_base = out_path+out_name+'_iter'+str(iterId)

        if not multithread:
            ####### single thread run ######################################
            exe.append('-o '+ output_base +'.root')
            print(' '.join(exe))
            run(exe, check=True)
            # ##### end of run  ############################################
        else:
            ####### multiple thread run  ###################################
            exe_list = []
            step = int(Nentries // split + 1)
            start = 0
            print('#######################################################')
            print('#  Multithread run may provide unexpected output. WIP #')
            print("#              Proceed with caution!!!                #")
            print('#######################################################')
            for thread in range(split):
                exe_step = exe.copy()
                exe_step.extend(['--start', str(start),
                                '--end', str(start+step),
                                '-o ' + output_base + '_part' + str(thread) + '.root'
                                ])
                exe_list.append(exe_step)
                start += step

            for exe_step in exe_list:
                print(' '.join(exe_step))

            procs_list = [Popen(cmd,
                                stdout=PIPE,
                                stderr=PIPE
                                ) for cmd in exe_list]
            for proc_id, proc in enumerate(procs_list):
                proc.wait()
                print(f'Done part {proc_id}')
            force = '-f ' if '-r' in flags else ''
            proc = Popen('hadd ' + force + output_base +'.root ' + output_base + '_part*.root',
                         shell=True,
                         stdout=PIPE,
                         stderr=PIPE
                         )
            proc.wait()
            print(f'Merge done iteration {iterId}')
            ####### end of run  ############################################


if __name__ == "__main__":
    main()
