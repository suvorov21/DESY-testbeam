#!/usr/bin/env python3

import os
import random
from sys import exit
import subprocess

def main():
    # input definition
    command = 'SpatialResol.exe.test'
    flags = ['-b', '-r']

    TEvent = True
    Niter = 20

    input_path  = "/eos/experiment/neutplatform/t2knd280/DESY_TPC/ROOT/v1/2019_06_14/"
    input_name  = "R2019_06_14-17_34_58-000.root"

    out_path    = "/eos/user/s/ssuvorov/DESY_testbeam/"
    out_name    = "low_gain"
    # end of input definition

    project_path = os.path.abspath(os.path.dirname(os.path.abspath(__file__)) + '/../bin/')
    #command = 'bin/'+command
    command = project_path + '/' + command
    comm = [command]
    for flag in flags:
        comm.append(flag)

    for i in range (0, Niter):
        exe = comm.copy()
        exe.append('-t'+str(i))
        # generate TEvent for the first time
        if TEvent:
            if i == 0:
                exe.append('-s')
                exe.append('-i '+input_path+input_name)
            else:
                exe.append('-i '+out_path+input_name)
        else:
            exe.append('-i '+input_path+input_name)

        # add iteration number to output file
        exe.append('-o '+out_path+out_name+'_iter'+str(i)+'.root')

        print(exe)
        subprocess.run(exe, check=True)


if __name__ == "__main__":
    main()
