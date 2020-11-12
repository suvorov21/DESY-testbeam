#!/usr/bin/env python3

import argparse
import subprocess
import shlex
import shutil
import os
import random
from sys import exit
from itertools import chain

def main():
    # **********************************************************************************
    #define input
    n_iter = 20
    doiter = True

    generate_TEvent = True
    submit = False
    launch = True

    bin_dir   = "/afs/cern.ch/work/s/ssuvorov/public/T2K_testbeam/DESY_TestBeam/build_911/"
    bin_name  = "SpatialResol.exe"
    bin_flag  = "-b -a"

    input_prefix  = ["/eos/experiment/neutplatform/t2knd280/DESY_TPC/ROOT/v1/",
                     "/eos/experiment/neutplatform/t2knd280/ERAM_2/",
                     "/eos/experiment/neutplatform/t2knd280/ERAM2-data/",
                     "/eos/experiment/neutplatform/t2knd280/ERAM_3/"]
    # input_prefix  = ["/eos/user/s/ssuvorov/DESY_testbeam/"]
    input_version = ""

    outpt_prefix  = "/eos/user/s/ssuvorov/DESY_testbeam/"
    outpt_version = "cosmic_2020"
    output_post   = ""

    # espresso     = 20 minutes
    # microcentury = 1 hour
    # longlunch    = 2 hours
    # workday      = 8 hours
    # tomorrow     = 1 day
    # testmatch    = 3 days
    # nextweek     = 1 week
    job_flavour = "longlunch"
    log_folder = "/afs/cern.ch/work/s/ssuvorov/public/T2K_testbeam/DESY_TestBeam/log/"
    # end of input definition
    # **********************************************************************************

    if (submit and launch):
        print("Submit and launch should not be executed together")
        return 0

    temp = "/temp_0/"

    parser = argparse.ArgumentParser(description='Submit jobs to condor at LXPLUS')
    parser.add_argument("-f", metavar="f", type=str,
        help='File list with input/output files names',
        required=True)

    args = parser.parse_args()

    project_path = os.path.abspath(os.path.dirname(os.path.abspath(__file__)) + "/../")
    if os.path.exists(project_path + "/script/" + temp):
        print("Rewriting the temp folder ", project_path + "/script/" + temp)
        shutil.rmtree(project_path + "/script/" + temp)
        os.mkdir(project_path + "/script/" + temp)
    else:
        print("Creating the temp folder ", project_path + "/script/" + temp)
        os.mkdir(project_path + "/script/" + temp)

    if not os.path.exists(outpt_prefix + outpt_version):
        print("Creating the output folder")
        os.mkdir(outpt_prefix + outpt_version)


    print("Creating tasks")
    with open(args.f) as file_list:
        i = 0
        # for each input file
        if launch:
            launcher = open(project_path + "/script/" + temp + "all.sh", "w")
        else:
            launcher = open(project_path + "/script/" + temp + "all.sh.bu", "w")
        for line in file_list:
            # skip line if comment-like
            line = line.split('#')[0]
            if line.isspace() or line == '':
                continue

            in_file  = line.split()[0]
            ot_file  = line.split()[1]

            if (in_file == "" or ot_file == ""):
                print("ERROR. No file names specified")
                exit(-1)

            # create a file list in case of existing subruns
            temp_filename = project_path + "/script/" + temp
            temp_filename += str(round(random.random()*1000)) + ".list"
            temp_file = open(temp_filename, "w")
            first_file_name = ""

            # path = input_prefix+"/"+input_version+"/"

            for root, dirt, find_file in chain.from_iterable(os.walk(path) for path in input_prefix):
                if "soft" in root:
                    continue
                for file in find_file:
                    if in_file[:21] in file and ".root" in file:
                        temp_file.write(os.path.join(root, file) + "\n")
                        if first_file_name == "":
                            first_file_name = file

            temp_file.close()

            print(first_file_name, "-->", ot_file)

            file_out = open(project_path + "/script/" + temp + str(i) + ".sh", "w")
            command = ""

            # fo each iteration
            for i_iter in range(0, n_iter):
                command += bin_dir + "/" + bin_name + " " + bin_flag + " -t " + str(i_iter)
                if (i_iter == 0 and generate_TEvent):
                    command += " -s "
                if (not generate_TEvent or i_iter == 0):
                    command += " -i " + temp_filename
                else:
                    command += " -i " + outpt_prefix+"/"+outpt_version+"/"+first_file_name

                command += " -o " + outpt_prefix+"/"+outpt_version+"/"+ot_file + output_post
                if doiter:
                    command += "_iter" + str(i_iter)
                command += ".root; "
                if not doiter:
                    break

            file_out.write("#!/bin/bash\n")
            # file_out.write("source /cvmfs/sft.cern.ch/lcg/contrib/gcc/7.3.0binutils/x86_64-centos7-gcc7-opt/setup.sh\n")
            # file_out.write("source /afs/cern.ch/work/s/ssuvorov/public/ROOT/root-6.18.00-build_gcc73/bin/thisroot.sh\n")
            file_out.write("cd " + bin_dir + "\n")
            file_out.write(command + "\n")

            file_out.close()

            launcher.write("chmod 765 ./" + str(i) + ".sh\n")
            launcher.write("./" + str(i) + ".sh\n")
            i+=1

        launcher.close()

    submit_file = open(project_path + "/script/" + temp + "/Submit.sub", "w")

    submit_file.write("executable              = $(filename)\n")

    submit_file.write("arguments               = $(ClusterId)$(ProcId)\n")
    submit_file.write("requirements            = (OpSysAndVer =?= \"CentOS7\")\n")
    submit_file.write("+job_flavour             = \"" + job_flavour + "\"\n")
    submit_file.write("output                  = " + log_folder + "/output/$(ClusterId).$(ProcId).out\n")
    submit_file.write("error                   = " + log_folder + "/error/$(ClusterId).$(ProcId).err\n")
    submit_file.write("log                     = " + log_folder + "/log/$(ClusterId).log\n")
    submit_file.write("queue filename matching files *.sh")

    submit_file.close()

    os.chdir(project_path + "/script/" + temp)
    if launch:
        subprocess.run(["chmod", "765", "./all.sh"], check=True)
        subprocess.run(["/bin/bash", "all.sh"], check=True)
    if submit:
        subprocess.run(["condor_submit",  "Submit.sub"], check=True)
    os.chdir(project_path + "/script/")
    return 0


if __name__ == "__main__":
    main()
