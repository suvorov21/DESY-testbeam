#!/usr/bin/env python3

import argparse
import subprocess
import shutil
import os
import random
from sys import exit
from itertools import chain

def main():
    # **********************************************************************************
    #define input
    n_iter = 0
    doiter = True
    #max_subruns = 3

    generate_TEvent = False
    submit = False
    launch = True

    bin_dir   = "/afs/cern.ch/work/u/uyevarou/private/Desy21/dev_v2/desy_testbeam/build"
    bin_name  = "SpatialResol.exe"
    bin_flag  = "-b -d"

    input_prefix  = ["/eos/experiment/neutplatform/t2knd280/DESY_TPC/ROOT/v1/",
                    "/eos/experiment/neutplatform/t2knd280/DESY_TPC/ROOT/v1/2019_06_16/",
                    "/eos/experiment/neutplatform/t2knd280/DESY_TPC/ROOT/v1/2019_06_14/",
                    "/eos/experiment/neutplatform/t2knd280/DESY_TPC/ROOT/v1/2019_06_15/",
                     # "/eos/experiment/neutplatform/t2knd280/ERAM_2/",
                     # "/eos/experiment/neutplatform/t2knd280/MM2_data/",
                     # "/eos/experiment/neutplatform/t2knd280/ERAM2-data/",
                     # "/eos/experiment/neutplatform/t2knd280/ERAM_3/",
                     # "/afs/cern.ch/work/s/shassani/public/Cosmics/",
                     # "/eos/experiment/neutplatform/t2knd280/HATPC-Cosmics-2021/ERAM1/",
                     "/eos/experiment/neutplatform/t2knd280/DESY-2021/ROOT/v1/",
                     "/eos/experiment/neutplatform/t2knd280/CERN-2021/ROOT/",
                     "/eos/experiment/neutplatform/np07/HAT/CERN22_MockUp/",
                     # "/eos/user/s/ssuvorov/DESY_testbeam/hs/"
                     ]

    #outpt_prefix  = "/eos/user/u/uyevarou/desy_testbeam/"
    outpt_prefix  = "/eos/experiment/neutplatform/np07/HAT/DESY_2021/vlada_reco/new/"
    #outpt_version = "z_phiFull_dEdx_fin"
    #outpt_version = "z_phi_SR_fin_parabola_Lenth"
    #outpt_version = "dev_v16_highStat"
    outpt_version = "CERN22_output_v4"
    #outpt_version = "DESY_2019"
    #outpt_version = "z_phi_Theta_add"
    #outpt_version = "z_phi_new_dEdx_add" 
    output_post   = ""

    # espresso     = 20 minutes
    # microcentury = 1 hour
    # longlunch    = 2 hours
    # workday      = 8 hours
    # tomorrow     = 1 day
    # testmatch    = 3 days
    # nextweek     = 1 week
    #job_flavour = "nextweek"
    job_flavour = "tomorrow"
    log_folder = "/afs/cern.ch/work/u/uyevarou/public/Desy21/log"
    # end of input definition
    # **********************************************************************************
    
    if (submit and launch):
        print("Submit and launch should not be executed together")
        return 0

    temp = "/temp/"

    parser = argparse.ArgumentParser(description='Submit jobs to condor at LXPLUS')
    parser.add_argument("-f", metavar="f", type=str,
        help='File list with input/output files names',
        required=True)

    args = parser.parse_args()

    project_path = os.path.abspath(os.path.dirname(os.path.abspath(__file__)) + "/../")
    script_path = project_path + "/script/" + temp
    if os.path.exists(script_path):
        resp = input('Temp folder with this name exists. Rewrite: (y/n): ')
        if resp != 'y':
            return 1
        print("Rewriting the temp folder ", script_path)
        shutil.rmtree(script_path)
        os.mkdir(script_path)
    else:
        print("Creating the temp folder ", script_path)
        os.mkdir(script_path)

    if not os.path.exists(outpt_prefix + outpt_version):
        print("Creating the output folder")
        os.mkdir(outpt_prefix + outpt_version)


    print("Creating tasks")
    with open(args.f) as file_list:
        i = 0
        # for each input file
        # create a bash script
        if launch:
            launcher = open(script_path + "all.sh", "w")
        else:
            launcher = open(script_path + "all.sh.bu", "w")
        for line in file_list:
            # skip line if comment-like
            line = line.split('#')[0]
            if line.isspace() or line == '':
                continue

            # parse inoput/output file names and additional flags
            in_file  = line.split()[0]
            ot_file  = line.split()[1]
            add_flags = []
            add_flags.extend(line.split()[2:])

            if (in_file == "" or ot_file == ""):
                print("ERROR. No file names specified")
                exit(-1)

            # create a file list in case of existing subruns
            temp_filename = script_path
            temp_filename += str(round(random.random()*1000)) + ".list"
            temp_file = open(temp_filename, "w")
            first_file_name = ""

            # search for inout files
            for root, dirt, find_file in chain.from_iterable(os.walk(path) for path in input_prefix):
                if "soft" in root:
                    continue
                for file in find_file:
                    if in_file[:21] in file and ".root" in file[-5:]:
                        temp_file.write(os.path.join(root, file) + "\n")
                        if first_file_name == "":
                            first_file_name = file

            temp_file.close()

            print(first_file_name, "--> ", f'{outpt_version}/{ot_file}', " with ", "".join(add_flags))

            file_out = open(script_path + str(i) + ".sh", "w")
            command = ""

            # fo each iteration
            for i_iter in range(0, n_iter):
                command += bin_dir + "/" + bin_name + " " + bin_flag
                command += " " + " ".join(add_flags)
                command += " -t " + str(i_iter)

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

            # create a bash file with a job
            file_out.write("#!/bin/bash\n")
            file_out.write("cd " + bin_dir + "\n")
            if submit:
                file_out.write("export PATH=/cvmfs/sft.cern.ch/lcg/contrib/CMake/3.18.3/Linux-x86_64/bin/:${PATH}\n")
                file_out.write(". /cvmfs/sft.cern.ch/lcg/views/LCG_100/x86_64-centos7-clang11-opt/setup.sh\n")
            file_out.write(command + "\n")

            file_out.close()

            launcher.write("chmod 765 ./" + str(i) + ".sh\n")
            launcher.write("./" + str(i) + ".sh\n")
            # launcher.write("email " + "\"Job done\" " + "\"ot_file\"\n")
            i+=1

        launcher.close()

    # create Condor submission file
    submit_file = open(script_path + "/Submit.sub", "w")

    submit_file.write("executable              = $(filename)\n")
    #submit_file.write("RequestCpus              = 4\n")
    submit_file.write("arguments               = $(ClusterId)$(ProcId)\n")
    submit_file.write("requirements            = (OpSysAndVer =?= \"CentOS7\")\n")
    submit_file.write("+JobFlavour             = \"" + job_flavour + "\"\n")
    submit_file.write("output                  = " + log_folder + "/output/$(ClusterId).$(ProcId).out\n")
    submit_file.write("error                   = " + log_folder + "/error/$(ClusterId).$(ProcId).err\n")
    submit_file.write("log                     = " + log_folder + "/log/$(ClusterId).log\n")
    submit_file.write("queue filename matching files *.sh")

    submit_file.close()

    os.chdir(script_path)
    if launch:
        subprocess.run(["chmod", "765", "./all.sh"], check=True)
        subprocess.run(["/bin/bash", "-c",
                       f". ~/.bashrc; cd {script_path}; . all.sh"
                       ],
                       check=True
                       )
    if submit:
        subprocess.run(["condor_submit",  "Submit.sub"], check=True)
    os.chdir(project_path + "/script/")
    return 0


if __name__ == "__main__":
    main()
