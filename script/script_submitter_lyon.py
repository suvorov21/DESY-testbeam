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
    n_iter = 1
    doiter = True

    generate_TEvent = True
    submit = False
    launch = True

    #bin_dir   = "/sps/t2k/giganti/desy_testbeam/build/"
    bin_dir_script = os.getcwd()
    bin_dir = bin_dir_script + "/../build"
    #print(bin_dir)

    #bin_name  = "SpatialResol.exe"
    bin_name  = "SpatialResol.exe"
    bin_flag  = "-b -r -d"

    input_prefix  = [#"/eos/experiment/neutplatform/t2knd280/DESY_TPC/ROOT/v1/",
                     #"/eos/experiment/neutplatform/t2knd280/ERAM_2/",
                     #"/eos/experiment/neutplatform/t2knd280/MM2_data/",
                     #"/eos/experiment/neutplatform/t2knd280/ERAM2-data/",
                     #"/eos/experiment/neutplatform/t2knd280/ERAM_3/",
                     #"/sps/t2k/testbeamdata/2019_06/root/v1",
                     "/sps/t2k/testbeamdata/2019_06/root/v1/2019_06_16/",
                     "/sps/t2k/testbeamdata/2019_06/root/v1/2019_06_14/",
                     "/sps/t2k/testbeamdata/2019_06/root/v1/2019_06_15/",
                     ]
    # input_prefix  = ["/eos/user/s/ssuvorov/DESY_testbeam/"]
    input_version = ""

    outpt_prefix  = bin_dir_script + "/../"
    outpt_version = "output_files_st2_v6"
    outpt_log     = "/logs"

    output_post   = ""

    # espresso     = 20 minutes
    # microcentury = 1 hour
    # longlunch    = 2 hours
    # workday      = 8 hours
    # tomorrow     = 1 day
    # testmatch    = 3 days
    # nextweek     = 1 week
    job_flavour = "tomorrow"
    log_folder = bin_dir_script + "/../log"
    # end of input definition
    # **********************************************************************************

    if (submit and launch):
        print("Submit and launch should not be executed together")
        return 0

    temp = "/temp_1/"

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
    if not os.path.exists(outpt_prefix + outpt_version + outpt_log):
        print("Creating the log folder")
        os.mkdir(outpt_prefix + outpt_version + outpt_log)


    print("Creating tasks")
    with open(args.f) as file_list:
        i = 0
        # for each input file
        # create a bash script
        if launch:
            launcher = open(project_path + "/script/" + temp + "all.sh", "w")
        else:
            launcher = open(project_path + "/script/" + temp + "all.sh.bu", "w")
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
            temp_filename = project_path + "/script/" + temp
            temp_filename += str(round(random.random()*1000)) + ".list"
            temp_file = open(temp_filename, "w")
            first_file_name = ""

            # search for inout files
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

            file_out = open(project_path + "/script/" + temp + "file_" + str(i) + ".sh", "w")
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
            if submit:
                #file_out.write("source /cvmfs/sft.cern.ch/lcg/contrib/gcc/7.3.0binutils/x86_64-centos7-gcc7-opt/setup.sh\n")
                file_out.write("source /pbs/software/centos-7-x86_64/root/6.18.00/bin/thisroot.sh\n")
            file_out.write("cd " + bin_dir + "\n")
            file_out.write("#$ -e " + outpt_prefix+"/"+outpt_version + outpt_log + "\n")
            file_out.write("#$ -o " + outpt_prefix+"/"+outpt_version + outpt_log + "\n")
            file_out.write("source /path/to/your/install/setup.sh\n")
            file_out.write(command + "\n")

            file_out.close()
            #print(file_out)

            launcher.write("chmod 765 ./file_" + str(i) + ".sh\n")
            launcher.write("./file_" + str(i) + ".sh\n")
            #launcher.write("email " + "\"Job done\" " + "\"ot_file\"\n")
            if submit:
               job = project_path + "/script/" + temp + "file_" + str(i) + ".sh\n"
               subprocess.call("qsub  -P P_t2k -l sps=1 -l ct=00:30:00  %s"%job,  shell=True)
            i+=1
            

        launcher.close()

    # create Condor submission file
    #submit_file = open(project_path + "/script/" + temp + "/Submit.sub", "w")

    #submit_file.write("executable              = $(filename)\n")

    #submit_file.write("arguments               = $(ClusterId)$(ProcId)\n")
    #submit_file.write("requirements            = (OpSysAndVer =?= \"CentOS7\")\n")
    #submit_file.write("+job_flavour             = \"" + job_flavour + "\"\n")
    #submit_file.write("output                  = " + log_folder + "/output/$(ClusterId).$(ProcId).out\n")
    #submit_file.write("error                   = " + log_folder + "/error/$(ClusterId).$(ProcId).err\n")
    #submit_file.write("log                     = " + log_folder + "/log/$(ClusterId).log\n")
    #submit_file.write("queue filename matching files *.sh")

    #submit_file.close()

    os.chdir(project_path + "/script/" + temp)
    if launch:
        subprocess.run(["chmod", "765", "./all.sh"], check=True)
        subprocess.run(["/bin/bash", "all.sh"], check=True)
    #if submit:
    #    subprocess.run(["qsub",  "Submit.sub"], check=True)
    os.chdir(project_path + "/script/")
    return 0


if __name__ == "__main__":
    main()
