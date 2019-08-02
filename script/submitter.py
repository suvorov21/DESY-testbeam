#! /usr/bin/python

import argparse
import subprocess
import shutil
import os
import random

if __name__ == "__main__":

  #define input
  Niter     = 5
  doiter    = True

  bin_dir   = "/afs/cern.ch/work/s/ssuvorov/public/T2K_testbeam/DESY_TestBeam/bin/"
  bin_name  = "SpatialResol_circ.exe"
  bin_flag  = "-b"

  input_prefix  = "/eos/experiment/neutplatform/t2knd280/DESY_TPC/ROOT/"
  input_version = "v1"

  outpt_prefix  = "/eos/user/s/ssuvorov/DESY_testbeam/"
  outpt_version = "v3"

  # espresso     = 20 minutes
  # microcentury = 1 hour
  # longlunch    = 2 hours
  # workday      = 8 hours
  # tomorrow     = 1 day
  # testmatch    = 3 days
  # nextweek     = 1 week
  JobFlavour    = "longlunch"
  log_folder    = "/afs/cern.ch/work/s/ssuvorov/public/T2K_testbeam/DESY_TestBeam/script/"
  # end of input definition

  parser = argparse.ArgumentParser(description='Submit jobs to condor at LXPLUS')
  parser.add_argument("-f", metavar="f", type=str,
    help='File list with input/output files names',
    required=True)

  args = parser.parse_args()

  project_path = os.path.abspath(os.path.dirname(os.path.abspath(__file__)) + "/../")
  if os.path.exists(project_path + "/script/temp"):
    shutil.rmtree(project_path + "/script/temp")
  os.mkdir(project_path + "/script/temp");

  # for each input file
  with open(args.f) as fl:
    i = 0
    for line in fl:
      in_file  = line.split()[0]
      ot_file  = line.split()[1]

      if (in_file == "" or ot_file == ""):
        Print("Error")

      # create a file list in case of existing subruns
      temp_filename = project_path + "/FileLists/temp" + str(round(random.random()*1000)) + ".list"
      temp_file = open(temp_filename, "w")
      path = input_prefix+"/"+input_version+"/"

      for r, d, f in os.walk(path):
        if "soft" in r:
          continue
        for file in f:
          if in_file[:21] in file and ".root" in file:
            temp_file.write(os.path.join(r, file) + "\n")
      temp_file.close()

      file_out = open(project_path + "/script/temp/" + str(i) + ".sh", "w")
      i+=1
      command = ""

      # fo each iteration
      for it in range(0, Niter):
        command += bin_dir + "/" + bin_name + " " + bin_flag + " -t " + str(it)
        command += " -i " + temp_filename
        command += " -o " + outpt_prefix+"/"+outpt_version+"/"+ot_file
        if (doiter):
           command += "_iter" + str(it)
        command += ".root; "
      # rm temp file list
      command += "rm " + temp_filename

      file_out.write("#!/bin/bash\n")
      file_out.write("source /cvmfs/sft.cern.ch/lcg/contrib/gcc/7.3.0binutils/x86_64-centos7-gcc7-opt/setup.sh\n")
      file_out.write("source /afs/cern.ch/work/s/ssuvorov/public/ROOT/root-6.18.00-build_gcc73/bin/thisroot.sh\n")
      file_out.write("cd /afs/cern.ch/work/s/ssuvorov/public/T2K_testbeam//DESY_TestBeam/bin/\n")
      file_out.write(command + "\n")

      file_out.close()

  submit_file = open(project_path + "/script/temp/Submit.sub", "w")

  submit_file.write("executable              = $(filename)\n")

  submit_file.write("arguments               = $(ClusterId)$(ProcId)\n")
  submit_file.write("requirements            = (OpSysAndVer =?= \"CentOS7\")\n")
  submit_file.write("+JobFlavour             = \"" + JobFlavour + "\"\n")
  submit_file.write("output                  = " + log_folder + "/output/$(ClusterId).$(ProcId).out\n")
  submit_file.write("error                   = " + log_folder + "/error/$(ClusterId).$(ProcId).err\n")
  submit_file.write("log                     = " + log_folder + "/log/$(ClusterId).log\n")
  submit_file.write("queue filename matching files *.sh")

  submit_file.close()

  os.chdir(project_path + "/script/temp")
  subprocess.run(["condor_submit",  "Submit.sub"])
  os.chdir(project_path + "/script/")
  shutil.rmtree(project_path + "/script/temp")
