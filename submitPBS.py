#!python

import os
from shutil import copytree, copy

#------------------------set input/output dirs here----------------------#

scriptdir = "/lustre/hades/user/kjedrzej/defallRoot6.24.02Hydra2-6.1.sh"
submissiondir="/lustre/hades/user/kjedrzej/submit/femtoiHKM"
inputdir="/lustre/hades/user/kjedrzej/iHKM/ic/output/14p5GeV_0to10cent"
outputdir="/lustre/hades/user/kjedrzej/iHKM/14p5GeV_0to10cent_femto"
jobfile="sendJob.pbs"

#------------------------------------------------------------------------#

def replaceFemto(path,ptype,nBin,events): # method for creating separate femto.ini files for each job
    fin = open("femto_template", "rt")
    fout = open("femto.ini", "wt")

    for line in fin:
        if "DIR" in line:
            fout.write(line.replace('DIR', str(path + "/")))
        elif "TYPE" in line:
            fout.write(line.replace('TYPE',str(ptype)))
        elif "TRANS" in line:
            fout.write(line.replace('TRANS',str(nBin)))
        elif "EVENTS" in line:
            fout.write(line.replace('EVENTS',str(events)))
        else:
            fout.write(line)
    fin.close()
    fout.close()

#----------------------change parameters here---------------------------#

kTlist = [19] # list of kT cuts, used in /therm2_femto, e.g kTlist = [0,1,2,3]
partType = ["proton-proton"] # list of particle types for analysis
NoFiles = 100 # amout of files in each direcotry created by divideIHKM.sh script; now deflaut is 100

mem = 16 # amout of RAM allocated for each job (in GB)
wall = 48 # time after which the job will be terminated (in hours)
name = "14p5GeVfemto" #name of the job

#-----------------------------------------------------------------------#