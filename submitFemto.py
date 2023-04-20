#!python

import os
from shutil import copytree, copy

#------------------------set input/output dirs here----------------------#

scriptdir = "/lustre/hades/user/kjedrzej/defallRoot6.24.02Hydra2-6.1.sh"
submissiondir="/lustre/hades/user/kjedrzej/submit/femtoiHKM"
inputdir="/lustre/hades/user/kjedrzej/iHKM/ic/output/14p5GeV_0to10cent"
outputdir="/lustre/hades/user/kjedrzej/iHKM/14p5GeV_0to10cent_femto"
jobarrayfile="jobarrayFemto.dat"

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

try:
    os.mkdir(submissiondir) # create directory for Slurm job submission
except Exception as e:
    print(e)
try:
    os.mkdir(outputdir) # create output directory
except Exception as e:
    print(e)
try:
    os.mkdir(outputdir+"/out") # create direcotry for Slurm output files
except Exception as e:
    print(e)
try:
    os.remove(jobarrayfile) # remove old Slurm job array file
except Exception as e:
    print(e)

file = open(jobarrayfile,"w")

#----------------------change parameters here---------------------------#

kTlist = [19] # list of kT cuts, used in /therm2_femto, e.g kTlist = [0,1,2,3]
dirlist = [dirnum for dirnum in range(1,6)] # list of output dirs created by divideIHKM.sh script; minimal value is included, maximal value is excluded
partType = ["proton-proton"] # list of particle types for analysis

NoFiles = 100 # amout of files in each direcotry created by divideIHKM.sh script; now deflaut is 100

#-----------------------------------------------------------------------#

if (not os.path.exists(inputdir)):
    raise Exception("Input direcotry does not exist")   # check if input direcory exists

for dir in dirlist:
    inpDir = inputdir + "/out_" + str(dir) # make sure you've used tge divideIHKM.sh script before runing this macro
    print(inpDir)
    
    for kT in kTlist:
        for pair in partType:
            if pair is "pion-pion":
                pions = "pipi"
            elif pair is "pionM-pionM":
                pions = "pimpim"
            elif pair is "proton-proton":
                pions = "pp"
            else:
                raise Exception("Unkonwn particle type")

            femtoDir = "{}/{}_{}_kT{}".format(outputdir,pions,dir,kT)
            print(femtoDir)
            try:
                os.mkdir(femtoDir)
            except Exception as e:
                print(e)
            try:
                copy("therm2_femto",femtoDir)
                replaceFemto(inpDir,pair,kT,NoFiles)
                copy("femto.ini",femtoDir)
                file.write(femtoDir+"\n")
            except Exception as e:
                print(e)
                
file.close()
copy(jobarrayfile, submissiondir)
copy("jobScriptFemto_SL.sh", submissiondir)
#copy("wrap.sh",submissiondir)
chmod_command = "chmod 777 " + submissiondir + "/jobScriptFemto_SL.sh" # this command might be old and unnecessary, but afaik it is recommended
print(chmod_command)
os.system(chmod_command)

arrElem = len(kTlist) * len(dirlist) * len(partType) # calculate the amount of jobs for submission

command="sbatch --array=1-{3} --mem=4000 --time=0-8:00:00 -D {0}  --output={1}/out/slurm-%A_%a.out -- {0}/jobScriptFemto_SL.sh {0}/jobarrayFemto.dat {1} {2}".format(submissiondir,outputdir,scriptdir,arrElem+1)
print(command)
os.system(command)


