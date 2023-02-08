#!python

import os
from shutil import copytree, copy

#----------------set input/output dirs here--------------#

scriptdir = "/lustre/hades/user/kjedrzej/defallRoot6.24.02Hydra2-6.1.sh"
submissiondir="/lustre/hades/user/kjedrzej/submit/femtoiHKM"
outputdir="/lustre/hades/user/kjedrzej/iHKM/ic/output/urqmd"
jobarrayfile="jobarrayFemto.dat"

#--------------------------------------------------------#

def replaceFemto(path,ptype,nBin,events):
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
    os.mkdir(submissiondir)
    #os.mkdir(outputdir)
    os.mkdir(outputdir+"/out")
    os.remove(jobarrayfile)
except Exception as e:
    print(e)

file = open(jobarrayfile,"w")

#--------------change parameters here-------------------#

kTlist = [19]
dirlist = [dirnum for dirnum in range(1,6)]
partType = ["proton-proton"]

NoFiles = 100

#-------------------------------------------------------#

for dir in dirlist:
    inpDir = outputdir + "/out_" + str(dir)
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
                copy("therm2_femto",femtoDir)
                replaceFemto(inpDir,pair,kT,NoFiles)
                copy("femto.ini",femtoDir)
                file.write(femtoDir+"\n")
            except Exception as e:
                print(e)
                
file.close()
copy(jobarrayfile, submissiondir)
copy("jobScriptFemto_SL.sh", submissiondir)
copy("wrap.sh",submissiondir)
chmod_command = "chmod 777 " + submissiondir + "/jobScriptFemto_SL.sh"
print(chmod_command)
os.system(chmod_command)

arrElem = len(kTlist) * len(dirlist) * len(partType)

command="sbatch --array=1-{3} --mem=4000 --time=0-8:00:00 -D {0}  --output={1}/out/slurm-%A_%a.out -- {0}/wrap.sh {0}/jobScriptFemto_SL.sh {0}/jobarrayFemto.dat {1} {2}".format(submissiondir,outputdir,scriptdir,arrElem+1)
print(command)
os.system(command)


