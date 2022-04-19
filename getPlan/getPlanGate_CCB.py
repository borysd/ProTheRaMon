#!/usr/bin/env python3
import sys
import os
import optparse
import pandas as pd
import pydicom as dicom
import numpy as np
from shutil import copyfile

# parse input
inputParse = optparse.OptionParser()
inputParse.add_option('--fieldToSave', '-f', default=None)
inputParse.add_option('--configPath', '-c', default=None)
inputParse.add_option('--Patientfile', '-p', default=None)
inputParse.add_option('--CTfile', '-t', default=None)
inputParse.add_option('--actorList', '-a', action ="append", default=None)
inputParse.add_option('--beammodel', '-b', default="default")
inputParse.add_option('--simLabel', '-l', default="Sim0000")
inputParse.add_option('--CPUs', '-u', default=1)
inputParse.add_option('--nparticles', '-n', default="0.1")
inputParse.add_option('--Nunit', '-N', default="primPr")  # "primPr"   OR "primNoPB"
inputParse.add_option('--dx', '-x', default="0.0")
inputParse.add_option('--dy', '-y', default="0.0")
inputParse.add_option('--dz', '-z', default="0.0")
inputParse.add_option('--HU', '-m', default="0.0")

options, fileList = inputParse.parse_args()

if options.fieldToSave == None:
    fieldToSave = options.fieldToSave
elif options.fieldToSave.isdigit():
    fieldToSave = int(options.fieldToSave)
else:
    fieldToSave = options.fieldToSave

# get Config Path
if options.configPath == None:
    # configPath = "current"
    configPath = os.getcwd()
else:
    configPath = options.configPath

# get Patient DICOM RTplan PATH
if options.Patientfile == None:
    print("Missing DICOM RT PLAN file path !")
    exit(0)
else:
    PATIENTfileList = os.path.abspath( options.Patientfile )

# get CT .mhd file name
# for example:  CT_ExtROICrop_2.0x2.0x2.0.mhd 
if options.CTfile == None:
    print("Missing MHD file path !")
    exit(0)
else:
    CTfileList = os.path.join( PATIENTfileList, "TPS", options.CTfile )

# take simLabel to create apropriate subfolder for your simulations results
SimLabelList = options.simLabel

# get number of CPUs (cores) - for appropriate calculating number of primaries per simulation
if options.CPUs.isdigit():
    CPUsNumber = int(options.CPUs)
else:
    CPUsNumber = int(1)

# get Actors List
# for example:  doseAll doseP LETAll LETP  
ActorsList = []
if options.actorList == None:
    print("Missing Actor List!")
    # ActorsList = ["doseAll"]
    exit(0)
else:
    ActorsList = options.actorList

print( ActorsList )


# get number of number of primaries in total per simulation [1E4  / 0.1 ]
if options.nparticles.isdigit():
    ParticleNum = float(options.nparticles)
else:
    ParticleNum = float(0.1)
print(ParticleNum)

# get units for Particle numbers [ primNoPB / primPr ]  = Number / Percent
if options.Nunit == None:
    ParticleNumUnit = "primPr"
else:
    ParticleNumUnit = options.Nunit


dx = float(options.dx)
dy = float(options.dy)
dz = float(options.dz)
dHU = float(options.HU)
robustXYZlist = [dx, dy, dz, dHU]

# #ASSUMPTIONS:
# #first arg is the RTplanFileName with full path
# # second arg is CT (in MHD/RAW format) with full path
# if len(fileList) == 1:
#     PATIENTfileList = os.path.abspath(fileList[0])
#     print("Missing MHD file path !")
# elif len(fileList) == 2:    
#     PATIENTfileList = os.path.abspath(fileList[0])
#     CTfileList = fileList[1]
#     # SimLabelList = "Sim0000"  #default Sim Label
# else:
#     PATIENTfileList = os.path.abspath(fileList[0])
#     # /home/user/patientData/G3P052E0
#     CTfileList = fileList[1]
#     # for example:  CT_ExtROICrop_2.0x2.0x2.0.mhd 
#     SimLabelList = fileList[2]
#     # for example  Sim0100


from getPlanGateLib_CCB import getPlan
getPlan(PatientPath=PATIENTfileList,
        CTpatientFileName=CTfileList,
        simLabel=SimLabelList,
        fieldToSave=fieldToSave,
        GATEconfigPath=configPath,
        cores=CPUsNumber,
        actors=ActorsList,
        Nparticles=ParticleNum,
        ParticleUnit=ParticleNumUnit,
        robustXYZ=robustXYZlist,
        GATErtplanFileName='current',
        defaultModelsFileName=options.beammodel,
        beamModelInterpolationMethod='slinear',
        displayInfo=True)
