import sys
import os
from os.path import basename
# import math
import optparse
from numpy.core.numeric import False_
import pandas as pd
import pydicom as dicom
import numpy as np
import SimpleITK as sitk
from shutil import copyfile
import scipy.io as sio
# import medpy
from multiprocessing import Process
import shutil
import subprocess
import glob
import time
from zipfile import ZipFile
import fredtools as ft

isotopes = ["O15", "O14", "N13", "C11", "C10", "P30", "K38"]


######  DEF ######  -------------------------  read DimSize from MHD text file ------------------------------------
def readMHDfile( CTfileMHD ):
    import SimpleITK as sitk

    MHDimage = sitk.ReadImage(CTfileMHD)        
    dims = MHDimage.GetSize()   # get Dimensions from MHD file
    spacing = MHDimage.GetSpacing()  # get Element Spacing from MHD file

    return [dims, spacing, MHDimage]

### --------------------------------------------------------------------------------------------------------------

# parse input
inputParse = optparse.OptionParser()
inputParse.add_option('--PatientPath', '-p', default=None)
inputParse.add_option('--CTfile', '-t', default=None)
inputParse.add_option('--simLabel', '-l', default="Sim0100")
# inputParse.add_option('--configPath', '-c', default=None)
options, fileList = inputParse.parse_args()

# get Patient DICOM RTplan PATH
if options.PatientPath == None:
    print("Missing DICOM RT PLAN file path !")
    # exit(0)
else:
    PATIENTfile = os.path.abspath( options.PatientPath )

print( PATIENTfile )
# PATIENTfile="/Users/dborys/Documents/POLSL/2020_IFJ_postdoc/getPlanTests/PETscripts/compare/myscript/production_maps/G3P052E0"

# get CT .mhd file name
# for example:  CT_ExtROICrop_2.0x2.0x2.0.mhd 
if options.CTfile == None:
    print("Missing MHD file path !")
    # exit(0)
else:
    CTfileList = os.path.join( PATIENTfile, "TPS", options.CTfile )
# CTfileList = "/Users/dborys/Documents/POLSL/2020_IFJ_postdoc/getPlanTests/patientData/G3P052E0/TPS/CT_ExtROICrop_2.5x2.5x2.5.mhd"


CTfilePath = os.path.join( PATIENTfile, "TPS", CTfileList ) 
# read x,y,z dims from MHD file
[MHDdims, MHDspacing, MHDimage] = readMHDfile( CTfilePath )

# take simLabel to create apropriate subfolder for your simulations results
SimLabel = options.simLabel
# SimLabel='Sim1000'

# FOR EVERY FIELD IN PATIENT/gate/SimXXX/FIELDname/Results
SimPath = os.path.join( PATIENTfile, "gate", SimLabel)

SimPathFRED = os.path.join( PATIENTfile, "fred", SimLabel)

isotopeList = ["C11","C10","N13","O15","O14","P30","K38"]
isotopeT1_2Values = [ 1220.04, 19.29, 597.9, 122.24, 70.598, 149.88, 458.16 ]   # T_iso


# ----------------------------- aggregateProd_Isotope ----------------------------------------------------------------


def aggregateProd_Isotope(SimPath, isotope, MHDdims, files, resultsPath):
    isoName = isotope + "_map"
    matrixSum = np.zeros(( MHDdims[1], MHDdims[2], MHDdims[0]),dtype = np.int16)
    # matrixSum = np.zeros(( MHDdims[1]*MHDdims[2], MHDdims[0]),dtype = np.int16)
    # matrixSum = np.zeros((140, 114, 121),dtype = int)
    isAnyProd = False

    for filename in files:
        # if "ISOTOPE_map" in filename and "-Prod.txt" in filename:
        if isoName in filename and "-Prod.txt" in filename:
            # print( filename )
            fileTxt = os.path.join(resultsPath, filename)
            matrix = np.loadtxt( fileTxt, dtype=np.int16 )
            matrix = np.reshape( matrix, (MHDdims[1], MHDdims[2], MHDdims[0]))
            matrixSum += matrix
            isAnyProd = True
    
    if isAnyProd:
    # matrixSumPermute = matrixSum.transpose(2, 0, 1)
        # matrixSumPermute = matrixSum.T
        matrixSumPermute = matrixSum
        fileOut = os.path.join(resultsPath, isotope + "_total.bin")
        matrixSumPermute.astype('int16').tofile(fileOut)


# ----------------------------- aggregateProd --------------------------------------------------------------------

def aggregateProd(SimPath, isotopes, MHDdims):
    subfolders= [f.path for f in os.scandir( SimPath ) if f.is_dir()]
    # print( subfolders )
    processes = []

    for dirname in list(subfolders):
        print( dirname )
        resultsPath = os.path.join( dirname, "Results")  #dirname is FIELDname
        #MERGE files for every isotope

        files = os.listdir( resultsPath )
        for isotope in isotopes:
            print(isotope)
            # aggregateProd_Isotope(SimPath, isotope, MHDdims, files, resultsPath) 
            p = Process(target=aggregateProd_Isotope, args=(SimPath, isotope, MHDdims, files, resultsPath))
            p.start()
            processes.append(p)
        
    for p in processes:
        p.join()


    time.sleep(30)
    # p.join()


# ----------------------------- compressProd_Isotope ----------------------------------------------------------------

def compress_Prod_isotope(SimPath, isotope):
    #compress isotopes (without dose)
    fileOutZIP = os.path.join(SimPath, isotope + "_map.zip")
    files_maps = os.path.join(SimPath, isotope + "_map*.txt")
    # bashCmd = ["zip -Tm "+ fileOutZIP + "  "+ files_maps]
    # process = subprocess.Popen(bashCmd, shell=True, stdout=subprocess.PIPE)

    subprocess.call(['zip', '-Tm', fileOutZIP] + glob.glob(files_maps))

    # (output, err) = process.communicate()
    # print(output)
    # print(err)


# ----------------------------- decompressProd_Isotope ----------------------------------------------------------------



def decompress_Prod_isotope(SimPath, isotope):
    #compress isotopes (without dose)
    fileOutZIP = os.path.join(SimPath,   isotope + "_map.zip")

    #unzip files
    if os.path.isfile(fileOutZIP):
        # bashstr = 'unzip \'' + fileOutZIP + '\' -d /'        
        # subprocess.call(['unzip',  fileOutZIP, '-d /'])

        with ZipFile(fileOutZIP, 'r') as f:
            #extract in different directory
            f.extractall(fileOutZIP[0])

        # remove .zip file
        subprocess.call(['rm', fileOutZIP])

# ----------------------------- compress_stats_Dose ----------------------------------------------------------------

def compress_stats_Dose(SimPath):
    #compress stats
    fileOutZIP = os.path.join(SimPath, "simStat.zip")
    files_maps = os.path.join(SimPath, "simStat*.out")
    bashCmd = ["zip -Tm "+ fileOutZIP + "  "+ files_maps]
    # process = subprocess.Popen(bashCmd, shell=True, stdout=subprocess.PIPE)    
    
    subprocess.call(['zip', '-Tm', fileOutZIP] + glob.glob(files_maps))

    #compress dose
    fileOutZIP = os.path.join(SimPath, "dose.zip")
    files_maps = os.path.join(SimPath, "*-Dose.*")
    bashCmd = ["zip -Tm "+ fileOutZIP + "  "+ files_maps]
    # process = subprocess.Popen(bashCmd, shell=True, stdout=subprocess.PIPE)     

    subprocess.call(['zip', '-Tm', fileOutZIP] + glob.glob(files_maps))
    
            
    
# ----------------------------- compressProd --------------------------------------------------------------------

def compressProd(SimPath, isotopes):
    subfolders= [f.path for f in os.scandir( SimPath ) if f.is_dir()]
    # print( subfolders )
    for dirname in list(subfolders):
        print( dirname )
        resultsPath = os.path.join( dirname, "Results")  #dirname is FIELDname
        #COMPRESS files for every isotope
        
        # files = os.listdir( resultsPath )
        for isotope in isotopes:
            print(isotope)
            # compress_Prod_isotope(resultsPath, isotope)
            p = Process(target=compress_Prod_isotope, args=(resultsPath, isotope))
            p.start()
            # p.join()

        # compress_stats_Dose(resultsPath)
        print("Compress stats and dose")        
        p = Process(target=compress_stats_Dose, args=(resultsPath,))
        p.start()
        

# ----------------------------- decompressProd --------------------------------------------------------------------

def decompressProd(SimPath, isotopes):
    subfolders= [f.path for f in os.scandir( SimPath ) if f.is_dir()]
    # print( subfolders )
    processes = []

    for dirname in list(subfolders):
        print( dirname )
        resultsPath = os.path.join( dirname, "Results")  #dirname is FIELDname
        #COMPRESS files for every isotope
        
        # files = os.listdir( resultsPath )
        for isotope in isotopes:
            print(isotope)
            # compress_Prod_isotope(resultsPath, isotope)
            p = Process(target=decompress_Prod_isotope, args=(resultsPath, isotope,))
            p.start()
            processes.append(p)
        

    for p in processes:
        p.join()
    time.sleep(20)
    # p.join()


# ----------------------------- compressDose --------------------------------------------------------------------

def compressDose(SimPath):
    subfolders= [f.path for f in os.scandir( SimPath ) if f.is_dir()]
    # print( subfolders )
    for dirname in list(subfolders):
        print( dirname )
        resultsPath = os.path.join( dirname, "Results")  #dirname is FIELDname
        #COMPRESS files for all doses
        files = os.listdir( resultsPath )
        
        p = Process(target=compressDose_mhd, args=(SimPath, files, resultsPath))
        p.start()

# ----------------------------- writeInterfile --------------------------------------------------------------------

def writeInterfile(fname, file_path, sizeXYZ, voxelSize, matrixActivity):

    # WRITE HEADER
    interfile_fileName = os.path.abspath(os.path.join(file_path, fname + ".h33"))
    interfile_header=open(interfile_fileName,'w')
    
    print('!INTERFILE :=', file=interfile_header)
    print('!GENERAL DATA :=', file=interfile_header)
    print('!data offset in bytes := 0', file=interfile_header)
    # fname_img = fname + ".i33"
    fname_img = os.path.abspath(os.path.join(file_path, fname + ".i33"))
    print('!name of data file := {:s}'.format(fname_img), file=interfile_header)
    print(';', file=interfile_header)
    print('!GENERAL IMAGE DATA :=', file=interfile_header)
    print('!total number of images := {:d}'.format(sizeXYZ[2]), file=interfile_header)
    print('imagedata byte order := LITTLEENDIAN', file=interfile_header)
    print(';', file=interfile_header)
    print('number of energy windows := 1', file=interfile_header)
    print(';', file=interfile_header)
    print('!number of images/energy window := {:d}'.format(sizeXYZ[2]), file=interfile_header)
    print('!matrix size [1] := {:d}'.format(sizeXYZ[0]), file=interfile_header)
    print('!matrix size [2] := {:d}'.format(sizeXYZ[1]), file=interfile_header)
    # print('!number format := UNSIGNED INTEGER', file=interfile_header)
    print('!number format := FLOAT', file=interfile_header)
    # print('!number of bytes per pixel := 2', file=interfile_header)
    print('!number of bytes per pixel := 4', file=interfile_header)
    print('scaling factor (mm/pixel) [1] := {:f}'.format(voxelSize[0]), file=interfile_header)
    print('scaling factor (mm/pixel) [2] := {:f}'.format(voxelSize[1]), file=interfile_header)
    print('!number of projections := {:d}'.format(sizeXYZ[2]), file=interfile_header)
    print('!number of slices := {:d}'.format(sizeXYZ[2]), file=interfile_header)
    print('slice thickness (pixels) := {:f}'.format(voxelSize[2]), file=interfile_header)
    print('!END OF INTERFILE :=', file=interfile_header)
    interfile_header.close()

    # WRITE RAW DATA
    interfile_fileNameRAW = os.path.abspath(os.path.join(file_path, fname + ".i33"))
    # interfile_raw=open(interfile_fileNameRAW,'wb')
    matrixActivity.astype('float32').tofile(interfile_fileNameRAW)
    # interfile_raw.close()



# -------------------------------   makeSource   --------------------------------------------------------------------

def makeSource( dimX, dimY, dimZ, voxelSize, simPath, isotopesList, isoT_12_Values ):
    # get number of Fields
    subfolders= [f.path for f in os.scandir( simPath ) if f.is_dir()]
    subfolders.sort()
    Nfields = len(subfolders)
    t_scan = 60


    # allBinFiles = 0

    # # iterate until all _total.bin files will be ready
    # while allBinFiles == 0: 
    #     allBinFiles = 1
    #     subfolders1= [f.path for f in os.scandir( simPath ) if f.is_dir()]
    #     # print( subfolders )
    #     for dirname in list(subfolders1):
    #         print( dirname )
    #         resultsPath = os.path.join( dirname, "Results")  #dirname is FIELDname
    #         for isotope in isotopes:            
    #             print(isotope)
    #             #compress isotopes (without dose)
    #             fileIsoBIN = os.path.join(resultsPath, isotope + "_total.bin")

    #             if os.path.isfile(fileIsoBIN):
    #                 allBinFiles = allBinFiles * 1
    #             else:
    #                 allBinFiles = allBinFiles * 0                


    times = []
    times.append(t_scan)   # times[0] == t_scan

    for n_iso in range( 0, Nfields-1 ):
        times.append( times[n_iso] + 90)
    times.sort(reverse=True)


    for idx_iso in range(0,len(isotopesList)):
        isAnyDir = False
        isoFile =  isotopesList[idx_iso] + "_scorer.mhd" 
        L_iso = np.log(2)/float(isoT_12_Values[idx_iso])
        
        
        matrixActivity = np.zeros((dimX*dimY*dimZ),dtype = np.float32)
        matrixEmission = np.zeros((dimX*dimY*dimZ),dtype = np.float32)
        
        for idx_dir in range( 0, len(subfolders) ):
        # for dirname in list( subfolders ):    #1st subfolder == 1Field ... etc
            dirname = subfolders[ idx_dir ]
            isAnyDir = True
            resultsPath = os.path.join( dirname, "out", "reg", "Phantom" )

            # matrix = np.zeros((dimX, dimY, dimZ),dtype = np.float32)
            # print( matrix.shape )

            fileIn = os.path.join( resultsPath, isoFile )
            # matrix = np.fromfile(fileIn, dtype='uint32')
            # matrix_tmp = np.fromfile(fileIn, dtype='uint16')
            
            matrix_tmp = ft.readMHD( fileIn )
            matrix = ft.arr( matrix_tmp ).ravel()
            print( matrix.shape )

            # calculate Activity map
            matrixActivityTemp = ( ( L_iso * np.exp( -L_iso*times[idx_dir] ) * matrix ) )#.astype( int )
            print( matrixActivityTemp.shape )
            matrixActivity += matrixActivityTemp

            # calculate Emission map
            matrixEmissionTemp =  ( np.exp( -L_iso*times[idx_dir] ) ) * matrix  * ( 1 - np.exp( -L_iso*t_scan ) )            
            matrixEmission += matrixEmissionTemp

        if isAnyDir:
            sizeXYZ = (dimX, dimY, dimZ)

            fnameActivity = isotopesList[idx_iso] + "_total_act"
            fnameEmission = isotopesList[idx_iso] + "_total_emi" 
            fileOutAct = os.path.join( simPath, fnameActivity+".bin" )
            fileOutEmis = os.path.join( simPath, fnameEmission+".bin" )
            # matrixActivity.astype('float32').tofile( fileOutAct )
            matrixEmission.astype('float32').tofile( fileOutEmis )
            writeInterfile(fnameActivity, simPath, sizeXYZ, voxelSize, matrixActivity)


## ____________________ new makeSourceMHD _____________________________________________________________________________



def makeSourceMHD( CTfilePath, simPath, isotopesList, isoT_12_Values ):
    # get number of Fields
    subfolders= [f.path for f in os.scandir( simPath ) if f.is_dir()]
    subfolders.sort()
    Nfields = len(subfolders)
    t_scan = 60

    times = []
    times.append(t_scan)   # times[0] == t_scan

    for n_iso in range( 0, Nfields-1 ):
        times.append( times[n_iso] + 90)
    times.sort(reverse=True)

    [MHDdims, MHDspacing, MHDimage] = readMHDfile( CTfilePath )
    dimX = MHDdims[0]
    dimY = MHDdims[1]
    dimZ = MHDdims[2]

    voxelSize = MHDspacing[0]

    for idx_iso in range(0,len(isotopesList)):
        isAnyDir = False
        isoFile = isotopesList[idx_iso] + "_total.bin" 
        L_iso = np.log(2)/float(isoT_12_Values[idx_iso])
        
        

        ct_scan = sitk.GetArrayFromImage(MHDimage)
        matrixActivity = np.zeros((ct_scan.shape))
        matrixEmission = np.zeros((ct_scan.shape))

        # matrixActivity = np.zeros((dimX*dimY*dimZ),dtype = np.float32)
        # matrixEmission = np.zeros((dimX*dimY*dimZ),dtype = np.float32)
        
        for idx_dir in range( 0, len(subfolders) ):
        # for dirname in list( subfolders ):    #1st subfolder == 1Field ... etc
            dirname = subfolders[ idx_dir ]
            isAnyDir = True
            resultsPath = os.path.join( dirname, "Results" )

            # matrix = np.zeros((dimX, dimY, dimZ),dtype = np.uint32)
            # matrix = np.zeros((ct_scan.shape),dtype='uint16')
            # print( matrix.shape )

            fileIn = os.path.join( resultsPath, isoFile )
            # matrix = np.fromfile(fileIn, dtype='uint32')
            # matrix = np.fromfile(fileIn, dtype='uint16').reshape( ct_scan.shape )
            matrix = np.fromfile(fileIn, dtype='uint16').reshape( ct_scan.shape[2],ct_scan.shape[1],ct_scan.shape[0] ).transpose( 2,1,0)
            # matrix = np.reshape( matrix, ct_scan.shape )
            print( matrix.shape )

            # calculate Activity map
            matrixActivityTemp = ( 10.0 * ( L_iso * np.exp( -L_iso*times[idx_dir] ) * matrix ) )#.astype( int )
            matrixActivity += matrixActivityTemp

            # calculate Emission map
            matrixEmissionTemp =  10.0 * ( np.exp( -L_iso*times[idx_dir] )) * matrix  * ( 1 - np.exp( -L_iso*t_scan ) )            
            matrixEmission += matrixEmissionTemp

        if isAnyDir:
            sizeXYZ = (dimX, dimY, dimZ)

            fnameActivity = isotopesList[idx_iso] + "_total_act"
            fnameEmission = isotopesList[idx_iso] + "_total_emi" 
            fileOutAct = os.path.join( simPath, fnameActivity+".mhd" )
            fileOutEmis = os.path.join( simPath, fnameEmission+".mhd" )

            actmap = sitk.GetImageFromArray(matrixActivity) #.astype('int16'))
            actmap.SetOrigin( MHDimage.GetOrigin() )
            actmap.SetSpacing( MHDimage.GetSpacing() )            

            emimap = sitk.GetImageFromArray(matrixEmission) #.astype('float32'))
            emimap.SetOrigin( MHDimage.GetOrigin() )
            emimap.SetSpacing( MHDimage.GetSpacing() )            

            writer = sitk.ImageFileWriter()
            writer.SetFileName(fileOutAct)
            writer.Execute(actmap)

            writer = sitk.ImageFileWriter()
            writer.SetFileName(fileOutEmis)
            writer.Execute(emimap)            


# ################### ----------------------------------------------------------------------------- ##################
# ################### ----------------------------------------------------------------------------- ##################
# ################### ----------------------------------------------------------------------------- ##################

# decompressProd(SimPath, isotopes)
# aggregateProd(SimPath, isotopes, MHDdims)
# compressProd(SimPath, isotopes)

dimX = MHDdims[0]
dimY = MHDdims[1]
dimZ = MHDdims[2]
# voxelSize = [2.0, 2.0, 2.0] # MHDspacing
voxelSize = [MHDspacing[0], MHDspacing[1], MHDspacing[2]]
simPath = SimPathFRED
# simPath="/Users/dborys/Documents/POLSL/2020_IFJ_postdoc/getPlanTests/PETscripts/compare/myscript/production_maps/G3P052E0/inroom"

isotopesList = isotopeList
isoT_12_Values = isotopeT1_2Values

makeSource( dimX, dimY, dimZ, voxelSize, simPath, isotopesList, isoT_12_Values )
# makeSourceMHD( CTfilePath, simPath, isotopesList, isoT_12_Values )

