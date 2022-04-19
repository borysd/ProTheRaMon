import fredtools as ft
import glob
import subprocess
import time
from zipfile import ZipFile
import os
import optparse
from multiprocessing import Process


# parse input
inputParse = optparse.OptionParser()
inputParse.add_option('--PatientPath', '-p', default=None)
inputParse.add_option('--patient', '-P', default="G3P007E1")
options, fileList = inputParse.parse_args()


#patPath = options.PatientPath
patPath = "/home/jpet/patientData/"
#SimLabel = options.simLabel

patientID = options.patient

CTfile="CT_ExtROICrop_2.5x2.5x2.5.mhd"

#patients= glob.glob("/home/jpet/patientData/G3P013E1")

#subfoldersPat= [f.path for f in os.scandir( patPath ) if f.is_dir()]

def readMHDfile( CTfileMHD ):
    import SimpleITK as sitk

    MHDimage = sitk.ReadImage(CTfileMHD)        
    dims = MHDimage.GetSize()   # get Dimensions from MHD file
    spacing = MHDimage.GetSpacing()  # get Element Spacing from MHD file

    return [dims, spacing]



if os.path.exists( os.path.join(patPath, patientID, "TPS", CTfile) ):
    fname = os.path.join(patPath, patientID, "TPS", CTfile)
    [dims, spacing] = readMHDfile(fname)
    print( dims )

# -----------------------------------------------------------------------




#### --------------------------------------------------------------------------------------------------

#processes = []

#for onePatient in patients:
#    print( onePatient )
#    p = Process(target=onePatientProcess, args=(onePatient,))
#    p.start()
#    processes.append(p)

#for p in processes:
#    p.join()



# -----------------------------------------------------------------------
