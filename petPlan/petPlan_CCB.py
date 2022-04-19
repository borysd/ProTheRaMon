import sys
import os
import math
import optparse
import pandas as pd
import pydicom as dicom
import numpy as np
import SimpleITK as sitk
from shutil import copyfile
import scipy.io as sio
import medpy



isotopes = ["O15", "O14", "N13", "C11", "C10", "P30", "K38"]


######  DEF ######  -------------------------  read DimSize from MHD text file
def readMHDfile( CTfileMHD ):
    import SimpleITK as sitk

    MHDimage = sitk.ReadImage(CTfileMHD)        
    dims = MHDimage.GetSize()   # get Dimensions from MHD file
    spacing = MHDimage.GetSpacing()  # get Element Spacing from MHD file

    offset = MHDimage.GetOrigin()  # get offset ?????    
    return [dims, spacing, offset]    

### --------------------------------------------------------------------------------------------------------------


# parse input
inputParse = optparse.OptionParser()
inputParse.add_option('--patientPath', '-p', default=None)
inputParse.add_option('--CTfile', '-t', default=None)
inputParse.add_option('--simLabel', '-l', default="Sim0000")
inputParse.add_option('--configPath', '-c', default=None)
inputParse.add_option('--simPath', '-s', default=None)
inputParse.add_option('--geometryName', '-g', default=None)
inputParse.add_option('--dx', '-x', default="0.0")
inputParse.add_option('--dy', '-y', default="0.0")
inputParse.add_option('--dz', '-z', default="0.0")
options, fileList = inputParse.parse_args()



# get Patient DICOM RTplan PATH
if options.patientPath == None:
    print("Missing DICOM RT PLAN file path !")
    # exit(0)
else:
    patientPath = os.path.abspath( options.patientPath )
# patientPath="/Users/dborys/Documents/POLSL/2020_IFJ_postdoc/getPlanTests/patientData/G3P072E0"
print(patientPath)

# get CT .mhd file name
# for example:  CT_ExtROICrop_2.0x2.0x2.0.mhd 
if options.CTfile == None:
    print("Missing MHD file path !")
    # exit(0)
else:
    CTpatientFileName = os.path.join( patientPath, "TPS", options.CTfile )
# CTpatientFileName = "CT_ExtROICrop_2.0x2.0x2.0.mhd"


# get Config Path
if options.configPath == None:
    # configPath = "current"
    print("Missing CONFIG path !")
    # exit(0)
    configPath = os.getcwd()
else:
    configPath = options.configPath
# configPath = "/Users/dborys/Documents/POLSL/2020_IFJ_postdoc/getPlanTests/GatePETSimConf"

# get Production Simulation Path
if options.simPath == None:
    print("Missing Production SIMULATION path !")
    # exit(0)
    # simPath = os.getcwd()
else:
    simPath = options.simPath
# simPath="/Users/dborys/Documents/POLSL/2020_IFJ_postdoc/getPlanTests/patientData/G3P072E0/gate/Sim1000"


# get Geometry name
# for example:  Modular_dual_1x6_simple 
if options.geometryName == None:
    print("Missing Geometry Name !")
    # exit(0)
else:
    geometryName = options.geometryName 
# geometryName  = "Modular_dual_1x6_simple"


# take simLabel to create apropriate subfolder for your simulations results
SimLabel = options.simLabel


dx = float(options.dx)
dy = float(options.dy)
dz = float(options.dz)
robustXYZlist = [dx, dy, dz]


localDir = os.path.join( patientPath, "gate", SimLabel)
if not os.path.exists(os.path.abspath( localDir )):
    os.makedirs(os.path.abspath( localDir ))


# isotopeList = ["C11","C10","N13","O15","O14","P30","K38"]
isotopeList = ["C11"]
isotopeT1_2Values = [ 122.04, 19.29, 597.9, 122.24, 70.598, 149.88, 458.16 ]   # T_iso




# ---------------------------------------------------------- DEFs ----------------------------------------------------------

def loadDicomRT(dicomRTFileName, displayInfo=False):
    import pydicom as dicom
    if displayInfo:
        print('# Loading RT plan from {:s}'.format(dicomRTFileName), end=' ')
    try:
        dicomRT=dicom.read_file(dicomRTFileName)
        if dicomRT.SOPClassUID=='1.2.840.10008.5.1.4.1.1.481.8':
            if displayInfo:
                print('...OK')
            return dicomRT
        else:
            print('\nError while loading RT plan from {:s}. The file does not contain RT plan.'.format(dicomRTFileName))
            exit(1)
    except:
        print('\nError while loading RT plan from {:s}. No such file.'.format(dicomRTFileName))
        exit(1)


# -----------------------------------------------------------------------------------------------------------------


def buildPLAN(RTplanFileName, displayInfo=False):
    import sys, os
    import pandas as pd
    import pydicom as dicom
    import numpy as np
    from shutil import copyfile
    
    if displayInfo:
        print('')

    dicomRT=loadDicomRT(RTplanFileName, displayInfo=displayInfo)
    planInfo={                          
              'RTplanFileName': RTplanFileName}

    # get field number and types
    NumberOfFields = dicomRT.FractionGroupSequence[0].NumberOfBeams
    NumberOfFieldsTreatment=0
    NumberOfFieldsSetup=0
    NumberOfFieldsOther=0
    for ifield in range(NumberOfFields):
        if dicomRT.IonBeamSequence[ifield].TreatmentDeliveryType=='TREATMENT':
            NumberOfFieldsTreatment+=1
            planInfo['TreatmentMachineName']=dicomRT.IonBeamSequence[ifield].TreatmentMachineName
        elif dicomRT.IonBeamSequence[ifield].TreatmentDeliveryType=='SETUP':
            NumberOfFieldsSetup+=1
        else:
            NumberOfFieldsOther+=1
    # get number of fractions        
    NumberOfFractionsPlanned = dicomRT.FractionGroupSequence[0].NumberOfFractionsPlanned

    #FractionID
    FractionID = dicomRT.FractionGroupSequence[0].FractionGroupNumber

    if displayInfo:
        print('# Building plan:')
        print('# Nominal gantry: {:s}'.format(planInfo['TreatmentMachineName']))
        print('# Number of fields: {:d} treatment, {:d} setup,  {:d} other'.format(NumberOfFieldsTreatment, NumberOfFieldsSetup, NumberOfFieldsOther))
        print('# Number of fraction planned: {:d}'.format(NumberOfFractionsPlanned))

    planInfo['pencilBeamNumberAll']=0
    planInfo['protonNumberAll']=0
    planInfo['planTotalEnergy']=0
    planInfo['fieldIDList']=[]
    planInfo['fieldNameList']=[]
    planInfo['PlanName']=dicomRT.RTPlanLabel

    # planInfo['NumberOfFractions']=fractionNo
    planInfo['NumberOfFields']=NumberOfFields
    planInfo['FractionID']=FractionID

    fieldsInfo=[]
    for ifield in range(NumberOfFields):
        # convert if is a TREATMENT field (not SETUP or OTHER)
        if dicomRT.IonBeamSequence[ifield].TreatmentDeliveryType!='TREATMENT':
            continue

        if displayInfo:
            print('Field {:d} / {:d}:'.format(ifield+1,NumberOfFieldsTreatment))

        IonBeamSequence = dicomRT.IonBeamSequence[ifield]

        # check if RadiationType is 'PROTON'
        if IonBeamSequence.RadiationType != 'PROTON':
            print('\t\tWarning: No PROTON plan. Field skipped.')
            # npbVec.append(0)
            continue

        # check if the field is MODULATED which means that some properties are the same for all pb
        if IonBeamSequence.ScanMode != 'MODULATED' or IonBeamSequence.BeamType != 'STATIC':
            print('Error this RTPLAN is not for Raster Scan')
            continue

        # check cumulative MU
        for ReferencedBeamSequence in dicomRT.FractionGroupSequence[0].ReferencedBeamSequence:
            if ReferencedBeamSequence.ReferencedBeamNumber==IonBeamSequence.BeamNumber:
                BeamMeterset=ReferencedBeamSequence.BeamMeterset
                break
            else:
                BeamMeterset=np.nan
        if np.isnan(BeamMeterset):
            print('Error: cannot find BeamMeterset for beam number ID {:d}'.format(IonBeamSequence.BeamNumber))
            exit(1)

        # check presence of Range Shifter and Range Shifter ID
        if IonBeamSequence.NumberOfRangeShifters==1:
            if 'RangeShifterSequence' in IonBeamSequence:
                RangeShifterID=IonBeamSequence.RangeShifterSequence[0].RangeShifterID
        else:
            RangeShifterID=''

        # check distance of scanning magnets to isocentre
        for LateralSpreadingDeviceSettingsSequence in IonBeamSequence.IonControlPointSequence[0].LateralSpreadingDeviceSettingsSequence:
            for LateralSpreadingDeviceSequence in IonBeamSequence.LateralSpreadingDeviceSequence:
                if LateralSpreadingDeviceSettingsSequence.ReferencedLateralSpreadingDeviceNumber==LateralSpreadingDeviceSequence.LateralSpreadingDeviceNumber:
                    if LateralSpreadingDeviceSequence.LateralSpreadingDeviceID=='MagnetX':
                        MagnetX_mm=LateralSpreadingDeviceSettingsSequence.IsocenterToLateralSpreadingDeviceDistance
                    elif LateralSpreadingDeviceSequence.LateralSpreadingDeviceID=='MagnetY':
                        MagnetY_mm=LateralSpreadingDeviceSettingsSequence.IsocenterToLateralSpreadingDeviceDistance                
        if MagnetX_mm<MagnetY_mm:
            print('Error: script assumes that the first scanning magnet is deflecting along X (MagnetX_mm={:.3f}, MagnetY_mm={:.3f},)'.format(MagnetX_mm,MagnetY_mm))
            exit(1)
        # define distance from field origin to isocentre            
        fieldOriginToIsoDist_mm=MagnetY_mm
        
        # collect field information
        fieldInfo={'BeamNumber': int(IonBeamSequence.BeamNumber),  # check beam number ID (field ID)
                   'BeamName': IonBeamSequence.BeamName,
                   'RadiationType': IonBeamSequence.RadiationType, # type of radiation (PROTON)
                   'BeamType': IonBeamSequence.BeamType, 
                   'ScanMode': IonBeamSequence.ScanMode,
                   'RangeShifterID': RangeShifterID,
                   'NumberOfControlPoints': int(IonBeamSequence.NumberOfControlPoints), # number of control points (2 x energy slice number) in field
                   'energySliceNumber': int(IonBeamSequence.NumberOfControlPoints/2), # number of energy slices in field
                   'MUsum': float(BeamMeterset), # cumulative MU for field
                   'WeightSum': float(IonBeamSequence.FinalCumulativeMetersetWeight), # cumulative weights for field (the sum of spot weights is used instead)
                   'GantryAngle_deg': float(IonBeamSequence.IonControlPointSequence[0].GantryAngle),
                   'CouchAngle_deg': float(IonBeamSequence.IonControlPointSequence[0].PatientSupportAngle),
                   'fieldOriginToIsoDist_mm': fieldOriginToIsoDist_mm, # distance of field origin point from Iso
                   'SnoutPosition': 368.5, # value taken from TPS eclipse (SnoutPosition 430 means 368.5 mm distance from Iso to RS surface) #IonBeamSequence.IonControlPointSequence[0].SnoutPosition, # snaut position
                   'IsocenterPosition_mm': np.array(IonBeamSequence.IonControlPointSequence[0].IsocenterPosition), # isocentre position
                   'VaccumExitToISODist_mm': 1680
                  }

        if displayInfo:
            print('\tBeam Type: {:s}'.format(fieldInfo['RadiationType']))
            print('\tBeam number ID = {:d}'.format(fieldInfo['BeamNumber']))
            print('\tBeam name = \'{:s}\''.format(fieldInfo['BeamName']))
            print('\tNumber of energy  = {:d}'.format(fieldInfo['energySliceNumber']))
            print('\tCumulative MU = {:.4f}'.format(fieldInfo['MUsum']))
            print('\tRange Shifter = {:s}'.format('None' if not fieldInfo['RangeShifterID'] else fieldInfo['RangeShifterID']))
            print('\tGantry Angle = {:.2f}°'.format(fieldInfo['GantryAngle_deg']))
            print('\tCouch Angle = {:.2f}°'.format(fieldInfo['CouchAngle_deg']))
            print('\tField Origin To Iso Dist = {:.2f} mm'.format(fieldInfo['fieldOriginToIsoDist_mm']))
            print('\tSnout position = {:.2f} mm'.format(fieldInfo['SnoutPosition']))            
            print('\tIso position Vec. =              [ {:.2f} {:.2f} {:.2f}] mm'.format(fieldInfo['IsocenterPosition_mm'][0],fieldInfo['IsocenterPosition_mm'][1],fieldInfo['IsocenterPosition_mm'][2]))

        fieldInfo['slicesInfo']=[]
        for IonControlPointSequence in IonBeamSequence.IonControlPointSequence:
            if np.array(IonControlPointSequence.ScanSpotMetersetWeights).sum()==0:
                continue
            NumberOfScanSpotPositions=int(IonControlPointSequence.NumberOfScanSpotPositions)
             
            sliceInfo={ 'NominalBeamEnergy': float(IonControlPointSequence.NominalBeamEnergy),
                        'NumberOfPaintings': int(IonControlPointSequence.NumberOfPaintings),
                        'ScanSpotTuneID': float(IonControlPointSequence.ScanSpotTuneID)}


            # calculate spots parameters
            sliceSpots=pd.DataFrame({'ScanSpotWeight': IonControlPointSequence.ScanSpotMetersetWeights,
                                     'ScanSpotPositionX_mm': IonControlPointSequence.ScanSpotPositionMap[0::2],
                                     'ScanSpotPositionY_mm': IonControlPointSequence.ScanSpotPositionMap[1::2]})
            sliceSpots['ScanSpotMU'] = sliceSpots.ScanSpotWeight/fieldInfo['WeightSum']*fieldInfo['MUsum']

            sliceInfo['PencilBeams']=len(sliceSpots.index)

            sliceInfo['sliceSpots']=sliceSpots

            fieldInfo['slicesInfo'].append(sliceInfo)

        fieldsInfo.append(fieldInfo)

        fieldInfo['PencilBeams']=0
        fieldInfo['PencilBeamsProtonNumber']=0
        fieldInfo['energyDelivered']=0
        for sliceInfo in fieldInfo['slicesInfo']:
            fieldInfo['PencilBeams']+=sliceInfo['PencilBeams']
            # fieldInfo['PencilBeamsProtonNumber']+=sliceInfo['sliceSpots'].FRED_protonNumber.sum()
            # fieldInfo['energyDelivered']=sliceInfo['sliceSpots'].FRED_energyDelivered.sum()

        planInfo['pencilBeamNumberAll']+=fieldInfo['PencilBeams']
        planInfo['protonNumberAll']+=fieldInfo['PencilBeamsProtonNumber']
        planInfo['planTotalEnergy']+=fieldInfo['energyDelivered']
        planInfo['fieldIDList'].append(fieldInfo['BeamNumber'])
        planInfo['fieldNameList'].append(fieldInfo['BeamName'])

        if displayInfo:
            print('\tPencil Beams loaded = {:d}'.format(fieldInfo['PencilBeams']))
    if displayInfo:
        print('# Total pencil beams loaded in the plan = {:d}'.format(planInfo['pencilBeamNumberAll']))
        print('# Total energy delivered by the plan = {:.20E}'.format(planInfo['planTotalEnergy']))
        
    return planInfo, fieldsInfo



# ----------------------------------------------------------------- PETPLAN ----------------------------------------------





def petPlan_write( localDir, confPath, prodSimPath, patientPath, CTpatientFileName, robustXYZ, geometryName ):
    # read x,y,z dims from MHD file
    [MHDdims, MHDspacing,MHDoffset] = readMHDfile( CTpatientFileName )
    # sourcePosition = [-MHDdims[0]*MHDspacing[0]/2, -MHDdims[1]*MHDspacing[1]/2, -MHDdims[2]*MHDspacing[2]/2]    

    RTplanFileName = os.path.join( patientPath, "TPS", "RN.dcm" )

    planInfo, fieldsInfo = buildPLAN(RTplanFileName, displayInfo=False)
    fieldIdx = len( planInfo['fieldIDList'] ) -1


    # sourcePosition = [-MHDdims[0]*MHDspacing[0]/2 - (MHDoffset[0]- fieldsInfo[fieldIdx]['IsocenterPosition_mm'][0]), -MHDdims[1]*MHDspacing[1]/2 - (MHDoffset[1]-fieldsInfo[fieldIdx]['IsocenterPosition_mm'][1]), -MHDdims[2]*MHDspacing[2]/2 - (MHDoffset[2]-fieldsInfo[fieldIdx]['IsocenterPosition_mm'][2]))]
    sourcePosition = [-MHDdims[0]*MHDspacing[0]/2 , -MHDdims[1]*MHDspacing[1]/2 , (MHDoffset[2]-fieldsInfo[fieldIdx]['IsocenterPosition_mm'][2])]


    if not os.path.exists(os.path.abspath(os.path.join(localDir, 'Results'))):
        os.makedirs(os.path.abspath(os.path.join(localDir, 'Results')))

    #mainSim.mac
    mainSim_fileName = os.path.abspath(os.path.join(localDir, "patient_PET.mac"))
    mainSim_file = open(mainSim_fileName,'w')
    print('#=====================================================', file=mainSim_file)
    print('# SIM PARAMETERS', file=mainSim_file)
    print('#=====================================================', file=mainSim_file)
    print('', file=mainSim_file)   # empty line 
    print('#=====================================================', file=mainSim_file)
    print('# VISUALISATION', file=mainSim_file)
    print('#=====================================================', file=mainSim_file)
    print('', file=mainSim_file)   # empty line 
    tempPath = os.path.abspath(os.path.join(confPath, "Visualisation", "visu_Disable.mac"))
    #tempPath = os.path.abspath(os.path.join(confPath, "Visualisation", "PETvisualisation.mac"))    
    print('/control/execute {:s}'.format(tempPath), file=mainSim_file)    
    # print('/vis/disable', file=mainSim_file)
    print('', file=mainSim_file)   # empty line 
    print('', file=mainSim_file)   # empty line 
    print('#=====================================================', file=mainSim_file)
    print('# GEOMETRY', file=mainSim_file)
    print('#=====================================================', file=mainSim_file)
    print('', file=mainSim_file)   # empty line 
    tempPath = os.path.abspath(os.path.join(confPath, "GateMaterials.db"))
    print('/gate/geometry/setMaterialDatabase    {:s}'.format(tempPath), file=mainSim_file)

    print('', file=mainSim_file)   # empty line 
    #local of from ConfigPath ??????
    # tempPath = os.path.abspath(os.path.join(confPath, "Geometry",  geometryName + ".mac"))    
    # print('/control/execute {:s}'.format(tempPath), file=mainSim_file)
    tempPath= os.path.abspath(os.path.join(localDir, "Geometry.mac"))
    print('/control/execute {:s}'.format(tempPath), file=mainSim_file)

    # print('#/control/execute fullRing.mac', file=mainSim_file)
    # print('/control/execute ../Geometry/Modular_dual_simple.mac ', file=mainSim_file)
    # print('#/control/execute ../Geometry/Modular_24_simple.mac', file=mainSim_file)
    print('#=====================================================', file=mainSim_file)
    print('# P H A N T O M', file=mainSim_file)
    print('#=====================================================', file=mainSim_file)
    print('', file=mainSim_file)   # empty line 
    tempPath= os.path.abspath(os.path.join(localDir, "PETvoxphantom.mac"))
    print('/control/execute {:s}'.format(tempPath), file=mainSim_file)
    # print('#/control/execute waterPH.mac', file=mainSim_file)
    # print('#/control/execute patient_box.mac', file=mainSim_file)
    print('#=====================================================', file=mainSim_file)
    print('# PHYSICS', file=mainSim_file)
    print('#=====================================================', file=mainSim_file)
    print('', file=mainSim_file)   # empty line 
    #tempPath = os.path.abspath(os.path.join(confPath, "Physics",  "PETphysics.mac"))
    tempPath = os.path.abspath(os.path.join(localDir, "PETphysics.mac"))
    print('/control/execute {:s}'.format(tempPath), file=mainSim_file)
    
    print('', file=mainSim_file)   # empty line 
    print('#=====================================================', file=mainSim_file)
    print('#	I N I T I A L I Z E ', file=mainSim_file)
    print('#=====================================================', file=mainSim_file)
    print('', file=mainSim_file)   # empty line 
    print('/gate/run/initialize', file=mainSim_file)
    print('', file=mainSim_file)   # empty line 
    print('#=====================================================', file=mainSim_file)
    print('#	DIGITIZER', file=mainSim_file)
    print('#=====================================================', file=mainSim_file)
    print('', file=mainSim_file)   # empty line 
    tempPath = os.path.abspath(os.path.join(confPath, "Digitizer",  "jpet_proton_beam.mac"))
    print('/control/execute {:s}'.format(tempPath), file=mainSim_file)
    # print('#/control/execute PETdigitizer.mac', file=mainSim_file)
    # print('/control/execute ../Digitizer/proton_beam.mac', file=mainSim_file)
    print('#=====================================================', file=mainSim_file)
    print('#	S O U R C E', file=mainSim_file)
    print('#=====================================================', file=mainSim_file)
    print('', file=mainSim_file)   # empty line
    tempPath= os.path.abspath(os.path.join(localDir, "PatientSource_PET.mac"))
    print('/control/execute    {:s}'.format(tempPath), file=mainSim_file)

    print('', file=mainSim_file)   # empty line 
    print('# =====================================================', file=mainSim_file)
    print('# 	DATA OUTPUT', file=mainSim_file)
    print('# =====================================================', file=mainSim_file)
    print('', file=mainSim_file)   # empty line 

    print('#ROOT Output format', file=mainSim_file)
    print('/gate/output/root/enable', file=mainSim_file)
    print('/gate/output/root/setFileName {:s}'.format(geometryName), file=mainSim_file)
    print('/gate/output/root/setRootHitFlag 0', file=mainSim_file)
    print('/gate/output/root/setRootSinglesFlag 1', file=mainSim_file)
    print('/gate/output/root/setRootCoincidencesFlag 1 ', file=mainSim_file)
    print('', file=mainSim_file)   # empty line 
    # print('#/gate/output/tree/enable', file=mainSim_file)
    # print('#/gate/output/tree/addFileName /tmp/p.npy', file=mainSim_file)
    # print('#/gate/output/tree/addFileName /tmp/p.root', file=mainSim_file)
    # print('#/gate/output/tree/hits/enable', file=mainSim_file)
    # print('#/gate/output/tree/addCollection Coincidences', file=mainSim_file)
    print('', file=mainSim_file)   # empty line 
    print('/gate/application/setTimeSlice     120  s', file=mainSim_file)
    print('/gate/application/setTimeStart     0.0  s', file=mainSim_file)
    print('/gate/application/setTimeStop      120  s', file=mainSim_file)
    print('', file=mainSim_file)   # empty line 
    print('#DOSE OUTPUT', file=mainSim_file)
    print('', file=mainSim_file)   # empty line 
    print('#/gate/patient/addOutput doseoutput', file=mainSim_file)
    print('#/gate/output/doseoutput/saveUncertainty false ', file=mainSim_file)
    print('#/gate/output/doseoutput/setFileName dosematrix.bin ', file=mainSim_file)

    print('', file=mainSim_file)   # empty line 
    print('#=====================================================', file=mainSim_file)
    print('# S T A R T  the A C Q U I S I T I O N', file=mainSim_file)
    print('#=====================================================', file=mainSim_file)
    print('/gate/application/startDAQ', file=mainSim_file)
    mainSim_file.close()

    # -------------------------------------PatientSource_PET.mac
    patientSource_fileName = os.path.abspath(os.path.join(localDir, "PatientSource_PET.mac"))
    patientSource_file = open(patientSource_fileName,'w')
    print('#=====================================================', file=patientSource_file)
    print('#	S O U R C E  ', file=patientSource_file)
    print('#=====================================================', file=patientSource_file)
    print('# Includes distributions of isotopes created in a patient', file=patientSource_file)

    tempPathO15 = os.path.abspath(os.path.join(localDir, "interO15.mac"))
    tempPathO14 = os.path.abspath(os.path.join(localDir, "interO14.mac"))
    tempPathN13 = os.path.abspath(os.path.join(localDir, "interN13.mac"))
    tempPathC11 = os.path.abspath(os.path.join(localDir, "interC11.mac"))
    tempPathC10 = os.path.abspath(os.path.join(localDir, "interC10.mac"))
    tempPathP30 = os.path.abspath(os.path.join(localDir, "interP30.mac"))
    tempPathK38 = os.path.abspath(os.path.join(localDir, "interK38.mac"))

    print('/control/execute {:s}'.format(tempPathO15), file=patientSource_file)
    print('/control/execute {:s}'.format(tempPathO14), file=patientSource_file)
    print('/control/execute {:s}'.format(tempPathN13), file=patientSource_file)
    print('/control/execute {:s}'.format(tempPathC11), file=patientSource_file)
    print('/control/execute {:s}'.format(tempPathC10), file=patientSource_file)
    print('/control/execute {:s}'.format(tempPathP30), file=patientSource_file)
    print('/control/execute {:s}'.format(tempPathK38), file=patientSource_file)
    
#    print('/control/execute interO15.mac       # inter file = patient #', file=patientSource_file)
#    print('/control/execute interN13.mac', file=patientSource_file)
#    print('/control/execute interC11.mac', file=patientSource_file)
#    print('/control/execute interO14.mac', file=patientSource_file)
#    print('/control/execute interC10.mac', file=patientSource_file)
#    print('/control/execute interK38.mac', file=patientSource_file)
#    print('/control/execute interP30.mac', file=patientSource_file)
    print('/gate/source/list', file=patientSource_file)
    patientSource_file.close()


# -------------------------------------   interO14.mac
    interO14_fileName = os.path.abspath(os.path.join(localDir, "interO14.mac"))
    interO14_file = open(interO14_fileName,'w')

    print('/gate/source/addSource voxelO14 voxel', file=interO14_file)
    print('/gate/source/voxelO14/reader/insert image', file=interO14_file)

    print('/gate/source/voxelO14/imageReader/translator/insert linear', file=interO14_file)
    print('/gate/source/voxelO14/imageReader/linearTranslator/setScale 1 Bq', file=interO14_file)
    tempPath = os.path.abspath(os.path.join(prodSimPath, "O14_total_act.h33"))
    # print('/gate/source/voxelO14/imageReader/readFile headerO14.h33', file=interO14_file)
    print('/gate/source/voxelO14/imageReader/readFile {:s}'.format( tempPath ), file=interO14_file)
    print('/gate/source/voxelO14/imageReader/verbose 1', file=interO14_file)

    print('# set the position of the source, such that the source is correctly positioned in the patient CT', file=interO14_file)
    print('# NB! this command sets the bottom left corner of the image!', file=interO14_file)
    print('# so we translate the source -sizex/2 -sizey/2 -sizez/2 mm for the center of the image to overlap with the coordinate origin', file=interO14_file)
    # print('/gate/source/voxelO14/setPosition -175. -127.5 -138.25 mm', file=interO14_file)
    print('/gate/source/voxelO14/setPosition {:.2f} {:.2f} {:.2f} mm'.format(sourcePosition[0], sourcePosition[1], sourcePosition[2]), file=interO14_file)


    print('# The following lines characterize the size', file=interO14_file)
    print('# no difference with an analytical source ', file=interO14_file)
    print('/gate/source/voxelO14/setType gps', file=interO14_file)
    print('/gate/source/voxelO14/gps/particle e+', file=interO14_file)
    print('/gate/source/voxelO14/setForcedUnstableFlag true', file=interO14_file)
    print('/gate/source/voxelO14/setForcedHalfLife 70.598 s ', file=interO14_file)
    print('/gate/source/voxelO14/gps/ene/type Arb', file=interO14_file)
    print('/gate/source/voxelO14/gps/hist/type arb', file=interO14_file)
    print('/gate/source/voxelO14/gps/hist/point 0.10305 0.0601725', file=interO14_file)
    print('/gate/source/voxelO14/gps/hist/point 0.3092 0.132906', file=interO14_file)
    print('/gate/source/voxelO14/gps/hist/point 0.51535 0.175516', file=interO14_file)
    print('/gate/source/voxelO14/gps/hist/point 0.7215 0.189495', file=interO14_file)
    print('/gate/source/voxelO14/gps/hist/point 0.92765 0.175249', file=interO14_file)
    print('/gate/source/voxelO14/gps/hist/point 1.1338 0.13751', file=interO14_file)
    print('/gate/source/voxelO14/gps/hist/point 1.33995 0.0858229', file=interO14_file)
    print('/gate/source/voxelO14/gps/hist/point 1.5461 0.0348707', file=interO14_file)
    print('/gate/source/voxelO14/gps/hist/point 1.75225 0.00454354', file=interO14_file)
    print('/gate/source/voxelO14/gps/hist/point 1.95835 0.000533185', file=interO14_file)
    print('/gate/source/voxelO14/gps/hist/point 2.1645 0.000505896', file=interO14_file)
    print('/gate/source/voxelO14/gps/hist/point 2.37065 0.000472327', file=interO14_file)
    print('/gate/source/voxelO14/gps/hist/point 2.5768 0.000423285', file=interO14_file)
    print('/gate/source/voxelO14/gps/hist/point 2.78295 0.000362572', file=interO14_file)
    print('/gate/source/voxelO14/gps/hist/point 2.9891 0.000293647', file=interO14_file)
    print('/gate/source/voxelO14/gps/hist/point 3.19525 0.000220826', file=interO14_file)
    print('/gate/source/voxelO14/gps/hist/point 3.4014 0.000148096', file=interO14_file)
    print('/gate/source/voxelO14/gps/hist/point 3.60755 8.49921e-05', file=interO14_file)
    print('/gate/source/voxelO14/gps/hist/point 3.8137 3.48888e-05', file=interO14_file)
    print('/gate/source/voxelO14/gps/hist/point 4.01985 5.31821e-06', file=interO14_file)
    print('/gate/source/voxelO14/gps/hist/point 4.01985 5.31821e-06', file=interO14_file)
    print('/gate/source/voxelO14/gps/hist/inter Lin', file=interO14_file)
    print('/gate/source/voxelO14/gps/type Volume', file=interO14_file)
    print('/gate/source/voxelO14/gps/angtype iso ', file=interO14_file)

    interO14_file.close()



# -------------------------------------   interO15.mac
    interO15_fileName = os.path.abspath(os.path.join(localDir, "interO15.mac"))
    interO15_file = open(interO15_fileName,'w')

    print('/gate/source/addSource voxelO15 voxel', file=interO15_file)
    print('/gate/source/voxelO15/reader/insert image', file=interO15_file)

    print('/gate/source/voxelO15/imageReader/translator/insert linear', file=interO15_file)
    print('/gate/source/voxelO15/imageReader/linearTranslator/setScale 1 Bq', file=interO15_file)

    tempPath = os.path.abspath(os.path.join(prodSimPath, "O15_total_act.h33"))
    print('/gate/source/voxelO15/imageReader/readFile {:s}'.format( tempPath ), file=interO15_file)
    # print('/gate/source/voxelO15/imageReader/readFile /scratch/scratch-ssd/jpet/Simulations_GATE/Sources/G3P052E0_PET/inroom/headerO15.h33', file=interO15_file)
    print('/gate/source/voxelO15/imageReader/verbose 1', file=interO15_file)

    print('# set the position of the source, such that the source is correctly positioned in the patient CT', file=interO15_file)
    print('# NB! this command sets the bottom left corner of the image!', file=interO15_file)
    print('# so we translate source -sizex/2 -sizey/2 -sizey/2 mm for the center of the image to overlap with the coordinate origin', file=interO15_file)
    # print('/gate/source/voxelO15/setPosition -175. -127.5 -138.25 mm', file=interO15_file)
    print('/gate/source/voxelO15/setPosition {:.2f} {:.2f} {:.2f} mm'.format(sourcePosition[0], sourcePosition[1], sourcePosition[2]), file=interO15_file)
    

    print('/gate/source/voxelO15/setType gps', file=interO15_file)
    print('/gate/source/voxelO15/gps/particle e+', file=interO15_file)
    print('/gate/source/voxelO15/setForcedUnstableFlag true', file=interO15_file)
    print('/gate/source/voxelO15/setForcedHalfLife 122.24 s', file=interO15_file)
    print('/gate/source/voxelO15/gps/energytype Oxygen15', file=interO15_file)
    print('/gate/source/voxelO15/gps/type Volume', file=interO15_file)
    print('/gate/source/voxelO15/gps/angtype iso ', file=interO15_file)
    interO15_file.close()


    # -------------------------------------   interP30.mac
    interP30_fileName = os.path.abspath(os.path.join(localDir, "interP30.mac"))
    interP30_file = open(interP30_fileName,'w')
    print('/gate/source/addSource voxelP30 voxel', file=interP30_file)
    print('/gate/source/voxelP30/reader/insert image', file=interP30_file)

    print('/gate/source/voxelP30/imageReader/translator/insert linear', file=interP30_file)
    print('/gate/source/voxelP30/imageReader/linearTranslator/setScale 1 Bq', file=interP30_file)

    tempPath = os.path.abspath(os.path.join(prodSimPath, "P30_total_act.h33"))
    print('/gate/source/voxelP30/imageReader/readFile {:s}'.format( tempPath ), file=interP30_file)
    # print('/gate/source/voxelP30/imageReader/readFile headerP30.h33', file=interP30_file)
    print('/gate/source/voxelP30/imageReader/verbose 1', file=interP30_file)

    print('# set the position of the source, such that the source is correctly positioned in the patient CT', file=interP30_file)
    print('# NB! this command sets the bottom left corner of the image!', file=interP30_file)
    print('# so we translate source -sizex/2 -sizey/2 -sizez/2 mm for the center of the image to overlap with the coordinate origin', file=interP30_file)
    # print('/gate/source/voxelP30/setPosition -175. -127.5 -138.25 mm', file=interP30_file)
    print('/gate/source/voxelP30/setPosition {:.2f} {:.2f} {:.2f} mm'.format(sourcePosition[0], sourcePosition[1], sourcePosition[2]), file=interP30_file)    

    print('/gate/source/voxelP30/setType gps', file=interP30_file)
    print('/gate/source/voxelP30/gps/particle e+', file=interP30_file)
    print('/gate/source/voxelP30/setForcedUnstableFlag true', file=interP30_file)
    print('/gate/source/voxelP30/setForcedHalfLife 149.88 s ', file=interP30_file)
    print('/gate/source/voxelP30/gps/ene/type Arb', file=interP30_file)
    print('/gate/source/voxelP30/gps/hist/type arb', file=interP30_file)
    print('/gate/source/voxelP30/gps/hist/point 0.08015 0.00910152', file=interP30_file)
    print('/gate/source/voxelP30/gps/hist/point 0.2404 0.024718', file=interP30_file)
    print('/gate/source/voxelP30/gps/hist/point 0.40065 0.0386907', file=interP30_file)
    print('/gate/source/voxelP30/gps/hist/point 0.5609 0.0515283', file=interP30_file)
    print('/gate/source/voxelP30/gps/hist/point 0.72115 0.0628475', file=interP30_file)
    print('/gate/source/voxelP30/gps/hist/point 0.8814 0.0722514', file=interP30_file)
    print('/gate/source/voxelP30/gps/hist/point 1.04165 0.0793988', file=interP30_file)
    print('/gate/source/voxelP30/gps/hist/point 1.2019 0.0840016', file=interP30_file)
    print('/gate/source/voxelP30/gps/hist/point 1.36215 0.0858624', file=interP30_file)
    print('/gate/source/voxelP30/gps/hist/point 1.5224 0.0850162', file=interP30_file)
    print('/gate/source/voxelP30/gps/hist/point 1.68265 0.0814421', file=interP30_file)
    print('/gate/source/voxelP30/gps/hist/point 1.8429 0.0754725', file=interP30_file)
    print('/gate/source/voxelP30/gps/hist/point 2.00315 0.0672826', file=interP30_file)
    print('/gate/source/voxelP30/gps/hist/point 2.1634 0.057359', file=interP30_file)
    print('/gate/source/voxelP30/gps/hist/point 2.32365 0.0462647', file=interP30_file)
    print('/gate/source/voxelP30/gps/hist/point 2.4839 0.03467', file=interP30_file)
    print('/gate/source/voxelP30/gps/hist/point 2.64415 0.023356', file=interP30_file)
    print('/gate/source/voxelP30/gps/hist/point 2.8044 0.0132229', file=interP30_file)
    print('/gate/source/voxelP30/gps/hist/point 2.96465 0.00528237', file=interP30_file)
    print('/gate/source/voxelP30/gps/hist/point 3.1249 0.000951513', file=interP30_file)
    print('/gate/source/voxelP30/gps/hist/point 3.1249 0.000951513', file=interP30_file)
    print('/gate/source/voxelP30/gps/hist/inter Lin', file=interP30_file)
    print('/gate/source/voxelP30/gps/type Volume', file=interP30_file)
    print('/gate/source/voxelP30/gps/angtype iso ', file=interP30_file)
    interP30_file.close()




# -------------------------------------   interN13.mac
    interN13_fileName = os.path.abspath(os.path.join(localDir, "interN13.mac"))
    interN13_file = open(interN13_fileName,'w')

    print('/gate/source/addSource voxelN13 voxel', file=interN13_file)
    print('/gate/source/voxelN13/reader/insert image', file=interN13_file)

    print('/gate/source/voxelN13/imageReader/translator/insert linear', file=interN13_file)
    print('/gate/source/voxelN13/imageReader/linearTranslator/setScale 1 Bq', file=interN13_file)

    tempPath = os.path.abspath(os.path.join(prodSimPath, "N13_total_act.h33"))
    print('/gate/source/voxelN13/imageReader/readFile {:s}'.format( tempPath ), file=interN13_file)
    # print('/gate/source/voxelN13/imageReader/readFile headerN13.h33', file=interN13_file)
    print('/gate/source/voxelN13/imageReader/verbose 1', file=interN13_file)

    print('# set the position of the source, such that the source is correctly positioned in the patient CT', file=interN13_file)
    print('# NB! this command sets the bottom left corner of the image!', file=interN13_file)
    print('# so we translate source -sizex/2 -sizey/2 -sizez/2 mm for the center of the image to overlap with the coordinate origin', file=interN13_file)
    # print('/gate/source/voxelN13/setPosition -175. -127.5 -138.25 mm', file=interN13_file)
    print('/gate/source/voxelN13/setPosition {:.2f} {:.2f} {:.2f} mm'.format(sourcePosition[0], sourcePosition[1], sourcePosition[2]), file=interN13_file)    

    print('/gate/source/voxelN13/setType gps', file=interN13_file)
    print('/gate/source/voxelN13/gps/particle e+', file=interN13_file)
    print('/gate/source/voxelN13/setForcedUnstableFlag true', file=interN13_file)
    print('/gate/source/voxelN13/setForcedHalfLife 597.9 s ', file=interN13_file)
    print('/gate/source/voxelN13/gps/ene/type Arb', file=interN13_file)
    print('/gate/source/voxelN13/gps/hist/type arb', file=interN13_file)
    print('/gate/source/voxelN13/gps/hist/point 0.0241936 0.0179594', file=interN13_file)
    print('/gate/source/voxelN13/gps/hist/point 0.08985 0.0400168', file=interN13_file)
    print('/gate/source/voxelN13/gps/hist/point 0.1498 0.0548392', file=interN13_file)
    print('/gate/source/voxelN13/gps/hist/point 0.20975 0.0660331', file=interN13_file)
    print('/gate/source/voxelN13/gps/hist/point 0.26965 0.0742828', file=interN13_file)
    print('/gate/source/voxelN13/gps/hist/point 0.32955 0.079895', file=interN13_file)
    print('/gate/source/voxelN13/gps/hist/point 0.3895 0.0829162', file=interN13_file)
    print('/gate/source/voxelN13/gps/hist/point 0.44945 0.083533', file=interN13_file)
    print('/gate/source/voxelN13/gps/hist/point 0.50935 0.0819092', file=interN13_file)
    print('/gate/source/voxelN13/gps/hist/point 0.56925 0.0781231', file=interN13_file)
    print('/gate/source/voxelN13/gps/hist/point 0.6292 0.0724624', file=interN13_file)
    print('/gate/source/voxelN13/gps/hist/point 0.68915 0.0651818', file=interN13_file)
    print('/gate/source/voxelN13/gps/hist/point 0.74905 0.0565922', file=interN13_file)
    print('/gate/source/voxelN13/gps/hist/point 0.80895 0.0471164', file=interN13_file)
    print('/gate/source/voxelN13/gps/hist/point 0.8689 0.0371908', file=interN13_file)
    print('/gate/source/voxelN13/gps/hist/point 0.92885 0.0273277', file=interN13_file)
    print('/gate/source/voxelN13/gps/hist/point 0.98875 0.0180523', file=interN13_file)
    print('/gate/source/voxelN13/gps/hist/point 1.04865 0.0099524', file=interN13_file)
    print('/gate/source/voxelN13/gps/hist/point 1.1086 0.00413211', file=interN13_file)
    print('/gate/source/voxelN13/gps/hist/point 1.16855 0.000884252', file=interN13_file)
    print('/gate/source/voxelN13/gps/hist/point 1.16855 0.000884252', file=interN13_file)
    print('/gate/source/voxelN13/gps/hist/inter Lin', file=interN13_file)
    print('/gate/source/voxelN13/gps/angtype iso ', file=interN13_file)
    interN13_file.close()        



# -------------------------------------   interC10.mac
    interC10_fileName = os.path.abspath(os.path.join(localDir, "interC10.mac"))
    interC10_file = open(interC10_fileName,'w')
    print('/gate/source/addSource voxelC10 voxel', file=interC10_file)
    print('/gate/source/voxelC10/reader/insert image', file=interC10_file)

    print('/gate/source/voxelC10/imageReader/translator/insert linear', file=interC10_file)
    print('/gate/source/voxelC10/imageReader/linearTranslator/setScale 1 Bq', file=interC10_file)

    tempPath = os.path.abspath(os.path.join(prodSimPath, "C10_total_act.h33"))
    print('/gate/source/voxelC10/imageReader/readFile {:s}'.format( tempPath ), file=interC10_file)
    # print('/gate/source/voxelC10/imageReader/readFile headerC10.h33', file=interC10_file)
    print('/gate/source/voxelC10/imageReader/verbose 1', file=interC10_file)

    print('# set the position of the source, such that the source is correctly positioned in the patient CT', file=interC10_file)
    print('# NB! this command sets the bottom left corner of the image!', file=interC10_file)
    print('# so we translate source -sizex/2 -sizey/2 -sizez/2 mm for the center of the image to overlap with the coordinate origin', file=interC10_file)
    # print('/gate/source/voxelC10/setPosition -175. -127.5 -138.25 mm', file=interC10_file)
    print('/gate/source/voxelC10/setPosition {:.2f} {:.2f} {:.2f} mm'.format(sourcePosition[0], sourcePosition[1], sourcePosition[2]), file=interC10_file)    

    print('/gate/source/voxelC10/setType gps', file=interC10_file)
    print('/gate/source/voxelC10/gps/particle e+', file=interC10_file)
    print('/gate/source/voxelC10/setForcedUnstableFlag true', file=interC10_file)
    print('/gate/source/voxelC10/setForcedHalfLife 19.290 s ', file=interC10_file)
    print('/gate/source/voxelC10/gps/ene/type Arb', file=interC10_file)
    print('/gate/source/voxelC10/gps/hist/type arb', file=interC10_file)
    print('/gate/source/voxelC10/gps/hist/point 0.10305 0.0601725', file=interC10_file)
    print('/gate/source/voxelC10/gps/hist/point 0.3092 0.132906', file=interC10_file)
    print('/gate/source/voxelC10/gps/hist/point 0.51535 0.175516', file=interC10_file)
    print('/gate/source/voxelC10/gps/hist/point 0.7215 0.189495', file=interC10_file)
    print('/gate/source/voxelC10/gps/hist/point 0.92765 0.175249', file=interC10_file)
    print('/gate/source/voxelC10/gps/hist/point 1.1338 0.13751', file=interC10_file)
    print('/gate/source/voxelC10/gps/hist/point 1.33995 0.0858229', file=interC10_file)
    print('/gate/source/voxelC10/gps/hist/point 1.5461 0.0348707', file=interC10_file)
    print('/gate/source/voxelC10/gps/hist/point 1.75225 0.00454354', file=interC10_file)
    print('/gate/source/voxelC10/gps/hist/point 1.95835 0.000533185', file=interC10_file)
    print('/gate/source/voxelC10/gps/hist/point 2.1645 0.000505896', file=interC10_file)
    print('/gate/source/voxelC10/gps/hist/point 2.37065 0.000472327', file=interC10_file)
    print('/gate/source/voxelC10/gps/hist/point 2.5768 0.000423285', file=interC10_file)
    print('/gate/source/voxelC10/gps/hist/point 2.78295 0.000362572', file=interC10_file)
    print('/gate/source/voxelC10/gps/hist/point 2.9891 0.000293647', file=interC10_file)
    print('/gate/source/voxelC10/gps/hist/point 3.19525 0.000220826', file=interC10_file)
    print('/gate/source/voxelC10/gps/hist/point 3.4014 0.000148096', file=interC10_file)
    print('/gate/source/voxelC10/gps/hist/point 3.60755 8.49921e-05', file=interC10_file)
    print('/gate/source/voxelC10/gps/hist/point 3.8137 3.48888e-05', file=interC10_file)
    print('/gate/source/voxelC10/gps/hist/point 4.01985 5.31821e-06', file=interC10_file)
    print('/gate/source/voxelC10/gps/hist/point 4.01985 5.31821e-06', file=interC10_file)
    print('/gate/source/voxelC10/gps/hist/inter Lin', file=interC10_file)
    print('/gate/source/voxelC10/gps/type Volume', file=interC10_file)
    print('/gate/source/voxelC10/gps/angtype iso ', file=interC10_file)
    interC10_file.close()




# -------------------------------------   interC11.mac
    interC11_fileName = os.path.abspath(os.path.join(localDir, "interC11.mac"))
    interC11_file = open(interC11_fileName,'w')
    print('/gate/source/addSource voxelC11 voxel', file=interC11_file)
    print('/gate/source/voxelC11/reader/insert image', file=interC11_file)

    print('/gate/source/voxelC11/imageReader/translator/insert linear', file=interC11_file)
    print('/gate/source/voxelC11/imageReader/linearTranslator/setScale 1 Bq', file=interC11_file)

    tempPath = os.path.abspath(os.path.join(prodSimPath, "C11_total_act.h33"))
    print('/gate/source/voxelC11/imageReader/readFile {:s}'.format( tempPath ), file=interC11_file)
    # print('/gate/source/voxelC11/imageReader/readFile headerC11.h33', file=interC11_file)
    print('/gate/source/voxelC11/imageReader/verbose 1', file=interC11_file)

    print('# set the position of the source, such that the source is correctly positioned in the patient CT', file=interC11_file)
    print('# NB! this command sets the bottom left corner of the image!', file=interC11_file)
    print('# so we translate source -sizex/2 -sizey/2 -sizey-2 mm for the center of the image to overlap with the coordinate origin', file=interC11_file)
    # print('/gate/source/voxelC11/setPosition -175. -127.5 -138.25 mm', file=interC11_file)
    print('/gate/source/voxelC11/setPosition {:.2f} {:.2f} {:.2f} mm'.format(sourcePosition[0], sourcePosition[1], sourcePosition[2]), file=interC11_file)    

    print('/gate/source/voxelC11/setType gps', file=interC11_file)
    print('/gate/source/voxelC11/gps/particle e+', file=interC11_file)
    print('/gate/source/voxelC11/setForcedUnstableFlag true', file=interC11_file)
    print('/gate/source/voxelC11/setForcedHalfLife 1220.04 s ', file=interC11_file)
    print('/gate/source/voxelC11/gps/energytype Carbon11', file=interC11_file)
    print('/gate/source/voxelC11/gps/type Volume', file=interC11_file)
    print('/gate/source/voxelC11/gps/angtype iso ', file=interC11_file)
    interC11_file.close()    


# -------------------------------------   interK38.mac
    interK38_fileName = os.path.abspath(os.path.join(localDir, "interK38.mac"))
    interK38_file = open(interK38_fileName,'w')
    print('/gate/source/addSource voxelK38 voxel', file=interK38_file)
    print('/gate/source/voxelK38/reader/insert image', file=interK38_file)

    print('/gate/source/voxelK38/imageReader/translator/insert linear', file=interK38_file)
    print('/gate/source/voxelK38/imageReader/linearTranslator/setScale 1 Bq', file=interK38_file)

    tempPath = os.path.abspath(os.path.join(prodSimPath, "K38_total_act.h33"))
    print('/gate/source/voxelK38/imageReader/readFile {:s}'.format( tempPath ), file=interK38_file)
    # print('/gate/source/voxelK38/imageReader/readFile headerK38.h33', file=interK38_file)
    print('/gate/source/voxelK38/imageReader/verbose 1', file=interK38_file)

    print('# set the position of the source, such that the source is correctly positioned in the patient CT', file=interK38_file)
    print('# NB! this command sets the bottom left corner of the image!', file=interK38_file)
    print('# so we translate source -sizex/2 -sizey/2 -sizez/2 mm for the center of the image to overlap with the coordinate origin', file=interK38_file)
    # print('/gate/source/voxelK38/setPosition -175. 127.5 -138.25 mm', file=interK38_file)
    print('/gate/source/voxelK38/setPosition {:.2f} {:.2f} {:.2f} mm'.format(sourcePosition[0], sourcePosition[1], sourcePosition[2]), file=interK38_file)    

    print('/gate/source/voxelK38/setType gps', file=interK38_file)
    print('/gate/source/voxelK38/gps/particle e+', file=interK38_file)
    print('/gate/source/voxelK38/setForcedUnstableFlag true', file=interK38_file)
    print('/gate/source/voxelK38/setForcedHalfLife 458.16 s ', file=interK38_file)
    print('/gate/source/voxelK38/gps/ene/type Arb', file=interK38_file)
    print('/gate/source/voxelK38/gps/hist/type arb', file=interK38_file)
    print('/gate/source/voxelK38/gps/hist/point 0.06835 0.00896571', file=interK38_file)
    print('/gate/source/voxelK38/gps/hist/point 0.205 0.0257479', file=interK38_file)
    print('/gate/source/voxelK38/gps/hist/point 0.34165 0.0402954', file=interK38_file)
    print('/gate/source/voxelK38/gps/hist/point 0.47835 0.0532211', file=interK38_file)
    print('/gate/source/voxelK38/gps/hist/point 0.61505 0.0642919', file=interK38_file)
    print('/gate/source/voxelK38/gps/hist/point 0.7517 0.0732261', file=interK38_file)
    print('/gate/source/voxelK38/gps/hist/point 0.88835 0.0798101', file=interK38_file)
    print('/gate/source/voxelK38/gps/hist/point 1.02505 0.0839162', file=interK38_file)
    print('/gate/source/voxelK38/gps/hist/point 1.1617 0.0852569', file=interK38_file)
    print('/gate/source/voxelK38/gps/hist/point 1.29835 0.0840583', file=interK38_file)
    print('/gate/source/voxelK38/gps/hist/point 1.43505 0.0801996', file=interK38_file)
    print('/gate/source/voxelK38/gps/hist/point 1.5717 0.0739922', file=interK38_file)
    print('/gate/source/voxelK38/gps/hist/point 1.70835 0.0656884', file=interK38_file)
    print('/gate/source/voxelK38/gps/hist/point 1.84505 0.0558936', file=interK38_file)
    print('/gate/source/voxelK38/gps/hist/point 1.98175 0.0449362', file=interK38_file)
    print('/gate/source/voxelK38/gps/hist/point 2.1184 0.0335648', file=interK38_file)
    print('/gate/source/voxelK38/gps/hist/point 2.25505 0.022564', file=interK38_file)
    print('/gate/source/voxelK38/gps/hist/point 2.39175 0.0127483', file=interK38_file)
    print('/gate/source/voxelK38/gps/hist/point 2.5284 0.00505686', file=interK38_file)
    print('/gate/source/voxelK38/gps/hist/point 2.66505 0.000916668', file=interK38_file)
    print('/gate/source/voxelK38/gps/hist/point 2.66505 0.000916668', file=interK38_file)
    print('/gate/source/voxelK38/gps/hist/inter Lin', file=interK38_file)
    print('/gate/source/voxelK38/gps/type Volume', file=interK38_file)
    print('/gate/source/voxelK38/gps/angtype iso ', file=interK38_file)
    interK38_file.close()    


    # kopiowanie plików z geometria !!!!!!!!!! GEOMETRY          
    copyfile(os.path.abspath(os.path.join(confPath, "Geometry", geometryName + ".mac")), os.path.abspath(os.path.join(localDir, "Geometry.mac")))            


    # -------------------------------------   PETvoxphantom.mac
    voxphantom_fileName = os.path.abspath(os.path.join(localDir, "PETvoxphantom.mac"))
    voxphantom_file = open(voxphantom_fileName,'w')    

    print('# patient', file=voxphantom_file)
    # print('/gate/cylindricalPET/daughters/name                patient', file=voxphantom_file)
    # print('/gate/cylindricalPET/daughters/insert              ImageNestedParametrisedVolume', file=voxphantom_file)
    print('/gate/world/daughters/name                patient', file=voxphantom_file)
    print('/gate/world/daughters/insert              ImageNestedParametrisedVolume', file=voxphantom_file)

    tempPath = os.path.abspath(os.path.join(confPath, "Materials",  "HUmaterials.db"))
    # print('/gate/geometry/setMaterialDatabase              /scratch/scratch-ssd/jpet/Simulations_GATE/proton_karol_test/HUmaterials.db #patient-HUmaterials.db', file=voxphantom_file)
    print('/gate/geometry/setMaterialDatabase              {:s}'.format( tempPath ), file=voxphantom_file)

    tempPath = os.path.abspath(os.path.join(confPath, "Materials",  "HU2mat.mac"))
    # print('/gate/patient/geometry/setHUToMaterialFile      /scratch/scratch-ssd/jpet/Simulations_GATE/proton_karol_test/HU2mat.mac #patient-HU2mat.txt', file=voxphantom_file)
    print('/gate/patient/geometry/setHUToMaterialFile      {:s}'.format( tempPath ), file=voxphantom_file)
    
    print('/gate/patient/attachPhantomSD', file=voxphantom_file)
    # print('/gate/patient/geometry/setImage                 /scratch/scratch-ssd/jpet/Simulations_GATE/CT_data/G3P052E0_CT_crop.mhd', file=voxphantom_file)
    print('/gate/patient/geometry/setImage                 {:s}'.format( CTpatientFileName ), file=voxphantom_file)

    sourcePosition = [-MHDdims[0]*MHDspacing[0]/2 , -MHDdims[1]*MHDspacing[1]/2 , (MHDoffset[2]-fieldsInfo[fieldIdx]['IsocenterPosition_mm'][2])]
 
    print(fieldsInfo[fieldIdx]['IsocenterPosition_mm'][2])
    CTPosition = [MHDdims[0]*MHDspacing[0]/2+MHDoffset[0] , MHDdims[1]*MHDspacing[1]/2+MHDoffset[1], fieldsInfo[fieldIdx]['IsocenterPosition_mm'][2]]

    # print('/gate/patient/geometry/TranslateTheImageAtThisIsoCenter  {:+.2f}\t{:+.2f}\t{:+.2f}\t mm # put here the isocentre from RN dicom plan'.format(fieldsInfo[fieldIdx]['IsocenterPosition_mm'][0]-robustXYZ[0],fieldsInfo[fieldIdx]['IsocenterPosition_mm'][1]-robustXYZ[1],fieldsInfo[fieldIdx]['IsocenterPosition_mm'][2]-robustXYZ[2]), file=voxphantom_file)         
    # print('/gate/patient/geometry/TranslateTheImageAtThisIsoCenter  {:+.2f}\t{:+.2f}\t{:+.2f}\t mm # put here the isocentre from RN dicom plan'.format(MHDoffset[0] + (MHDdims[0]*MHDspacing[0])/2.0 , MHDoffset[1] + (MHDdims[1]*MHDspacing[1])/2.0, Disocenter_z ), file=voxphantom_file)         
    print('/gate/patient/geometry/TranslateTheImageAtThisIsoCenter  {:+.2f}\t{:+.2f}\t{:+.2f}\t mm # put here the isocentre from RN dicom plan'.format(CTPosition[0], CTPosition[1], CTPosition[2]), file=voxphantom_file)  
    
    # print('/gate/patient/geometry/TranslateTheImageAtThisIsoCenter -18.60       -258.10 -633.50 mm#for G3P052E0_CT_Crop', file=voxphantom_file)

    print('/gate/patient/placement/setRotationAxis           0 1 0', file=voxphantom_file)
    print('/gate/patient/placement/setRotationAngle          -0.0 deg', file=voxphantom_file)



    # -------------------------------------   PETphysics.mac
#    physics_fileName = os.path.abspath(os.path.join(localDir, "PETphysics.mac"))
#    physics_file = open(physics_fileName,'w')
#    print('#=====================================================', file=physics_file)
#    print('# PHYSICS', file=physics_file)
#    print('#=====================================================', file=physics_file)
#    print('/gate/physics/addPhysicsList emstandard', file=physics_file)
#    print('', file=physics_file)

#    print('/gate/physics/addProcess PhotoElectric gamma', file=physics_file)
#    print('/gate/physics/processes/PhotoElectric/setModel StandardModel ', file=physics_file)
#    print('', file=physics_file)

#    print('/gate/physics/addProcess Compton gamma', file=physics_file)
#    print('/gate/physics/processes/Compton/setModel StandardModel ', file=physics_file)
#    print('', file=physics_file)

#    print('/gate/physics/addProcess GammaConversion', file=physics_file)
#    print('/gate/physics/processes/GammaConversion/setModel StandardModel ', file=physics_file)
#    print('', file=physics_file)

#    print('/gate/physics/addProcess ElectronIonisation e+', file=physics_file)
#    print('/gate/physics/addProcess ElectronIonisation e-', file=physics_file)
#    print('/gate/physics/processes/ElectronIonisation/setModel StandardModel e+', file=physics_file)
#    print('/gate/physics/processes/ElectronIonisation/setModel StandardModel e-', file=physics_file)
#    print('', file=physics_file)

#    print('/gate/physics/addProcess Bremsstrahlung e+', file=physics_file)
#    print('/gate/physics/addProcess Bremsstrahlung e-', file=physics_file)
#    print('/gate/physics/processes/Bremsstrahlung/setModel StandardModel e+', file=physics_file)
#    print('/gate/physics/processes/Bremsstrahlung/setModel StandardModel e-', file=physics_file)
#    print('', file=physics_file)

#    print('/gate/physics/addProcess G4PositronAnnihilation e+', file=physics_file)
#    print('/gate/physics/processes/G4PositronAnnihilation/setModel StandardModel ', file=physics_file)
#    print('', file=physics_file)

#    print('/gate/physics/addProcess eMultipleScattering e-', file=physics_file)
#    print('/gate/physics/addProcess eMultipleScattering e+', file=physics_file)
#    physics_file.close()    



    # kopiowanie plików z fizyka !!!!!!!!!! PHYSICS          
    copyfile(os.path.abspath(os.path.join(confPath, "Physics", "PETphysics.mac")), os.path.abspath(os.path.join(localDir, "PETphysics.mac")))            


# -------------------------######## _______________ petPlan _____________ ######## --------------------
# _____________________________________________________________________________________________________


prodSimPath = simPath
confPath = configPath

petPlan_write( localDir, confPath, prodSimPath, patientPath, CTpatientFileName, robustXYZlist, geometryName )

