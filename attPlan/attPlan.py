
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
# import medpy
from multiprocessing import Process
from multiprocessing import Pool
import glob
from scipy.interpolate import RegularGridInterpolator, LinearNDInterpolator



# VOXEL size in mm (scalling factor: mm/pixel)
voxel_x = 2.5
voxel_y = 2.5
voxel_z = 2.5 

# FOV size in mm (needed to calculate the matrix size)
fov_x = 400. 
fov_y = 400. 
fov_z = 500. 
fov_z_dualhead_3x4 = 400. 

# MATRIX size 
mat_x = int(fov_x/voxel_x) 
mat_y = int(fov_y/voxel_y)
mat_z = int(fov_z/voxel_z)
# mat_z_dualhead_3x4 = int(fov_z_dualhead_3x4/voxel_z)



# ---------------------------------------------------------- DEFs ----------------------------------------------------------


######  DEF ######  -------------------------  read DimSize from MHD text file ------------------------------------
def readMHDfile( CTfileMHD ):
    # import SimpleITK as sitk

    MHDimage = sitk.ReadImage(CTfileMHD)        
    dims = MHDimage.GetSize()   # get Dimensions from MHD file
    spacing = MHDimage.GetSpacing()  # get Element Spacing from MHD file
    offset = MHDimage.GetOrigin()  # get offset ?????    
    return [dims, spacing, offset, MHDimage]


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


# ---------------------------------------------------------- DEFs ----------------------------------------------------------


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


## --


def load_and_convertCT_Jakub( CTfilePath ):
    [MHDdims, MHDspacing, MHDoffset,  MHDimage] = readMHDfile( CTfilePath )
    
    # #CT DATA
    ct_voxel_x = MHDspacing[0] 
    ct_voxel_y = MHDspacing[1]
    ct_voxel_z = MHDspacing[2]

    ct_mat_x = MHDdims[0] 
    ct_mat_y = MHDdims[1]
    ct_mat_z = MHDdims[2]

    ct_offset = -1000
    break_point = 47
    
    below_a = 9.6E-5
    below_b = 0.0
    above_a = 5.10E-5
    above_b = 4.71E-2
    
    # # Convert the image to a  numpy array first 
    ct_scan = sitk.GetArrayFromImage(MHDimage)
    ct_scan = np.transpose(ct_scan, (2, 1, 0))
    umap = np.zeros((ct_scan.shape))
    # print(ct_scan.shape)

    # # Jakub

    # # CT voltage tube 120kVp
    # # All factors are taken from the Carney et al., Transforming CT images for attenuation correction in PET/CT, Medical Physics, Vol. 33, No. 4, 2006
    # # Below the breaking point: u = below_a*(HU+1000) + below_b 
    # # Above the breaking point: u = above_a*(HU+1000) + above_b

    umap1 = below_a*(ct_scan-ct_offset)+below_b
    umap2 = above_a*(ct_scan-ct_offset)+above_b
    umap[ct_scan<break_point] = umap1[ct_scan<break_point]
    umap[ct_scan>=break_point] = umap2[ct_scan>=break_point]
    print('... done!')

    return umap



def load_and_convertCT_Karol( CTfilePath ):
    [MHDdims, MHDspacing, MHDoffset, MHDimage] = readMHDfile( CTfilePath )
    
    # #CT DATA
    ct_voxel_x = MHDspacing[0] 
    ct_voxel_y = MHDspacing[1]
    ct_voxel_z = MHDspacing[2]

    ct_mat_x = MHDdims[0] 
    ct_mat_y = MHDdims[1]
    ct_mat_z = MHDdims[2]

    # ct_offset = -1000
    break_point = 0
    
    mu_PET_H2O = 0.096
    mu_PET_Bone = 0.172
    mu_CT_H2O = 0.184
    mu_CT_Bone = 0.428
    F = ( mu_CT_H2O*(mu_PET_Bone-mu_PET_H2O) / (1000*(mu_CT_Bone-mu_CT_H2O)) )

    # # Convert the image to a  numpy array first 
    ct_scan = sitk.GetArrayFromImage(MHDimage)
    ct_scan = np.transpose(ct_scan, (2, 1, 0))
    umap = np.zeros((ct_scan.shape))
    # print(ct_scan.shape)

    # # Karol

    # # CT voltage tube 120kVp
    # # All factors are taken from the Burger EJNM 2002

    umap1 = mu_PET_H2O*((ct_scan+1000))/1000
    umap2 = mu_PET_H2O + (ct_scan)*F
    umap[ct_scan<break_point] = umap1[ct_scan<break_point]
    umap[ct_scan>=break_point] = umap2[ct_scan>=break_point]

    # umap = ct_scan   #for control purposes
    print('... done!')

    # img = sitk.GetImageFromArray(umap)
    # diiims = img.GetSize()

    return umap




def save_rescaled_umap_to_header( fpath, fname, mat_x, mat_y, mat_z, voxel_x, voxel_y, voxel_z ):
    bytes_per_pixel = "4"
    temppathRAW = os.path.abspath(os.path.join(fpath, fname + '.raw'))

    with open(os.path.abspath(os.path.join(fpath,fname + '.hdr')), "w+") as f:
        f.write("!name of data file := {0}\n".format(temppathRAW))
        f.write("!total number of images := 1\n")
        f.write("imagedata byte order := LITTLEENDIAN\n")
        f.write("number of dimensions := 3\n")
        f.write("!matrix size [1] := {0}\n".format(mat_x))
        f.write("!matrix size [2] := {0}\n".format(mat_y))
        f.write("!matrix size [3] := {0}\n".format(mat_z))
        f.write("!number format := float\n")
        f.write("!number of bytes per pixel := {0}\n".format(bytes_per_pixel))
        f.write("scaling factor (mm/pixel) [1] := {0}\n".format(voxel_x))
        f.write("scaling factor (mm/pixel) [2] := {0}\n".format(voxel_y))
        f.write("scaling factor (mm/pixel) [3] := {0}\n".format(voxel_z))
        f.write("image duration (sec) := 1\n")


# -----------------------------


def save_umap_to_binary ( PatientPath, CTfilePath, CTfileName, RTplanFileName, umap, fname):

    [MHDdims, MHDspacing, MHDoffset, MHDimage] = readMHDfile( CTfilePath )

    planInfo, fieldsInfo = buildPLAN(RTplanFileName, displayInfo=False)
    # fieldIdx = len( planInfo['fieldIDList'] ) -1
    fieldIdx = 0    # take the first field


    # #CT DATA
    ct_voxel_x = MHDspacing[0] 
    ct_voxel_y = MHDspacing[1]
    ct_voxel_z = MHDspacing[2]

    ct_mat_x = MHDdims[0] 
    ct_mat_y = MHDdims[1]
    ct_mat_z = MHDdims[2]

    offset_x = MHDoffset[0]
    offset_y = MHDoffset[1]
    offset_z = MHDoffset[2]


    isocenter_x = fieldsInfo[fieldIdx]['IsocenterPosition_mm'][0]
    isocenter_y = fieldsInfo[fieldIdx]['IsocenterPosition_mm'][1]
    isocenter_z = fieldsInfo[fieldIdx]['IsocenterPosition_mm'][2]


    # CALCULATING ORIGINAL CT COORDINATES
    x=np.linspace( offset_x-isocenter_x , offset_x-isocenter_x+(ct_mat_x-1)*ct_voxel_x , ct_mat_x)
    y=np.linspace( offset_y-isocenter_y , offset_y-isocenter_y+(ct_mat_y-1)*ct_voxel_y , ct_mat_y)
    z=np.linspace( offset_z-isocenter_z , offset_z-isocenter_z+(ct_mat_z-1)*ct_voxel_z , ct_mat_z)

    linear_interp = RegularGridInterpolator((x, y, z), umap, bounds_error=False, fill_value=0.0)

    # BARRELS, DUALHEAD_1x6, DUALHEAD_2x6 SETUPS

    print('Processing umap for barrels and two dualheads ...')
    # CALCULATING NEW umap COORDINATES

    pts_normal = np.zeros((mat_x * mat_y * mat_z, 3))
    rescaled_umap_normal = np.zeros((mat_x, mat_y, mat_z))

    for i in range (0, rescaled_umap_normal.shape[2]):
        for j in range (0, rescaled_umap_normal.shape[1]):
            for k in range (0, rescaled_umap_normal.shape[0]):
                ind = k + j*rescaled_umap_normal.shape[0] + i*rescaled_umap_normal.shape[0]*rescaled_umap_normal.shape[1]
                pts_normal[ind,0] = -fov_x/2 + voxel_x*(1 + k) - voxel_x/2
                pts_normal[ind,1] = -fov_y/2 + voxel_y*(1 + j) - voxel_y/2
                pts_normal[ind,2] = -fov_z/2 + voxel_z*(1 + i) - voxel_z/2

    #print pts_normal
    # TRILINEAR INTERPOLATION
    rescaled_umap = linear_interp(pts_normal).reshape(rescaled_umap_normal.shape, order = 'F')

    # SAVING NEW umap
    fnameN = 'umap' + CTfileName + '_' + fname
    temppath = os.path.join( PatientPath, "TPS",  fnameN + '.raw' ) 
    rescaled_umap.flatten('F').astype('float32').tofile( temppath )
    print('... done!')

    fpath = os.path.join( PatientPath, "TPS" ) 
    save_rescaled_umap_to_header( fpath, fnameN, mat_x, mat_y, mat_z, voxel_x, voxel_y, voxel_z )


# ------------------------------


def save_umap_to_binary_CTspace ( PatientPath, CTfilePath, CTfileName, RTplanFileName, umap, fname):

    [MHDdims, MHDspacing, MHDoffset, MHDimage] = readMHDfile( CTfilePath )

    planInfo, fieldsInfo = buildPLAN(RTplanFileName, displayInfo=False)
    # fieldIdx = len( planInfo['fieldIDList'] ) -1
    fieldIdx = 0    # take the first field


    # #CT DATA
    # ct_voxel_x = MHDspacing[0] 
    # ct_voxel_y = MHDspacing[1]
    # ct_voxel_z = MHDspacing[2]

    # ct_mat_x = MHDdims[0] 
    # ct_mat_y = MHDdims[1]
    # ct_mat_z = MHDdims[2]

    mat_x = MHDdims[0] 
    mat_y = MHDdims[1]
    mat_z = MHDdims[2]

    voxel_x = MHDspacing[0] 
    voxel_y = MHDspacing[1]
    voxel_z = MHDspacing[2]


    # offset_x = MHDoffset[0]
    # offset_y = MHDoffset[1]
    # offset_z = MHDoffset[2]


    # isocenter_x = fieldsInfo[fieldIdx]['IsocenterPosition_mm'][0]
    # isocenter_y = fieldsInfo[fieldIdx]['IsocenterPosition_mm'][1]
    # isocenter_z = fieldsInfo[fieldIdx]['IsocenterPosition_mm'][2]


    # # CALCULATING ORIGINAL CT COORDINATES
    # x=np.linspace( offset_x-isocenter_x , offset_x-isocenter_x+(ct_mat_x-1)*ct_voxel_x , ct_mat_x)
    # y=np.linspace( offset_y-isocenter_y , offset_y-isocenter_y+(ct_mat_y-1)*ct_voxel_y , ct_mat_y)
    # z=np.linspace( offset_z-isocenter_z , offset_z-isocenter_z+(ct_mat_z-1)*ct_voxel_z , ct_mat_z)

    # linear_interp = RegularGridInterpolator((x, y, z), umap, bounds_error=False, fill_value=0.0)

    # # BARRELS, DUALHEAD_1x6, DUALHEAD_2x6 SETUPS

    # print('Processing umap for barrels and two dualheads ...')
    # # CALCULATING NEW umap COORDINATES

    # pts_normal = np.zeros((mat_x * mat_y * mat_z, 3))
    # rescaled_umap_normal = np.zeros((mat_x, mat_y, mat_z))

    # for i in range (0, rescaled_umap_normal.shape[2]):
    #     for j in range (0, rescaled_umap_normal.shape[1]):
    #         for k in range (0, rescaled_umap_normal.shape[0]):
    #             ind = k + j*rescaled_umap_normal.shape[0] + i*rescaled_umap_normal.shape[0]*rescaled_umap_normal.shape[1]
    #             pts_normal[ind,0] = -fov_x/2 + voxel_x*(1 + k) - voxel_x/2
    #             pts_normal[ind,1] = -fov_y/2 + voxel_y*(1 + j) - voxel_y/2
    #             pts_normal[ind,2] = -fov_z/2 + voxel_z*(1 + i) - voxel_z/2

    # #print pts_normal
    # # TRILINEAR INTERPOLATION
    # rescaled_umap = linear_interp(pts_normal).reshape(rescaled_umap_normal.shape, order = 'F')

    # SAVING NEW umap
    fnameN = 'umap_CT_' + CTfileName + '_' + fname
    temppath = os.path.join( PatientPath, "TPS",  fnameN + '.raw' ) 
    umap.flatten('F').astype('float32').tofile( temppath )
    print('...CTspace done!')

    fpath = os.path.join( PatientPath, "TPS" ) 
    save_rescaled_umap_to_header( fpath, fnameN, mat_x, mat_y, mat_z, voxel_x, voxel_y, voxel_z )


# ------------------------------

def save_umap_to_binary_CTspace_interpolated ( PatientPath, CTfilePath, CTfileName, RTplanFileName, umap, fname):

    [MHDdims, MHDspacing, MHDoffset, MHDimage] = readMHDfile( CTfilePath )

    planInfo, fieldsInfo = buildPLAN(RTplanFileName, displayInfo=False)
    # fieldIdx = len( planInfo['fieldIDList'] ) -1
    fieldIdx = 0    # take the first field


    # #CT DATA
    ct_voxel_x = MHDspacing[0] 
    ct_voxel_y = MHDspacing[1]
    ct_voxel_z = MHDspacing[2]

    ct_mat_x = MHDdims[0] 
    ct_mat_y = MHDdims[1]
    ct_mat_z = MHDdims[2]

    mat_x = MHDdims[0] 
    mat_y = MHDdims[1]
    mat_z = MHDdims[2]

    voxel_x = MHDspacing[0] 
    voxel_y = MHDspacing[1]
    voxel_z = MHDspacing[2]

    offset_x = MHDoffset[0]
    offset_y = MHDoffset[1]
    offset_z = MHDoffset[2]


    # isocenter_x = fieldsInfo[fieldIdx]['IsocenterPosition_mm'][0]
    isocenter_x = 0
    # isocenter_y = fieldsInfo[fieldIdx]['IsocenterPosition_mm'][1]
    isocenter_y = 0
    isocenter_z = fieldsInfo[fieldIdx]['IsocenterPosition_mm'][2]


    # CALCULATING ORIGINAL CT COORDINATES
    x=np.linspace( offset_x-isocenter_x , offset_x-isocenter_x+(ct_mat_x-1)*ct_voxel_x , ct_mat_x)
    y=np.linspace( offset_y-isocenter_y , offset_y-isocenter_y+(ct_mat_y-1)*ct_voxel_y , ct_mat_y)
    z=np.linspace( offset_z-isocenter_z , offset_z-isocenter_z+(ct_mat_z-1)*ct_voxel_z , ct_mat_z)

    linear_interp = RegularGridInterpolator((x, y, z), umap, bounds_error=False, fill_value=0.0)

    # BARRELS, DUALHEAD_1x6, DUALHEAD_2x6 SETUPS

    print('Processing umap for barrels and two dualheads ...')
    # CALCULATING NEW umap COORDINATES

    pts_normal = np.zeros((mat_x * mat_y * mat_z, 3))
    rescaled_umap_normal = np.zeros((mat_x, mat_y, mat_z))

    for i in range (0, rescaled_umap_normal.shape[2]):
        for j in range (0, rescaled_umap_normal.shape[1]):
            for k in range (0, rescaled_umap_normal.shape[0]):
                ind = k + j*rescaled_umap_normal.shape[0] + i*rescaled_umap_normal.shape[0]*rescaled_umap_normal.shape[1]
                pts_normal[ind,0] = -fov_x/2 + voxel_x*(1 + k) - voxel_x/2
                pts_normal[ind,1] = -fov_y/2 + voxel_y*(1 + j) - voxel_y/2
                pts_normal[ind,2] = -fov_z/2 + voxel_z*(1 + i) - voxel_z/2

    #print pts_normal
    # TRILINEAR INTERPOLATION
    rescaled_umap = linear_interp(pts_normal).reshape(rescaled_umap_normal.shape, order = 'F')

    # SAVING NEW umap
    fnameN = 'umap_CT_' + CTfileName + '_' + fname
    temppath = os.path.join( PatientPath, "TPS",  fnameN + '.raw' ) 
    umap.flatten('F').astype('float32').tofile( temppath )
    print('...CTspace done!')

    fpath = os.path.join( PatientPath, "TPS" ) 
    save_rescaled_umap_to_header( fpath, fnameN, mat_x, mat_y, mat_z, voxel_x, voxel_y, voxel_z )


# ------------------------------




def save_umap_to_binary_CTspace_interpolatedDB ( PatientPath, CTfilePath, CTfileName, RTplanFileName, umap, fname):

    [MHDdims, MHDspacing, MHDoffset, MHDimage] = readMHDfile( CTfilePath )

    planInfo, fieldsInfo = buildPLAN(RTplanFileName, displayInfo=False)
    # fieldIdx = len( planInfo['fieldIDList'] ) -1
    fieldIdx = 0    # take the first field


    # #CT DATA
    ct_voxel_x = MHDspacing[0] 
    ct_voxel_y = MHDspacing[1]
    ct_voxel_z = MHDspacing[2]

    ct_mat_x = MHDdims[0] 
    ct_mat_y = MHDdims[1]
    ct_mat_z = MHDdims[2]

    mat_x = MHDdims[0] 
    mat_y = MHDdims[1]
    mat_z = MHDdims[2]

    voxel_x = MHDspacing[0] 
    voxel_y = MHDspacing[1]
    voxel_z = MHDspacing[2]

    offset_x = MHDoffset[0]
    offset_y = MHDoffset[1]
    offset_z = MHDoffset[2]


    # isocenter_x = fieldsInfo[fieldIdx]['IsocenterPosition_mm'][0]
    isocenter_x = 0
    # isocenter_y = fieldsInfo[fieldIdx]['IsocenterPosition_mm'][1]
    isocenter_y = 0
    isocenter_z = fieldsInfo[fieldIdx]['IsocenterPosition_mm'][2]
    Disocenter_z = offset_z - isocenter_z + (ct_mat_z/2)*ct_voxel_z

    # CALCULATING ORIGINAL CT COORDINATES
    # x0=np.linspace( offset_x , offset_x+(ct_mat_x-1)*ct_voxel_x , ct_mat_x)
    # y0=np.linspace( offset_y , offset_y+(ct_mat_y-1)*ct_voxel_y , ct_mat_y)
    # z0=np.linspace( offset_z , offset_z+(ct_mat_z-1)*ct_voxel_z , ct_mat_z)

    # x0=np.linspace( -(ct_mat_x/2.0)*ct_voxel_x , -(ct_mat_x/2.0)*ct_voxel_x+(ct_mat_x-1)*ct_voxel_x , ct_mat_x)
    # y0=np.linspace( -(ct_mat_y/2.0)*ct_voxel_y , -(ct_mat_y/2.0)*ct_voxel_y+(ct_mat_y-1)*ct_voxel_y , ct_mat_y)
    # z0=np.linspace( -(ct_mat_z/2.0)*ct_voxel_z , -(ct_mat_z/2.0)*ct_voxel_z+(ct_mat_z-1)*ct_voxel_z , ct_mat_z)

    # x=np.linspace( -(ct_mat_x/2.0)*ct_voxel_x , -(ct_mat_x/2.0)*ct_voxel_x+(ct_mat_x-1)*ct_voxel_x , ct_mat_x)
    # y=np.linspace( -(ct_mat_y/2.0)*ct_voxel_y , -(ct_mat_y/2.0)*ct_voxel_y+(ct_mat_y-1)*ct_voxel_y , ct_mat_y)
    # z=np.linspace( -(ct_mat_z/2.0)*ct_voxel_z-Disocenter_z, -(ct_mat_z/2.0)*ct_voxel_z-Disocenter_z+(ct_mat_z-1)*ct_voxel_z , ct_mat_z)




    x0=np.linspace( offset_x*0.5 , offset_x*0.5+(ct_mat_x-1)*ct_voxel_x , ct_mat_x)
    y0=np.linspace( offset_y*0.5 , offset_y*0.5+(ct_mat_y-1)*ct_voxel_y , ct_mat_y)
    z0=np.linspace( offset_z*0.5 , offset_z*0.5+(ct_mat_z-1)*ct_voxel_z , ct_mat_z)

    x=np.linspace( offset_x*0.5 , offset_x*0.5+(ct_mat_x-1)*ct_voxel_x , ct_mat_x)
    y=np.linspace( offset_y*0.5 , offset_y*0.5+(ct_mat_y-1)*ct_voxel_y , ct_mat_y)
    z=np.linspace( offset_z*0.5-Disocenter_z , offset_z*0.5+(ct_mat_z-1)*ct_voxel_z-Disocenter_z , ct_mat_z)




    # x0=np.linspace( 0 , (ct_mat_x-1) , ct_mat_x)
    # y0=np.linspace( 0 , (ct_mat_y-1) , ct_mat_y)
    # z0=np.linspace( 0 , (ct_mat_z-1) , ct_mat_z)

    # x=np.linspace( 0 , (ct_mat_x-1) , ct_mat_x)
    # y=np.linspace( 0 , (ct_mat_y-1) , ct_mat_y)
    # z=np.linspace( 0-(Disocenter_z/ct_voxel_z)*0.5 , (ct_mat_z-1)-(Disocenter_z/ct_voxel_z)*0.5 , ct_mat_z)



    # x=np.linspace( offset_x-isocenter_x , offset_x-isocenter_x+(ct_mat_x-1)*ct_voxel_x , ct_mat_x)
    # y=np.linspace( offset_y-isocenter_y , offset_y-isocenter_y+(ct_mat_y-1)*ct_voxel_y , ct_mat_y)
    # z=np.linspace( offset_z-isocenter_z , offset_z-isocenter_z+(ct_mat_z-1)*ct_voxel_z , ct_mat_z)

    # CALCULATING NEW umap COORDINATES
    # x=np.linspace( offset_x , offset_x+(ct_mat_x-1)*ct_voxel_x , ct_mat_x)
    # y=np.linspace( offset_y , offset_y+(ct_mat_y-1)*ct_voxel_y , ct_mat_y)
    # # x=np.linspace( isocenter_x , isocenter_x+(ct_mat_x-1)*ct_voxel_x , ct_mat_x)
    # # y=np.linspace( isocenter_y , isocenter_y+(ct_mat_y-1)*ct_voxel_y , ct_mat_y)
    # z=np.linspace( isocenter_z , isocenter_z+(ct_mat_z-1)*ct_voxel_z , ct_mat_z)


    # x = np.linspace( -fov_x / 2 + ct_voxel_x / 2, -fov_x / 2 + ct_voxel_x / 2 + (ct_mat_x - 1) * ct_voxel_x, ct_mat_x)
    # y = np.linspace( -fov_y / 2 + ct_voxel_y / 2, -fov_y / 2 + ct_voxel_y / 2 + (ct_mat_y - 1) * ct_voxel_y, ct_mat_y  )
    # z = np.linspace( -fov_z / 2 + ct_voxel_z / 2, -fov_z / 2 + ct_voxel_z / 2 + (ct_mat_z - 1) * ct_voxel_z, ct_mat_z)



    # xx0, yy0, zz0 = np.meshgrid(x0, y0, z0)
    xx, yy, zz = np.meshgrid(x, y, z)
    xx = np.transpose(xx, (1, 0, 2))
    yy = np.transpose(yy, (1, 0, 2))
    zz = np.transpose(zz, (1, 0, 2))

    # linear_interp = RegularGridInterpolator((xo, yo, zo), umap, bounds_error=False, fill_value=0.0)
    linear_interp = RegularGridInterpolator((x0, y0, z0), umap, bounds_error=False, fill_value=0.0)

    # TRILINEAR INTERPOLATION
    # rescaled_umap = linear_interp(pts_normal).reshape(rescaled_umap_normal.shape, order = 'F')
    rescaled_umap = linear_interp((xx, yy, zz))

    # SAVING NEW umap
    fnameN = 'umap_CT_' + CTfileName + '_' + fname
    temppath = os.path.join( PatientPath, "TPS", fnameN + '.raw' ) 
    rescaled_umap.flatten('F').astype('float32').tofile( temppath )
    print('...CTspace done!')

    fpath = os.path.join( PatientPath, "TPS" ) 
    save_rescaled_umap_to_header( fpath, fnameN, mat_x, mat_y, mat_z, voxel_x, voxel_y, voxel_z )




################ ----------------------------- new version -------------------------------------------------------


def load_and_convertCT_Jakub_2img( CTfilePath ):
    [MHDdims, MHDspacing, MHDoffset,  MHDimage] = readMHDfile( CTfilePath )
    
    # #CT DATA
    ct_voxel_x = MHDspacing[0] 
    ct_voxel_y = MHDspacing[1]
    ct_voxel_z = MHDspacing[2]

    ct_mat_x = MHDdims[0] 
    ct_mat_y = MHDdims[1]
    ct_mat_z = MHDdims[2]

    ct_offset = -1000
    break_point = 47
    
    below_a = 9.6E-5
    below_b = 0.0
    above_a = 5.10E-5
    above_b = 4.71E-2
    
    # # Convert the image to a  numpy array first 
    ct_scan = sitk.GetArrayFromImage(MHDimage)
    # ct_scan = np.transpose(ct_scan, (2, 1, 0))
    umap = np.zeros((ct_scan.shape))
    # print(ct_scan.shape)

    # # Jakub

    # # CT voltage tube 120kVp
    # # All factors are taken from the Carney et al., Transforming CT images for attenuation correction in PET/CT, Medical Physics, Vol. 33, No. 4, 2006
    # # Below the breaking point: u = below_a*(HU+1000) + below_b 
    # # Above the breaking point: u = above_a*(HU+1000) + above_b

    umap1 = below_a*(ct_scan-ct_offset)+below_b
    umap2 = above_a*(ct_scan-ct_offset)+above_b
    umap[ct_scan<break_point] = umap1[ct_scan<break_point]
    umap[ct_scan>=break_point] = umap2[ct_scan>=break_point]
    print('... done!')

    umap_img = sitk.GetImageFromArray(umap)
    umap_img.SetOrigin( MHDimage.GetOrigin() )
    umap_img.SetSpacing( MHDimage.GetSpacing() )
    
    return umap_img



def load_and_convertCT_Karol_2img( CTfilePath ):
    [MHDdims, MHDspacing, MHDoffset, MHDimage] = readMHDfile( CTfilePath )
    
    # #CT DATA
    ct_voxel_x = MHDspacing[0] 
    ct_voxel_y = MHDspacing[1]
    ct_voxel_z = MHDspacing[2]

    ct_mat_x = MHDdims[0] 
    ct_mat_y = MHDdims[1]
    ct_mat_z = MHDdims[2]

    # ct_offset = -1000
    break_point = 0
    
    mu_PET_H2O = 0.096
    mu_PET_Bone = 0.172
    mu_CT_H2O = 0.184
    mu_CT_Bone = 0.428
    F = ( mu_CT_H2O*(mu_PET_Bone-mu_PET_H2O) / (1000*(mu_CT_Bone-mu_CT_H2O)) )

    # # Convert the image to a  numpy array first 
    ct_scan = sitk.GetArrayFromImage(MHDimage)
    # ct_scan = np.transpose(ct_scan, (2, 1, 0))
    umap = np.zeros((ct_scan.shape))
    # print(ct_scan.shape)

    # # Karol

    # # CT voltage tube 120kVp
    # # All factors are taken from the Burger EJNM 2002

    umap1 = mu_PET_H2O*((ct_scan+1000))/1000
    umap2 = mu_PET_H2O + (ct_scan)*F
    umap[ct_scan<break_point] = umap1[ct_scan<break_point]
    umap[ct_scan>=break_point] = umap2[ct_scan>=break_point]

    # umap = ct_scan   #for control purposes
    print('... done!')

    umap_img = sitk.GetImageFromArray(umap)
    umap_img.SetOrigin( MHDimage.GetOrigin() )
    umap_img.SetSpacing( MHDimage.GetSpacing() )
    
    return umap_img



def save_umap_to_MHD_CTspace( PatientPath, CTfilePath, CTfileName, RTplanFileName, umap_img, fname):

    [MHDdims, MHDspacing, MHDoffset, MHDimage] = readMHDfile( CTfilePath )

    planInfo, fieldsInfo = buildPLAN(RTplanFileName, displayInfo=False)
    # fieldIdx = len( planInfo['fieldIDList'] ) -1
    fieldIdx = 0    # take the first field


    # # #CT DATA
    # ct_voxel_x = MHDspacing[0] 
    # ct_voxel_y = MHDspacing[1]
    # ct_voxel_z = MHDspacing[2]

    # ct_mat_x = MHDdims[0] 
    # ct_mat_y = MHDdims[1]
    # ct_mat_z = MHDdims[2]

    # mat_x = MHDdims[0] 
    # mat_y = MHDdims[1]
    # mat_z = MHDdims[2]

    # voxel_x = MHDspacing[0] 
    # voxel_y = MHDspacing[1]
    # voxel_z = MHDspacing[2]


    # offset_x = MHDoffset[0]
    # offset_y = MHDoffset[1]
    # offset_z = MHDoffset[2]

    # isocenter_x = fieldsInfo[fieldIdx]['IsocenterPosition_mm'][0]
    # isocenter_y = fieldsInfo[fieldIdx]['IsocenterPosition_mm'][1]
    # isocenter_z = fieldsInfo[fieldIdx]['IsocenterPosition_mm'][2]

    # SAVING NEW umap_img

    fnameN = 'umap_CT_' + CTfileName + '_' + fname
    temppath = os.path.join( PatientPath, "TPS",  fnameN + '.raw' ) 
    # umap.flatten('F').astype('float32').tofile( temppath )

    outputImageFileName = os.path.join( PatientPath, "TPS",  fnameN + '.mhd' ) 
    writer = sitk.ImageFileWriter()
    writer.SetFileName(outputImageFileName)
    writer.Execute(umap_img)
    print('...CTspace done!')

    # fpath = os.path.join( PatientPath, "TPS" ) 
    # save_rescaled_umap_to_header( fpath, fnameN, mat_x, mat_y, mat_z, voxel_x, voxel_y, voxel_z )


### --------------------------------------------------------------------------------------------------------------
### --------------------------------------------------------------------------------------------------------------
### --------------------------------------------------------------------------------------------------------------
### --------------------------------------------------------------------------------------------------------------

# parse input
inputParse = optparse.OptionParser()
inputParse.add_option('--PatientPath', '-p', default=None)
inputParse.add_option('--CTfile', '-t', default=None)
# inputParse.add_option('--simLabel', '-l', default="Sim0100")
# inputParse.add_option('--configPath', '-c', default=None)
options, fileList = inputParse.parse_args()

# get Patient DICOM RTplan PATH
if options.PatientPath == None:
    print("Missing DICOM RT PLAN file path !")
    # exit(0)
else:
    PatientPath = os.path.abspath( options.PatientPath )

# PATIENTfile="/Users/dborys/Documents/POLSL/2020_IFJ_postdoc/getPlanTests/PETscripts/compare/myscript/production_maps/G3P052E0"
print( PatientPath )

# get CT .mhd file name
# for example:  CT_ExtROICrop_2.0x2.0x2.0.mhd 
if options.CTfile == None:
    print("Missing MHD file path !")
    # exit(0)
else:
    CTfileList = os.path.join( PatientPath, "TPS", options.CTfile )    
    CTfileName = options.CTfile
    CTfileName = CTfileName[:-4]

# CTfileList = "/Users/dborys/Documents/POLSL/2020_IFJ_postdoc/getPlanTests/patientData/G3P052E0/TPS/CT_ExtROICrop_2.5x2.5x2.5.mhd"


CTfilePath = os.path.join( PatientPath, "TPS", CTfileList ) 
RTfilePath = os.path.join( PatientPath, "TPS", "RN.dcm" ) 
# read x,y,z dims from MHD file
# [MHDdims, MHDspacing] = readMHDfile( CTfilePath )



# ------------------------------------------------------------------------------------------------------------


umapJakub = load_and_convertCT_Jakub( CTfileList )
umapKarol = load_and_convertCT_Karol( CTfileList )

# save_umap_to_binary( PatientPath, CTfileList, CTfileName, RTfilePath, umapJakub, 'J')
# save_umap_to_binary( PatientPath, CTfileList, CTfileName, RTfilePath, umapKarol, 'K')

# save_umap_to_binary_CTspace( PatientPath, CTfileList, CTfileName, RTfilePath, umapJakub, 'J')
# save_umap_to_binary_CTspace( PatientPath, CTfileList, CTfileName, RTfilePath, umapKarol, 'K')


save_umap_to_binary_CTspace_interpolatedDB( PatientPath, CTfileList, CTfileName, RTfilePath, umapJakub, 'J')
save_umap_to_binary_CTspace_interpolatedDB( PatientPath, CTfileList, CTfileName, RTfilePath, umapKarol, 'K')

# umap_imgJakub = load_and_convertCT_Jakub_2img( CTfileList )
# umap_imgKarol = load_and_convertCT_Karol_2img( CTfileList )


# save_umap_to_MHD_CTspace( PatientPath, CTfileList, CTfileName, RTfilePath, umap_imgJakub, 'J')
# save_umap_to_MHD_CTspace( PatientPath, CTfileList, CTfileName, RTfilePath, umap_imgKarol, 'K')