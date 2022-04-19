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


def buildPLAN(RTplanFileName, beamModelFileName, RSModelFileName, materialDefinitionFileName, fractionNo, beamModelInterpolationMethod='slinear', displayInfo=False):
    import sys, os
    import pandas as pd
    import pydicom as dicom
    import numpy as np
    from shutil import copyfile
    
    if displayInfo:
        print('')

    dicomRT=loadDicomRT(RTplanFileName, displayInfo=displayInfo)
    planInfo={'beamModelFileName': beamModelFileName,
              'RSModelFileName': RSModelFileName,
              'materialDefinitionFileName': materialDefinitionFileName,
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

    planInfo['NumberOfFractions']=fractionNo
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
            # beamModelEnergy=interpolateBeamModel(beamModel, float(IonControlPointSequence.NominalBeamEnergy), method=beamModelInterpolationMethod)
            
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

def mergePlanInfoAndFieldsInfo(planInfoList, fieldsInfoList):
    # merge planInfoList to planInfo
    pencilBeamNumberAll=0
    protonNumberAll=0
    planTotalEnergy=0
    RTplanFileNames=[]
    fieldIDList=[]
    for planInfo in planInfoList:
        pencilBeamNumberAll+=planInfo['pencilBeamNumberAll']
        protonNumberAll+=planInfo['protonNumberAll']
        planTotalEnergy+=planInfo['planTotalEnergy']
        fieldIDList.extend(planInfo['fieldIDList'])
        RTplanFileNames.append(planInfo['RTplanFileName'])
        
    planInfo={'beamModelFileName': planInfoList[0]['beamModelFileName'],
              'RSModelFileName': planInfoList[0]['RSModelFileName'],
              'materialDefinitionFileName': planInfoList[0]['materialDefinitionFileName'],
              'RTplanFileName': RTplanFileNames,
              'TreatmentMachineName': planInfoList[0]['TreatmentMachineName'],
              'pencilBeamNumberAll': pencilBeamNumberAll,
              'protonNumberAll': protonNumberAll,
              'planTotalEnergy': planTotalEnergy,
              'fieldIDList': list(range(1,len(fieldIDList)+1)),
              'fieldNameList': ['fieldMulti {:d}'.format(fieldID) for fieldID in list(range(1,len(fieldIDList)+1))]}

    # merge fieldsInfoList to fieldsInfo
    fieldsInfo=[]
    BeamNumber=1
    for fieldsInfoSingle in fieldsInfoList:
        for fieldInfo in fieldsInfoSingle:
            fieldInfo['BeamNumber']=BeamNumber
            fieldInfo['BeamName']='fieldMulti {:d}'.format(BeamNumber)
            BeamNumber+=1
        fieldsInfo.extend(fieldsInfoSingle)
    return planInfo, fieldsInfo

# ----------------------------------- DEF: ------------ calculateGATEparticleNum ---------------------------------------

def calculateGATEparticleNum( csv_filePath, cores, Nparticles, ParticleUnit ):
    #read info from fieldStat.csv file loacted in rtplanfileName/TPS path
    import pandas as pd
    import sys, os
    stats_file = os.path.join( csv_filePath, "TPS", "fieldStat.csv")

    df = pd.read_csv( stats_file, delim_whitespace=True )

    # how many fields ?
    fieldsNumber = df.shape[0]    
    PencilBeams = []
    ProtonNo = []

    for fieldIdx in range(0, fieldsNumber):
        PencilBeams.append(int(df.loc[fieldIdx,'PencilBeams']))
        # print(PencilBeams[fieldIdx])
        ProtonNo.append(int(df.loc[fieldIdx,'ProtonNo']))
        print(ProtonNo[fieldIdx])

    # calculate final number of protons including PencilBeams OR ProtonNo AND dividing by Number of cores/CPUs    
    GateProtons = []
    if ParticleUnit == "primPr":
        for fieldIdx in range( 0, len(PencilBeams) ):
            GateProtons.append( int( ( Nparticles*ProtonNo[fieldIdx])/cores ) )  
    else:
        for fieldIdx in range( 0, len(PencilBeams) ):
            GateProtons.append( int( ( Nparticles*PencilBeams[fieldIdx])/cores ) )

    # return the LIST with Number of Protons for GATE (DIVIDED by Number of CORES - READY for GATE job splitter)
    # list for each field in order of USAGE  
    return GateProtons


def writeRtplanGate( rtplanfileName, CTpatientFileName, simLabel, GATEconfigPath, cores, actors, Nparticles, ParticleUnit, robustXYZ, planInfo, fieldsInfo, getPlanVersion='', fieldToSave=None, displayInfo=True):
    import sys, os, glob
    import pandas as pd
    import pydicom as dicom
    import numpy as np
    from shutil import copyfile
    # Xsec='emittance'

    rtplanfileName=os.path.abspath(rtplanfileName)

    pathPatient =  os.path.dirname( os.path.dirname( CTpatientFileName ) )
    # now we should obtain PatientPath
    simPath = os.path.abspath( os.path.join( pathPatient, "gate", simLabel ) )

    # currDir = os.getcwd()
    currDir = pathPatient
    if not os.path.exists(os.path.join(currDir, 'gate')):
        os.makedirs(os.path.join(currDir, 'gate'))

    #create subdirectory in patient/gate  with simLabel :  Sim0001 etc
    if not os.path.exists( simPath ):
        os.makedirs( simPath )

    # list of primaries numbers - for each field (with the order of field usage (fieldIdx))
    NBparts = calculateGATEparticleNum( pathPatient, cores, Nparticles, ParticleUnit )

    # should be given at the command line, if not $PWD is used
    confPath = GATEconfigPath

    ######  DEF ######   read DimSize from MHD text file
    def readMHDfile( CTfileMHD ):
        import SimpleITK as sitk

        MHDimage = sitk.ReadImage(CTfileMHD)        
        dims = MHDimage.GetSize()   # get Dimensions from MHD file
        spacing = MHDimage.GetSpacing()  # get Element Spacing from MHD file

        return [dims, spacing]

    # read x,y,z dims from MHD file
    [MHDdims, MHDspacing] = readMHDfile( CTpatientFileName )

    
    for fieldIdx in range(0,len(fieldsInfo)):
        #create new subfolder 
        #subfolder name
        # example : 1_Field2
        sub_name_field = str(fieldIdx+1) + "_Field" + str(fieldsInfo[fieldIdx]['BeamNumber'])
        if not os.path.exists(os.path.abspath(os.path.join( simPath, sub_name_field ))):
            os.makedirs(os.path.join( simPath, sub_name_field))
        
        localDir = os.path.abspath(os.path.join( simPath, sub_name_field))

        if not os.path.exists(os.path.abspath(os.path.join(localDir, 'Results'))):
            os.makedirs(os.path.abspath(os.path.join(localDir, 'Results')))

        # create new Filename with Field ID added
        # rtplanfileName_Field = localDir + "/" + rtplanFname + "_field" + str(fieldsInfo[fieldIdx]['BeamNumber']) + rtplanFnameExt; 
        rtplanfileName_Field = os.path.abspath(os.path.join(localDir, "rtplan.mac"))

        rtplan_file=open(rtplanfileName_Field,'w')

        # write information to a .mac file
        #GENERAL PLAN INFO SECTION !!!!
        print('#TREATMENT-PLAN-DESCRIPTION', file=rtplan_file)

        print('#PlanName', file=rtplan_file)
        #print('{:s}'.format(planInfo['PlanName']), file=rtplan_file)
        # RTplanFileName   instead of   PlanName   because it was EMPTY in some cases
        print('{:s}'.format(planInfo['RTplanFileName']), file=rtplan_file)

        print('#NumberOfFractions', file=rtplan_file)
        print('{:d}'.format(1), file=rtplan_file)   # THERE CAN BE ONLY ONE
        # print('{:d}'.format(planInfo['NumberOfFractions']), file=rtplan_file)

        print('##FractionID', file=rtplan_file)
        print('{:d}'.format(planInfo['FractionID']), file=rtplan_file)
        print('##NumberOfFields', file=rtplan_file)
        print('{:d}'.format(1), file=rtplan_file)  # THERE CAN BE ONLY ONE
        # print('{:d}'.format(planInfo['NumberOfFields']), file=rtplan_file)

        # Only current FieldsID
        print('##FieldsID', file=rtplan_file)
        print('{:d}'.format(fieldsInfo[fieldIdx]['BeamNumber']), file=rtplan_file)
        
        # calculate (sum) TotalMetersetWeightOfAllFields
        # UWAGA : NIEZGODNE Z PRZYKLADEM !!!!!!!!!!!!!!!!!!!!! 
        TotalMetersetWeightOfAllFields=0.0
        for fieldIdx_ in range(0,len(fieldsInfo)):
            weight = fieldsInfo[fieldIdx_]['WeightSum']
            TotalMetersetWeightOfAllFields = TotalMetersetWeightOfAllFields + weight

        print('##TotalMetersetWeightOfAllFields', file=rtplan_file)
        print('{:f}'.format(TotalMetersetWeightOfAllFields), file=rtplan_file) 

        CumulativeMetersetWeight = 0.0

        #SPECIFIC FIELD INFO SECTION !!!!
    
        # print(' ', file=rtplan_file) # Empty line
        print('#FIELD-DESCRIPTION', file=rtplan_file)
        print('###FieldID', file=rtplan_file)
        print('{:d}'.format(fieldsInfo[fieldIdx]['BeamNumber']), file=rtplan_file)


        print('###FinalCumulativeMeterSetWeight', file=rtplan_file)
        print('{:f}'.format(fieldsInfo[fieldIdx]['WeightSum']), file=rtplan_file)

        print('###GantryAngle', file=rtplan_file)
        print('{:+.1f}'.format(fieldsInfo[fieldIdx]['GantryAngle_deg']), file=rtplan_file)

        print('###PatientSupportAngle', file=rtplan_file)
        print('{:+.1f}'.format(fieldsInfo[fieldIdx]['CouchAngle_deg']), file=rtplan_file)

        print('###IsocenterPosition', file=rtplan_file)
        print('{:+.1f}\t{:+.1f}\t{:+.1f}\t'.format(fieldsInfo[fieldIdx]['IsocenterPosition_mm'][0],fieldsInfo[fieldIdx]['IsocenterPosition_mm'][1],fieldsInfo[fieldIdx]['IsocenterPosition_mm'][2]), file=rtplan_file)

        print('###NumberOfControlPoints', file=rtplan_file)
        # print('{:d}'.format(fieldsInfo[fieldIdx]['NumberOfControlPoints']), file=rtplan_file)
        print('{:d}'.format(int(fieldsInfo[fieldIdx]['NumberOfControlPoints']/2)), file=rtplan_file)
        
        # print(' ', file=rtplan_file) # Empty line
        print('#SPOTS-DESCRIPTION', file=rtplan_file)
        for controlptsIdx in range(0, len(fieldsInfo[fieldIdx]['slicesInfo'])):
            print('####ControlPointIndex', file=rtplan_file)            
            print('{:d}'.format(controlptsIdx), file=rtplan_file)
            print('####SpotTunnedID', file=rtplan_file)
            print('{:.1f}'.format(fieldsInfo[fieldIdx]['slicesInfo'][controlptsIdx]['ScanSpotTuneID']), file=rtplan_file)            

            print('####CumulativeMetersetWeight', file=rtplan_file)
            print('{:f}'.format(CumulativeMetersetWeight), file=rtplan_file)

            print('####Energy (MeV)', file=rtplan_file)
            print('{:f}'.format(fieldsInfo[fieldIdx]['slicesInfo'][controlptsIdx]['NominalBeamEnergy']), file=rtplan_file)

            print('####NbOfScannedSpots', file=rtplan_file)
            print('{:d}'.format(fieldsInfo[fieldIdx]['slicesInfo'][controlptsIdx]['PencilBeams']), file=rtplan_file)

            sliceSpots = fieldsInfo[fieldIdx]['slicesInfo'][controlptsIdx]['sliceSpots']
            
            print('####X Y Weight', file=rtplan_file)
            for _,sliceSpot in sliceSpots.iterrows():
                spot_x = sliceSpot['ScanSpotPositionX_mm']
                spot_y = sliceSpot['ScanSpotPositionY_mm']
                spot_weight = sliceSpot['ScanSpotWeight']
                print('{:+.1f}\t{:+.1f}\t{:f}\t'.format(spot_x, spot_y, spot_weight), file=rtplan_file)
                CumulativeMetersetWeight = CumulativeMetersetWeight + spot_weight

        rtplan_file.close()

        # get machineName : GTR3 or GTR4
        TreatmentMachineName = planInfo['TreatmentMachineName'] # getGantryName(RTplanFileName)

        RSname = TreatmentMachineName + "_RangeShifter.mac"
        RSfileN = os.path.abspath(os.path.join(confPath, "Geometry", RSname))
        

        ######  DEF ###### read RS translation
        def readRangeShifterTranslation( RSfile ):
            # Using readlines() 
            RSfileName = open(RSfile, 'r') 
            Lines = RSfileName.readlines() 

            RStranslation = 0.0        
            for line in Lines: 
                if line.find("gate/rangeshifter/placement/setTranslation") >= 0:
                    # print( line )
                    #remove multiple spaces from line
                    line_cleared = " ".join(line.split())
                    # separete by space
                    separated = line_cleared.split(sep=" ")
                    RStranslation = float( separated[3] ) 
            RSfileName.close()
            return RStranslation        

        RStranslation = readRangeShifterTranslation( RSfileN )

        RSfileName = os.path.abspath(os.path.join(localDir , "RS.mac"))
        RS_file=open(RSfileName,'w')
        print('# leave this file empty if no RS is used', file=RS_file)

        #check if there is 'RS_Block' name in RangeShifterID field
        if len(fieldsInfo[fieldIdx]['RangeShifterID']) > 0: 
            # if RS is used             
            # print('/control/execute {:s}Geometry/{:s}_RangeShifter.mac # choose GTR3 or GTR4'.format(confPath,TreatmentMachineName), file=RS_file)
            tempPath = os.path.abspath(os.path.join(confPath, "Geometry", TreatmentMachineName + "_RangeShifter.mac"))
            print('/control/execute {:s} # choose GTR3 or GTR4'.format(tempPath), file=RS_file)
            print('/gate/rangeshifter/placement/setRotationAxis         0 0 1', file=RS_file)
            print('/gate/rangeshifter/placement/setRotationAngle        {:+.1f} deg # set here the gantry rotation'.format(fieldsInfo[fieldIdx]['GantryAngle_deg']), file=RS_file)        
            # print('/gate/rangeshifter/placement/setTranslation          0 -386.825 0 mm # check in GTR3_RangeShifter/GTR4_RangeShifter for the proper value. It is fixed for GTR3 and GTR4.', file=RS_file)
            print('/gate/rangeshifter/placement/setTranslation          0 {:+f} 0 mm # check in GTR3_RangeShifter/GTR4_RangeShifter for the proper value. It is fixed for GTR3 and GTR4.'.format(RStranslation), file=RS_file)            
            print('/control/add  rangeshifterTranslationAngle           -90   {:+.1f} # change the second value (90) by the gantry rotation'.format(fieldsInfo[fieldIdx]['GantryAngle_deg']), file=RS_file)
            print('/gate/rangeshifter/placement/setPhiOfTranslation     {rangeshifterTranslationAngle} deg', file=RS_file)
        RS_file.close()

        source_fileName = os.path.abspath(os.path.join(localDir, "source.mac"))
        source_file=open(source_fileName,'w')   
        print('# beam model and plan description', file=source_file)
        print('/gate/source/addSource PBS                                      TPSPencilBeam', file=source_file)
        print('/gate/source/PBS/setParticleType                                proton', file=source_file)
        tempPath = os.path.abspath(os.path.join(localDir ,"rtplan.mac"))
        print('/gate/source/PBS/setPlan                                        {:s}'.format(tempPath), file=source_file)
        print('/gate/source/PBS/setFlatGenerationFlag                          false', file=source_file)
        print('/gate/source/PBS/setSpotIntensityAsNbIons                       false', file=source_file)
        print('/gate/source/PBS/setSigmaEnergyInMeVFlag                        true', file=source_file)
        print('/gate/source/PBS/setSortedSpotGenerationFlag                    false', file=source_file)
        print('/gate/source/PBS/setTestFlag                                    false', file=source_file)
        print('/gate/source/PBS/setBeamConvergence                             true', file=source_file)
        tempPath = os.path.abspath(os.path.join(confPath, "Sources", "beamModel", TreatmentMachineName + "_20190919_beamModel_sourceDescription.mac"))
        # print('/gate/source/PBS/setSourceDescriptionFile                       {:s}Sources/BeamModel/{:s}_20190919_beamModel_sourceDescription.mac # change this according to beam model for GTR3 and GTR4'.format(confPath,TreatmentMachineName), file=source_file)
        print('/gate/source/PBS/setSourceDescriptionFile                       {:s}'.format(tempPath), file=source_file)        
        source_file.close()

        #doseActor.mac
        doseActor_fileName = os.path.abspath(os.path.join(localDir, "Actors.mac"))
        doseActor_file=open(doseActor_fileName,'w')
        
        if "doseAll" in actors:
            print('# Dose actor for the patientCT', file=doseActor_file)
            print('# doseAll', file=doseActor_file)
            print('/gate/actor/addActor                    DoseActor dose', file=doseActor_file)
            print('/gate/actor/dose/attachTo               patientCT', file=doseActor_file)
            print('/gate/actor/dose/stepHitType            random', file=doseActor_file)        
            # print('/gate/actor/dose/setResolution          175 144 121 # put here the pixel resolution from CT mhd.', file=doseActor_file)
            print('/gate/actor/dose/setResolution          {:d} {:d} {:d} # put here the pixel resolution from CT mhd.'.format(MHDdims[0],MHDdims[1],MHDdims[2]), file=doseActor_file)        
            print('/gate/actor/dose/enableEdep             false', file=doseActor_file)
            print('/gate/actor/dose/enableUncertaintyEdep  false', file=doseActor_file)
            print('/gate/actor/dose/enableSquaredEdep      false', file=doseActor_file)
            print('/gate/actor/dose/enableDose             true', file=doseActor_file)
            print('/gate/actor/dose/enableUncertaintyDose  false', file=doseActor_file)
            print('/gate/actor/dose/enableSquaredDose      false', file=doseActor_file)
            tempPath = os.path.join("Results", "dose.mhd")
            # print('/gate/actor/dose/save                   Results/dose.mhd', file=doseActor_file)
            print('/gate/actor/dose/save                   {:s}'.format(tempPath), file=doseActor_file)
            print('', file=doseActor_file)
        if "doseP" in actors:
            print('# doseP', file=doseActor_file)
            print('/gate/actor/addActor DoseActor               doseP', file=doseActor_file)
            print('/gate/actor/doseP/attachTo                   patientCT', file=doseActor_file)
            print('/gate/actor/doseP/stepHitType                random', file=doseActor_file)
            print('/gate/actor/doseP/setResolution              {:d} {:d} {:d} # put here the pixel resolution from CT mhd.'.format(MHDdims[0],MHDdims[1],MHDdims[2]), file=doseActor_file)            
            print('/gate/actor/doseP/addFilter                  particleFilter', file=doseActor_file)
            print('/gate/actor/doseP/particleFilter/addParticle proton', file=doseActor_file)
            print('/gate/actor/doseP/enableEdep                 false', file=doseActor_file)
            print('/gate/actor/doseP/enableUncertaintyEdep      false', file=doseActor_file)
            print('/gate/actor/doseP/enableSquaredEdep          false', file=doseActor_file)
            print('/gate/actor/doseP/enableDose                 true', file=doseActor_file)
            print('/gate/actor/doseP/enableUncertaintyDose      false', file=doseActor_file)
            print('/gate/actor/doseP/enableSquaredDose          false', file=doseActor_file)
            print('/gate/actor/doseP/normaliseDoseToMax         false', file=doseActor_file)
            print('/gate/actor/doseP/normaliseDoseToIntegral    false', file=doseActor_file)
            print('/gate/actor/doseP/saveEveryNSeconds          600', file=doseActor_file)
            tempPath = os.path.join("Results", "doseP.mhd")
            print('/gate/actor/doseP/save                       {:s}'.format(tempPath), file=doseActor_file)
            print('', file=doseActor_file)
        if "LETdAll" in actors:
            print('# LETdAll', file=doseActor_file)
            print('/gate/actor/addActor    LETActor               LETdAll', file=doseActor_file)
            print('/gate/actor/LETdAll/attachTo                   patientCT', file=doseActor_file)
            print('/gate/actor/LETdAll/setResolution              {:d} {:d} {:d} # put here the pixel resolution from CT mhd.'.format(MHDdims[0],MHDdims[1],MHDdims[2]), file=doseActor_file)            
            print('/gate/actor/LETdAll/setType                    DoseAveraged', file=doseActor_file)
            print('/gate/actor/LETdAll/doParallelCalculation      true', file=doseActor_file)
            print('/gate/actor/LETdAll/saveEveryNSeconds          600', file=doseActor_file)
            tempPath = os.path.join("Results", "LETdAll.mhd")
            print('/gate/actor/LETdAll/save                       {:s}'.format(tempPath), file=doseActor_file)
            print('', file=doseActor_file)
        if "LETdP" in actors:
            print('# LETdP', file=doseActor_file)
            print('/gate/actor/addActor    LETActor             LETdP', file=doseActor_file)
            print('/gate/actor/LETdP/attachTo                   patientCT', file=doseActor_file)
            print('/gate/actor/LETdP/setResolution              {:d} {:d} {:d} # put here the pixel resolution from CT mhd.'.format(MHDdims[0],MHDdims[1],MHDdims[2]), file=doseActor_file)            
            print('/gate/actor/LETdP/setType                    DoseAveraged', file=doseActor_file)
            print('/gate/actor/LETdP/addFilter                  particleFilter', file=doseActor_file)
            print('/gate/actor/LETdP/particleFilter/addParticle proton', file=doseActor_file)
            print('/gate/actor/LETdP/doParallelCalculation      true', file=doseActor_file)
            print('/gate/actor/LETdP/saveEveryNSeconds          600', file=doseActor_file)
            tempPath = os.path.join("Results", "LETdP.mhd")
            print('/gate/actor/LETdAll/save                       {:s}'.format(tempPath), file=doseActor_file)
            print('', file=doseActor_file)
        if "prodAll" in actors: 
            # O15
            print('# prodAll_O15', file=doseActor_file)
            print('/gate/actor/addActor ProductionAndStoppingActor  betaO15', file=doseActor_file)
            tempPath = os.path.join("Results", "O15_map.txt")
            print('/gate/actor/betaO15/save                         {:s}'.format(tempPath), file=doseActor_file)
            print('/gate/actor/betaO15/attachTo                     patientCT', file=doseActor_file)
            print('/gate/actor/betaO15/stepHitType                  post', file=doseActor_file)
            print('/gate/actor/betaO15/setVoxelSize                 {:.1f} {:.1f} {:.1f} # put here the pixel spacing from CT mhd.'.format(MHDspacing[0],MHDspacing[1],MHDspacing[2]), file=doseActor_file)
            print('/gate/actor/betaO15/addFilter                    particleFilter', file=doseActor_file)
            print('/gate/actor/betaO15/particleFilter/addParticle   O15', file=doseActor_file)
            print('', file=doseActor_file)
            # O14
            print('# prodAll_O14', file=doseActor_file)
            print('/gate/actor/addActor ProductionAndStoppingActor  betaO14', file=doseActor_file)
            tempPath = os.path.join("Results", "O14_map.txt")
            print('/gate/actor/betaO14/save                         {:s}'.format(tempPath), file=doseActor_file)
            print('/gate/actor/betaO14/attachTo                     patientCT', file=doseActor_file)
            print('/gate/actor/betaO14/stepHitType                  post', file=doseActor_file)
            print('/gate/actor/betaO14/setVoxelSize                 {:.1f} {:.1f} {:.1f} # put here the pixel spacing from CT mhd.'.format(MHDspacing[0],MHDspacing[1],MHDspacing[2]), file=doseActor_file)
            print('/gate/actor/betaO14/addFilter                    particleFilter', file=doseActor_file)
            print('/gate/actor/betaO14/particleFilter/addParticle   O14', file=doseActor_file)
            print('', file=doseActor_file)
            # N13
            print('# prodAll_N13', file=doseActor_file)
            print('/gate/actor/addActor ProductionAndStoppingActor  betaN13', file=doseActor_file)
            tempPath = os.path.join("Results", "N13_map.txt")
            print('/gate/actor/betaN13/save                         {:s}'.format(tempPath), file=doseActor_file)
            print('/gate/actor/betaN13/attachTo                     patientCT', file=doseActor_file)
            print('/gate/actor/betaN13/stepHitType                  post', file=doseActor_file)
            print('/gate/actor/betaN13/setVoxelSize                 {:.1f} {:.1f} {:.1f} # put here the pixel spacing from CT mhd.'.format(MHDspacing[0],MHDspacing[1],MHDspacing[2]), file=doseActor_file)
            print('/gate/actor/betaN13/addFilter                    particleFilter', file=doseActor_file)
            print('/gate/actor/betaN13/particleFilter/addParticle   N13', file=doseActor_file)
            print('', file=doseActor_file)
            # C11
            print('# prodAll_C11', file=doseActor_file)
            print('/gate/actor/addActor ProductionAndStoppingActor  betaC11', file=doseActor_file)
            tempPath = os.path.join("Results", "C11_map.txt")
            print('/gate/actor/betaC11/save                         {:s}'.format(tempPath), file=doseActor_file)
            print('/gate/actor/betaC11/attachTo                     patientCT', file=doseActor_file)
            print('/gate/actor/betaC11/stepHitType                  post', file=doseActor_file)
            print('/gate/actor/betaC11/setVoxelSize                 {:.1f} {:.1f} {:.1f} # put here the pixel spacing from CT mhd.'.format(MHDspacing[0],MHDspacing[1],MHDspacing[2]), file=doseActor_file)
            print('/gate/actor/betaC11/addFilter                    particleFilter', file=doseActor_file)
            print('/gate/actor/betaC11/particleFilter/addParticle   C11', file=doseActor_file)
            print('', file=doseActor_file)
            # C10
            print('# prodAll_C10', file=doseActor_file)
            print('/gate/actor/addActor ProductionAndStoppingActor  betaC10', file=doseActor_file)
            tempPath = os.path.join("Results", "C10_map.txt")
            print('/gate/actor/betaC10/save                         {:s}'.format(tempPath), file=doseActor_file)
            print('/gate/actor/betaC10/attachTo                     patientCT', file=doseActor_file)
            print('/gate/actor/betaC10/stepHitType                  post', file=doseActor_file)
            print('/gate/actor/betaC10/setVoxelSize                 {:.1f} {:.1f} {:.1f} # put here the pixel spacing from CT mhd.'.format(MHDspacing[0],MHDspacing[1],MHDspacing[2]), file=doseActor_file)
            print('/gate/actor/betaC10/addFilter                    particleFilter', file=doseActor_file)
            print('/gate/actor/betaC10/particleFilter/addParticle   C10', file=doseActor_file)
            print('', file=doseActor_file)
            # P30
            print('# prodAll_P30', file=doseActor_file)
            print('/gate/actor/addActor ProductionAndStoppingActor  betaP30', file=doseActor_file)
            tempPath = os.path.join("Results", "P30_map.txt")
            print('/gate/actor/betaP30/save                         {:s}'.format(tempPath), file=doseActor_file)
            print('/gate/actor/betaP30/attachTo                     patientCT', file=doseActor_file)
            print('/gate/actor/betaP30/stepHitType                  post', file=doseActor_file)
            print('/gate/actor/betaP30/setVoxelSize                 {:.1f} {:.1f} {:.1f} # put here the pixel spacing from CT mhd.'.format(MHDspacing[0],MHDspacing[1],MHDspacing[2]), file=doseActor_file)
            print('/gate/actor/betaP30/addFilter                    particleFilter', file=doseActor_file)
            print('/gate/actor/betaP30/particleFilter/addParticle   P30', file=doseActor_file)
            print('', file=doseActor_file)
            # K38
            print('# prodAll_K38', file=doseActor_file)
            print('/gate/actor/addActor ProductionAndStoppingActor  betaK38', file=doseActor_file)
            tempPath = os.path.join("Results", "K38_map.txt")
            print('/gate/actor/betaK38/save                         {:s}'.format(tempPath), file=doseActor_file)
            print('/gate/actor/betaK38/attachTo                     patientCT', file=doseActor_file)
            print('/gate/actor/betaK38/stepHitType                  post', file=doseActor_file)
            print('/gate/actor/betaK38/setVoxelSize                 {:.1f} {:.1f} {:.1f} # put here the pixel spacing from CT mhd.'.format(MHDspacing[0],MHDspacing[1],MHDspacing[2]), file=doseActor_file)
            print('/gate/actor/betaK38/addFilter                    particleFilter', file=doseActor_file)
            print('/gate/actor/betaK38/particleFilter/addParticle   K38', file=doseActor_file)
            print('', file=doseActor_file)
        doseActor_file.close()

        # patientCT.mac
        patientCT_fileName = os.path.abspath(os.path.join(localDir, "patientCT.mac"))
        patientCT_file=open(patientCT_fileName,'w')
        print('# CT calibration and patient CT', file=patientCT_file)
        print('/gate/world/daughters/name            patientCT', file=patientCT_file)
        print('/gate/world/daughters/insert          ImageNestedParametrisedVolume', file=patientCT_file)
        
        tempPath = os.path.abspath(os.path.join(localDir ,"HUmaterials.db"))
        print('/gate/geometry/setMaterialDatabase            {:s}'.format(tempPath), file=patientCT_file)
        tempPath = os.path.abspath(os.path.join(localDir ,"HU2mat.mac"))
        print('/gate/patientCT/geometry/setHUToMaterialFile  {:s}'.format(tempPath), file=patientCT_file)
        
        tempPath = os.path.abspath(os.path.join(localDir ,"MATCT_Ipot.mac"))
        print('/control/execute         {:s}'.format(tempPath), file=patientCT_file)

        # print('/gate/patientCT/geometry/setImage             ./../../TPS/CT_ExtROICrop_2.0x2.0x2.0.mhd', file=patientCT_file)
        print('/gate/patientCT/geometry/setImage             {:s}'.format(CTpatientFileName), file=patientCT_file)
        print('', file=patientCT_file)
        print('# rotate and translate patientCT', file=patientCT_file)

        print('/gate/patientCT/geometry/TranslateTheImageAtThisIsoCenter  {:+.2f}\t{:+.2f}\t{:+.2f}\t mm # put here the isocentre from RN dicom plan'.format(fieldsInfo[fieldIdx]['IsocenterPosition_mm'][0]-robustXYZ[0],fieldsInfo[fieldIdx]['IsocenterPosition_mm'][1]-robustXYZ[1],fieldsInfo[fieldIdx]['IsocenterPosition_mm'][2]-robustXYZ[2]), file=patientCT_file)        
        print('/gate/patientCT/geometry/setRotationAroundPixelIsoCenter true', file=patientCT_file)
        print('/gate/patientCT/placement/setRotationAxis           0 1 0', file=patientCT_file)
        rotationAngle = fieldsInfo[fieldIdx]['CouchAngle_deg'] * -1.0  #see comment below (MINUS sign required)
        # print('/gate/patientCT/placement/setRotationAngle          {:+.1f} deg # put here a minus table rotation (negative value) from RN dicom plan'.format(fieldsInfo[fieldIdx]['CouchAngle_deg']), file=patientCT_file)
        print('/gate/patientCT/placement/setRotationAngle          {:+.1f} deg # put here a minus table rotation (negative value) from RN dicom plan'.format(rotationAngle), file=patientCT_file)        
        patientCT_file.close()

        
        # copy materials & physics
        if robustXYZ[3]== 0:
            copyfile(os.path.abspath(os.path.join(confPath, "Materials", "HU2Mat_CCB_CTCTA_simplified_Air0.mac")), os.path.abspath(os.path.join(localDir, "HU2mat.mac")))
            copyfile(os.path.abspath(os.path.join(confPath, "Materials", "HUmaterials_CCB_CTCTA_simplified_Air0HU.db")), os.path.abspath(os.path.join(localDir ,"HUmaterials.db")))
            copyfile(os.path.abspath(os.path.join(confPath, "Materials", "MATCT_Ipot.mac")), os.path.abspath(os.path.join(localDir ,"MATCT_Ipot.mac")))
        elif robustXYZ[3]== -2:
            copyfile(os.path.abspath(os.path.join(confPath, "Materials", "HU2Mat_CCB_CTCTA_simplified_Air0_-2.0pr.mac")), os.path.abspath(os.path.join(localDir, "HU2mat.mac")))
            copyfile(os.path.abspath(os.path.join(confPath, "Materials", "HUmaterials_CCB_CTCTA_simplified_Air0HU_-2.0pr.db")), os.path.abspath(os.path.join(localDir ,"HUmaterials.db")))
            copyfile(os.path.abspath(os.path.join(confPath, "Materials", "MATCT_Ipot.mac")), os.path.abspath(os.path.join(localDir ,"MATCT_Ipot.mac")))
        elif robustXYZ[3]== 2:
            copyfile(os.path.abspath(os.path.join(confPath, "Materials", "HU2Mat_CCB_CTCTA_simplified_Air0_+2.0pr.mac")), os.path.abspath(os.path.join(localDir, "HU2mat.mac")))
            copyfile(os.path.abspath(os.path.join(confPath, "Materials", "HUmaterials_CCB_CTCTA_simplified_Air0HU_+2.0pr.db")), os.path.abspath(os.path.join(localDir ,"HUmaterials.db")))
            copyfile(os.path.abspath(os.path.join(confPath, "Materials", "MATCT_Ipot.mac")), os.path.abspath(os.path.join(localDir ,"MATCT_Ipot.mac")))


        #check if there is 'RS_Block' name in RangeShifterID field
        if len(fieldsInfo[fieldIdx]['RangeShifterID']) > 0 : 
            #if RS is in the field - copy appropriate physics into physics.mac
            copyfile(os.path.abspath(os.path.join(confPath, "Physics", "protonRTRS_patient.mac")), os.path.abspath(os.path.join(localDir, "physics.mac")) )
            # copyfile(os.path.abspath(os.path.join(confPath, "Physics", "protonRTRS_Phantom_step=0.1mm.mac")), os.path.abspath(os.path.join(localDir, "physics.mac")) )        
        else :
            #if NO RS in the field - copy appropriate physics into physics.mac
            copyfile(os.path.abspath(os.path.join(confPath, "Physics", "protonRT_patient.mac")), os.path.abspath(os.path.join(localDir, "physics.mac")) )
            # copyfile(os.path.abspath(os.path.join(confPath, "Physics", "protonRT_Phantom_step=0.1mm.mac")), os.path.abspath(os.path.join(localDir, "physics.mac")) )
        


        #mainSim.mac
        mainSim_fileName = os.path.abspath(os.path.join(localDir, "mainSim.mac"))
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
        # print('/control/execute {:s}Visualisation/visu_Disable.mac'.format(confPath), file=mainSim_file)
        print('/control/execute {:s}'.format(tempPath), file=mainSim_file)
        # print('/control/execute {:s}Visualisation/visu_Enable.mac'.format(confPath), file=mainSim_file)
        print('', file=mainSim_file)   # empty line 
        print('#=====================================================', file=mainSim_file)
        print('# VERBOSITY', file=mainSim_file)
        print('#=====================================================', file=mainSim_file)
        print('', file=mainSim_file)   # empty line 
        tempPath = os.path.abspath(os.path.join(confPath, "Verbose", "verbose.mac"))
        # print('/control/execute {:s}/Verbose/verbose.mac'.format(confPath), file=mainSim_file)
        print('/control/execute {:s}'.format(tempPath), file=mainSim_file)
        print('', file=mainSim_file)   # empty line 
        print('#=====================================================', file=mainSim_file)
        print('# MATERIALS', file=mainSim_file)
        print('#=====================================================', file=mainSim_file)
        print('', file=mainSim_file)   # empty line 
        tempPath = os.path.abspath(os.path.join(confPath, "Materials", "GateMaterials.db"))
        # print('/gate/geometry/setMaterialDatabase {:s}Materials/GateMaterials.db'.format(confPath), file=mainSim_file)
        print('/gate/geometry/setMaterialDatabase {:s}'.format(tempPath), file=mainSim_file)
        print('', file=mainSim_file)   # empty line 
        # print('# Generate materials from Hounsfield units', file=mainSim_file)
        # print('/gate/HounsfieldMaterialGenerator/SetMaterialTable                  {:s}Materials/HU2Mat_CCB_CTCTA_simplified_Air0HU_MaterialsTable.mac'.format(confPath), file=mainSim_file)
        # print('/gate/HounsfieldMaterialGenerator/SetDensityTable                   {:s}Materials/HU2Mat_CCB_CTCTA_simplified_Air0HU_DensitiesTable.mac'.format(confPath), file=mainSim_file)
        # print('/gate/HounsfieldMaterialGenerator/SetDensityTolerance               0.05 g/cm3', file=mainSim_file)
        # print('/gate/HounsfieldMaterialGenerator/SetOutputMaterialDatabaseFilename HUmaterials.db', file=mainSim_file)
        # print('/gate/HounsfieldMaterialGenerator/SetOutputHUMaterialFilename       HU2mat.mac', file=mainSim_file)
        # print('/gate/HounsfieldMaterialGenerator/Generate', file=mainSim_file)
        print('', file=mainSim_file)   # empty line 
        print('#=====================================================', file=mainSim_file)
        print('# GEOMETRY', file=mainSim_file)
        print('#=====================================================', file=mainSim_file)
        print('', file=mainSim_file)   # empty line 
        tempPath = os.path.abspath(os.path.join(confPath, "Geometry", "emptyWorld.mac"))
        # print('/control/execute {:s}/Geometry/emptyWorld.mac'.format(confPath), file=mainSim_file)
        print('/control/execute {:s}'.format(tempPath), file=mainSim_file)
        print('', file=mainSim_file)   # empty line 
        tempPath = os.path.abspath(os.path.join(localDir ,"RS.mac"))
        print('/control/execute {:s}'.format(tempPath), file=mainSim_file)
        tempPath = os.path.abspath(os.path.join(localDir ,"patientCT.mac"))
        print('/control/execute {:s}'.format(tempPath), file=mainSim_file)
        print('', file=mainSim_file)   # empty line 
        print('#=====================================================', file=mainSim_file)
        print('# PHYSICS', file=mainSim_file)
        print('#=====================================================', file=mainSim_file)
        print('', file=mainSim_file)   # empty line 
        # print('# /control/execute {:s}/Physics/protonRTRS_Phantom.mac'.format(confPath), file=mainSim_file)
        tempPath = os.path.abspath(os.path.join(localDir ,"physics.mac"))
        print('/control/execute {:s}'.format(tempPath), file=mainSim_file)
        print('', file=mainSim_file)   # empty line 
        print('#=====================================================', file=mainSim_file)
        print('# ACTORS', file=mainSim_file)
        print('#=====================================================', file=mainSim_file)
        print('', file=mainSim_file)   # empty line 
        print('# Statistics for world', file=mainSim_file)
        tempPath = os.path.abspath(os.path.join( confPath, "Actors", "World_statistics.mac"))
        # print('/control/execute {:s}Actors/World_statistics.mac'.format(confPath), file=mainSim_file)
        print('/control/execute {:s}'.format(tempPath), file=mainSim_file)
        print('', file=mainSim_file)   # empty line 
        print('# Dose and production actors for the patientCT', file=mainSim_file)
        tempPath = os.path.abspath(os.path.join(localDir ,"Actors.mac"))
        print('/control/execute {:s}'.format(tempPath), file=mainSim_file)
        print('', file=mainSim_file)   # empty line 
        print('#=====================================================', file=mainSim_file)
        print('# INITIALISATION', file=mainSim_file)
        print('#=====================================================', file=mainSim_file)
        print('', file=mainSim_file)   # empty line 
        print('/gate/run/initialize', file=mainSim_file)
        # print('#/gate/physics/print Results/physics.txt', file=mainSim_file)
        # tempPath = os.path.join( confPath, "Results", "physics.txt")
        # print('#/gate/physics/print {:s}'.format(tempPath), file=mainSim_file)
        print('', file=mainSim_file)   # empty line 
        print('#=====================================================', file=mainSim_file)
        print('# SOURCE', file=mainSim_file)
        print('#=====================================================', file=mainSim_file)
        print('', file=mainSim_file)   # empty line 
        tempPath = os.path.abspath(os.path.join(localDir ,"source.mac"))
        print('/control/execute {:s}'.format(tempPath), file=mainSim_file)
        print('', file=mainSim_file)   # empty line 
        print('#=====================================================', file=mainSim_file)
        print('# RANDOM', file=mainSim_file)
        print('#=====================================================', file=mainSim_file)
        print('', file=mainSim_file)   # empty line 
        print('/gate/random/setEngineName MersenneTwister', file=mainSim_file)
        print('/gate/random/setEngineSeed auto', file=mainSim_file)
        print('/gate/random/verbose 0', file=mainSim_file)
        print('', file=mainSim_file)   # empty line 
        print('#=====================================================', file=mainSim_file)
        print('# ACQUISITION SETTINGS', file=mainSim_file)
        print('#=====================================================', file=mainSim_file)
        print('', file=mainSim_file)   # empty line 
        # print('/gate/application/setTotalNumberOfPrimaries 1E5', file=mainSim_file)
        NbPrimaries = NBparts[fieldIdx]  # taken from fieldStat.csv FILE
        print('/gate/application/setTotalNumberOfPrimaries {:d}'.format(NbPrimaries), file=mainSim_file)
        print('/gate/application/setTimeStart 0. s', file=mainSim_file)
        print('/gate/application/setTimeStop 100. s', file=mainSim_file)
        print('/gate/application/startDAQ', file=mainSim_file)
        print('', file=mainSim_file)   # empty line 
        print('#=====================================================', file=mainSim_file)
        print('# POSTPROCESSING', file=mainSim_file)
        print('#=====================================================', file=mainSim_file)
        # print('# /control/shell echo "\e[36m$(cat Results/simStat.out | grep '.*\(PPS\|ElapsedTime\)')\e[0m"', file=mainSim_file)
        mainSim_file.close()




def writeRtplan(rtplanfileName, planInfo, fieldsInfo, getPlanVersion='', fieldToSave=None, displayInfo=True):
    import sys, os, glob
    import pandas as pd
    import pydicom as dicom
    import numpy as np
    from shutil import copyfile
    Xsec='emittance'

    rtplanfileName=os.path.abspath(rtplanfileName)
    
    # remove old beam model file, RS model file and material definition file from simulation folder
    if glob.glob(os.path.join(os.path.split(os.path.abspath(rtplanfileName))[0],'*beamModel.txt')):
        os.remove(glob.glob(os.path.join(os.path.split(os.path.abspath(rtplanfileName))[0],'*beamModel.txt'))[0])
    if glob.glob(os.path.join(os.path.split(os.path.abspath(rtplanfileName))[0],'*RSModel.txt')):
        os.remove(glob.glob(os.path.join(os.path.split(os.path.abspath(rtplanfileName))[0],'*RSModel.txt'))[0])
    if glob.glob(os.path.join(os.path.split(os.path.abspath(rtplanfileName))[0],'*materialDefinition.txt')):
        os.remove(glob.glob(os.path.join(os.path.split(os.path.abspath(rtplanfileName))[0],'*materialDefinition.txt'))[0])

    # copy current beam model file, RS model file and material definition file to simulation folder
    copyfile(planInfo['beamModelFileName'], os.path.join(os.path.split(os.path.abspath(rtplanfileName))[0],os.path.basename(planInfo['beamModelFileName'])))
    copyfile(planInfo['RSModelFileName'], os.path.join(os.path.split(os.path.abspath(rtplanfileName))[0],os.path.basename(planInfo['RSModelFileName'])))
    copyfile(planInfo['RSModelFileName'], os.path.join(os.path.split(os.path.abspath(rtplanfileName))[0],'regions.inp'))
    copyfile(planInfo['materialDefinitionFileName'], os.path.join(os.path.split(os.path.abspath(rtplanfileName))[0],os.path.basename(planInfo['materialDefinitionFileName'])))
    copyfile(planInfo['materialDefinitionFileName'], os.path.join(os.path.split(os.path.abspath(rtplanfileName))[0],'materials.inp'))
    rtplan_H=open(rtplanfileName,'w')
    
    # choose field to save or save all fields
    if not fieldToSave==None:
        if   isinstance(fieldToSave, int):
            if fieldToSave in planInfo['fieldIDList']:
                fieldToSaveIdx=planInfo['fieldIDList'].index(fieldToSave)
                print('# Requested to save the rtplan only for the field ID {:d}'.format(fieldToSave))
            else:
                print('# Error: cannot find field ID {:d} in the treatment plan. Available field IDs are '.format(fieldToSave)+str(planInfo['fieldIDList']))
                exit(-1)
        elif isinstance(fieldToSave, str):
            if fieldToSave in planInfo['fieldNameList']:    
                fieldToSaveIdx=planInfo['fieldNameList'].index(fieldToSave)
                
                print('# Requested to save the rtplan only for the field name \'{:s}\''.format(fieldToSave))
            else:
                print('# Error: cannot find field name \'{:s}\' in the treatment plan. Available field names are '.format(fieldToSave)+str(planInfo['fieldNameList']))
                exit(-1)
        fieldsInfo=[fieldsInfo[fieldToSaveIdx]]   
        planInfo['pencilBeamNumberAll']=fieldsInfo[0]['PencilBeams']
        planInfo['protonNumberAll']=fieldsInfo[0]['PencilBeamsProtonNumber']
        
    # write information header
    print('###############################################################################################################################################', file=rtplan_H)
    print('###### Beam plan for {:s} machine prepared in getPlan v.{:s}'.format(planInfo['TreatmentMachineName'], getPlanVersion), file=rtplan_H)
    if isinstance(planInfo['RTplanFileName'],list):
        for RTplanFileName_idx, RTplanFileName in enumerate(planInfo['RTplanFileName']):
            print('###### DICOM plan {:d}: {:s}'.format(RTplanFileName_idx+1,os.path.abspath(RTplanFileName)), file=rtplan_H)
    else:
        print('###### DICOM plan: {:s}'.format(os.path.abspath(planInfo['RTplanFileName'])), file=rtplan_H)
    print('###### PB no.: {:d}, proton no.: {:.20E}'.format(planInfo['pencilBeamNumberAll'],planInfo['protonNumberAll']), file=rtplan_H)
    print('###############################################################################################################################################', file=rtplan_H)    
    
    # write field for G0
    for fieldInfo in fieldsInfo:
        print('field: {:d} ; '.format(fieldInfo['BeamNumber']), file=rtplan_H, end='')
        print('O=[ {:+.4f}, {:+.4f}, {:+.4f} ];\t'.format(0,-fieldInfo['fieldOriginToIsoDist_mm']/10,0), file=rtplan_H, end='')
        print('f=[ {:+.3f}, {:+.3f}, {:+.3f} ];\t'.format(0, +1, 0), file=rtplan_H, end='')
        print('u=[ {:+.3f}, {:+.3f}, {:+.3f} ];\t'.format(0, 0, +1), file=rtplan_H, end='')
        print('exitWindowPlane={:.4f};'.format((fieldInfo['fieldOriginToIsoDist_mm']-fieldInfo['VaccumExitToISODist_mm'])/10), file=rtplan_H)
    
    #write master for field
    for fieldInfo in fieldsInfo:
        print('pbmaster: {:d} ; '.format(fieldInfo['BeamNumber']), file=rtplan_H, end='')
        print('particle={:s};\t'.format(fieldInfo['RadiationType'].lower()), file=rtplan_H, end='')
        print('Xsec={:s};\t'.format(Xsec), file=rtplan_H, end='')
        print('emittanceRefPlaneDistance={:.3f};\t'.format(fieldInfo['fieldOriginToIsoDist_mm']/10), file=rtplan_H, end='')
        print('columns=[P.x,P.y,P.z,v.x,v.y,v.z,T,pSpread,N,twissAlphaX,twissBetaX,emittanceX,twissAlphaY,twissBetaY,emittanceY] ', file=rtplan_H)

    # group fields, prepare gantry region, save the default regions FoR
    print('group: fields {:s}'.format(' '.join('field_{:d}'.format(fieldInfo['BeamNumber'])  for fieldInfo in fieldsInfo)), file=rtplan_H)
    print('region: gantry ; O = [0,0,0] ; L = [1,1,1] ; lAdaptiveSize=t ; material = vacuum', file=rtplan_H) 
    print('set_parent: gantry fields NozzleGroup', file=rtplan_H)    
    print('save_regions: 0', file=rtplan_H)
    
    print('##################################################', file=rtplan_H)
    print('###### Start of setup and delivery sequence ######', file=rtplan_H)
    print('##################################################', file=rtplan_H)    
    for fieldInfo in fieldsInfo:
        print('###### setup sequence for field {:d}'.format(fieldInfo['BeamNumber']), file=rtplan_H)
        
        print('# translate nozzle regions by snaut position', file=rtplan_H)
        print('transform: NozzleGroup shift_by {:.4f} {:.4f} {:.4f}'.format(0,-fieldInfo['SnoutPosition']/10,0), file=rtplan_H)    
        
        print('# rotate field and NozzleGroup by gantry angle', file=rtplan_H)
        print('transform: gantry rotate {:s} {:.4f}'.format('z',fieldInfo['GantryAngle_deg']), file=rtplan_H)
        
        print('# translate and rotate phantom by iso position and couch rotation', file=rtplan_H)
        print('transform: phantom shift_by {:+.4f} {:+.4f} {:+.4f}'.format(-fieldInfo['IsocenterPosition_mm'][0]/10,-fieldInfo['IsocenterPosition_mm'][1]/10,-fieldInfo['IsocenterPosition_mm'][2]/10), file=rtplan_H)
        print('transform: phantom rotate {:s} {:+.4f}'.format('y',-fieldInfo['CouchAngle_deg']), file=rtplan_H)
        
        print('# delivery sequences for slices in field {:d}'.format(fieldInfo['BeamNumber']), file=rtplan_H)
        print('deactivate: fields, NozzleGroup', file=rtplan_H)
        print('activate: field_{:d} {:s}'.format(fieldInfo['BeamNumber'],fieldInfo['RangeShifterID']), file=rtplan_H)
        print('deliver: field_{:d}'.format(fieldInfo['BeamNumber']), file=rtplan_H)
        
        print('# reset field, nozzle regions and phantom to default FoR', file=rtplan_H)
        print('restore: 0', file=rtplan_H)

    print('##################################################', file=rtplan_H)
    print('###### End of setup and delivery sequence ########', file=rtplan_H)
    print('##################################################', file=rtplan_H)


    spotNumber=0
    for fieldInfo in fieldsInfo:
        for sliceInfo in fieldInfo['slicesInfo']:
            for _,sliceSpot in sliceInfo['sliceSpots'].iterrows():
                spotNumber+=1
                print('pb: {:d}\t{:d}\t'.format(spotNumber,fieldInfo['BeamNumber']), file=rtplan_H, end='')
                print('{:+.10f}\t{:+.10f}\t{:+.10f}\t'.format(sliceSpot['FRED_Px'],sliceSpot['FRED_Py'],sliceSpot['FRED_Pz']), file=rtplan_H, end='')
                print('{:+.10f}\t{:+.10f}\t{:+.10f}\t'.format(sliceSpot['FRED_Vx'],sliceSpot['FRED_Vy'],sliceSpot['FRED_Vz']), file=rtplan_H, end='')
                print('{:.3f}\t'.format(sliceInfo['FRED_beamEnergy']), file=rtplan_H, end='')
                print('{:.5E}\t'.format(sliceInfo['FRED_pSpread']), file=rtplan_H, end='')
                print('{:.10E}\t'.format(sliceSpot['FRED_protonNumber']), file=rtplan_H, end='')
                print('{:+.10E}\t{:+.10E}\t{:+.10E}\t'.format(sliceInfo['FRED_alphaX'],sliceInfo['FRED_betaX_cm'],sliceInfo['FRED_epsilonX_cm']), file=rtplan_H, end='')
                print('{:+.10E}\t{:+.10E}\t{:+.10E}\t'.format(sliceInfo['FRED_alphaY'],sliceInfo['FRED_betaY_cm'],sliceInfo['FRED_epsilonY_cm']), file=rtplan_H)
    if displayInfo:
        print('# RT plan saved to {:s}'.format(rtplanfileName))
    rtplan_H.close()

def getPlan(PatientPath, CTpatientFileName, simLabel, GATEconfigPath, cores, actors, Nparticles, ParticleUnit, robustXYZ, GATErtplanFileName='current', fieldToSave=None, defaultModelsFileName='default', beamModelInterpolationMethod='slinear', displayInfo=False):
# def getPlan(PatientPath, CTpatientFileName, simLabel, fieldToSave=None, GATEconfigPath='current', cores, actors, Nparticles, ParticleUnit, GATErtplanFileName='current', defaultModelsFileName='default', beamModelInterpolationMethod='slinear', displayInfo=False):
    import os, re    
    # version of getPlan
    version='3.7b'
    if displayInfo:
        print('############################################################')
        print('# Importer of RTPLAN dicom for GATE at CCB v.{:s}'.format(version))
        print('# Implementation works for GATE version 8.2 or higher')
        print('############################################################')
    
    # !!!!!!!!!!!!!!!!!!!! 
    RTplanFileName = os.path.abspath(os.path.join( PatientPath, "TPS", "RN.dcm" ))
    # WE assume that PatientPath contains subdiractory TPS  with RN.dcm file inside

    CTpatientFileName_Path = os.path.abspath(os.path.join( PatientPath, "TPS", CTpatientFileName ))

    # print(GATEconfigPath)

    # find rtplan file name
    if GATErtplanFileName=='current':
        GATErtplanFileName=os.path.join(os.getcwd(),'rtplan.mac')

    # get gantry name function
    def getGantryName(RTplanFileName):
        dicomRT=loadDicomRT(RTplanFileName, displayInfo=False)
        for ifield in range(dicomRT.FractionGroupSequence[0].NumberOfBeams):
            if dicomRT.IonBeamSequence[ifield].TreatmentDeliveryType=='TREATMENT':
                TreatmentMachineName=dicomRT.IonBeamSequence[ifield].TreatmentMachineName
        if not TreatmentMachineName in ['GTR3','GTR4']:
            print('Room name \'{:s}\' in {:s} file not recognised'.format(TreatmentMachineName, RTplanFileName))
            exit(1)
        return TreatmentMachineName
    
    # get beam model, RS model and material definition file names function
    def getBeamModelFileNames(defaultModelsFileName,TreatmentMachineName):
        with open(defaultModelsFileName) as f:
            for num, line in enumerate(f, 1):    
                beamModelFileName_re = re.search('{:s}beamModelFileName\W+=\W+\'(\S+)\''.format(TreatmentMachineName), line)
                if beamModelFileName_re:
                    beamModelFileName=beamModelFileName_re.group(1)
                RSModelFileName_re = re.search('{:s}RSModelFileName\W+=\W+\'(\S+)\''.format(TreatmentMachineName), line)
                if RSModelFileName_re:
                    RSModelFileName=RSModelFileName_re.group(1)
                materialDefinitionFileName_re = re.search('{:s}materialDefinitionFileName\W+=\W+\'(\S+)\''.format(TreatmentMachineName), line)
                if materialDefinitionFileName_re:
                    materialDefinitionFileName=materialDefinitionFileName_re.group(1)
        return beamModelFileName, RSModelFileName, materialDefinitionFileName 

    # get default beam model
    if defaultModelsFileName=='default':
        # defaultModelsFileName='/home/shared/fred/beamModel/fred/defaultBeamModels/CCB_defaultModels.txt'
        defaultModelsFileName='/Users/dborys/Documents/POLSL/2020_IFJ_postdoc/getPlanTests/example1/fred/GTR3_20190919_beamModel.txt'
    else:
        if not os.path.isfile(defaultModelsFileName):
            print('Could not load file with default beam models (file {:s} does not exist)'.format(defaultModelsFileName))
            exit(1)

    # check if RTplanFileName is a single string or list of strings
    if isinstance(RTplanFileName, list):
        ### list of RTplanFileNames ###
        if all(isinstance(item, str) for item in RTplanFileName):
            if displayInfo:
                print('# Requested a list of RT plans. The plan will be saved for all TR plans multiplied by corresponding number of fractions.')

            RTplanFileNames=RTplanFileName
            
            # get TreatmentMachineNames for all plans and check if they are the same
            TreatmentMachineNames=[]
            for RTplanFileName in RTplanFileNames:
                TreatmentMachineNames.append(getGantryName(RTplanFileName))
            if len(set(TreatmentMachineNames))==1:
                TreatmentMachineName=TreatmentMachineNames[0]
            else:
                print('TreatmentMachineName is not the same for all RT plans in list. ({:s})'.format(' '.join(TreatmentMachineNames)))
                exit(1)
                
            # get beam model, RS model and material definition file names
            # beamModelFileName, RSModelFileName, materialDefinitionFileName=getBeamModelFileNames(defaultModelsFileName,TreatmentMachineName)
            
            # build plan for each RT plan in list
            planInfoList=[]
            fieldsInfoList=[]
            for RTplanFileName in RTplanFileNames:
                # get fractionNo
                dicomRT=loadDicomRT(RTplanFileName, displayInfo=False)
                if 'FractionGroupSequence' in dicomRT:
                    if 'NumberOfFractionsPlanned' in dicomRT.FractionGroupSequence[0]:
                        fractionNo=int(dicomRT.FractionGroupSequence[0].NumberOfFractionsPlanned)
                    else:
                        print('No \'NumberOfFractionsPlanned\' in FractionGroupSequence[0] in RT plan {:s}'.format(RTplanFileName))
                        exit(1)
                else:
                    print('No \'FractionGroupSequence\' in RT plan {:s}'.format(RTplanFileName))
                    exit(1)
                
                beamModelFileName = " "
                RSModelFileName = " " 
                materialDefinitionFileName = " "

                planInfo, fieldsInfo=buildPLAN(RTplanFileName, beamModelFileName, RSModelFileName, materialDefinitionFileName, fractionNo = fractionNo, beamModelInterpolationMethod=beamModelInterpolationMethod, displayInfo=displayInfo)
                planInfoList.append(planInfo)
                fieldsInfoList.append(fieldsInfo)
            
            # merge planInfoList to planInfo and fieldsInfoList to fieldsInfo
            planInfo, fieldsInfo=mergePlanInfoAndFieldsInfo(planInfoList, fieldsInfoList)
            
            # write plan
            # writeRtplan(FREDrtpalnFileName, planInfo, fieldsInfo, fieldToSave=None, getPlanVersion=version, displayInfo=displayInfo)

            writeRtplanGate(GATErtplanFileName, CTpatientFileName_Path, simLabel, GATEconfigPath, cores, actors, Nparticles, ParticleUnit, robustXYZ, planInfo, fieldsInfo, fieldToSave=None, getPlanVersion=version, displayInfo=displayInfo)
               
        else:
            print('RTplanFileName is not a list of strings.')
            exit(1)
        
    elif isinstance(RTplanFileName, str):
        ### single RTplanFileName ###
        if displayInfo:
            print('# Requested a single RT plan. The plan will be saved for single fraction.')
    
        # get gantry name function
        TreatmentMachineName=getGantryName(RTplanFileName)
    
        # get beam model, RS model and material definition file names
        # beamModelFileName, RSModelFileName, materialDefinitionFileName=getBeamModelFileNames(defaultModelsFileName,TreatmentMachineName)

        beamModelFileName = " "
        RSModelFileName = " " 
        materialDefinitionFileName = " "

        # get fractionNo
        dicomRT=loadDicomRT(RTplanFileName, displayInfo=False)
        if 'FractionGroupSequence' in dicomRT:
            if 'NumberOfFractionsPlanned' in dicomRT.FractionGroupSequence[0]:
                fractionNo=int(dicomRT.FractionGroupSequence[0].NumberOfFractionsPlanned)
            else:
                print('No \'NumberOfFractionsPlanned\' in FractionGroupSequence[0] in RT plan {:s}'.format(RTplanFileName))
                exit(1)
        else:
            print('No \'FractionGroupSequence\' in RT plan {:s}'.format(RTplanFileName))
            exit(1)


        # build Plan
        planInfo, fieldsInfo=buildPLAN(RTplanFileName, beamModelFileName, RSModelFileName, materialDefinitionFileName, fractionNo = fractionNo, beamModelInterpolationMethod=beamModelInterpolationMethod, displayInfo=displayInfo)
        
        # write plan
        writeRtplanGate(GATErtplanFileName, CTpatientFileName_Path, simLabel, GATEconfigPath, cores, actors, Nparticles, ParticleUnit, robustXYZ, planInfo, fieldsInfo, fieldToSave=fieldToSave, getPlanVersion=version, displayInfo=displayInfo)
    else:
        print('RTplanFileName is not a string.')
        exit(1)
        
