"""
.. py:module:: TuMAGtools.utils
.. module:: utils
        :platform: Unix
        :synopsis: function for reading TuMAG images and headers
.. moduleauthor:: Pablo Santamarina  (SPG - IAA) <psanta@iaa.es>
"""

# ============================ IMPORTS ====================================== #
 
from email.mime import image
import numpy as np
import struct
import matplotlib.pyplot as plt

# ============================= CONFIG ====================================== #

# ------------------------------ SIZE IN BYTES ------------------------------ #

nBytesMpsHeader = 20
nBytesLengthImageStruct = 2
nBytesLengthSSt = 2
nBytesLengthTail = 8  # SSC + CRC

# COMMON sizes
nBytesCameraId        = 1           # cameraId                uint8_t     --> struct.unpack('B') 
nBytesTimeStamp_start = 8           # timeStamp_start         Uint64_t    --> struct.unpack('Q')
nBytesTimeStamp_end   = 8           # timeStamp_end           Uint64_t    --> struct.unpack('Q')
nBytesImageSize       = 4           # imageSize               uint32_t    --> struct.unpack('L')
nBytesObservationMode = 1           # observationMode         uint8_t     --> struct.unpack('B')
nBytesComponentID     = 1           # componentID             uint8_t     --> struct.unpack('B')
nBytesPipelineConfig  = 1           # pipelineConfig          uint8_t     --> struct.unpack('B')
nBytesnAcc            = 2           # NAcc                    uint16_t    --> struct.unpack('H')
nBytesImageIndex      = 4           # imageIndex              uint32_t    --> struct.unpack('L')
nBytesImageIndex_end  = 4           # imageIndex_end          uint32_t    --> struct.unpack('L')
nBytesRoi_x_offset    = 2           # roi_x_offset            uint16_t    --> struct.unpack('H')
nBytesRoi_x_size      = 2           # roi_x_size              uint16_t    --> struct.unpack('H')
nBytesRoi_y_offset    = 2           # roi_y_offset            uint16_t    --> struct.unpack('H')
nBytesRoi_y_size      = 2           # roi_y_size              uint16_t    --> struct.unpack('H')
nBytesImageType       = 1           # Image type              uint8_t     --> struct.unpack('B')
nBytesOMCounter       = 2           # OM_counter              uint16_t    --> struct.unpack('H')

# TuMAG sizes
nBytesFW1               = 1           # position FW1            uint8_t     --> struct.unpack('B')
nBytesFW2               = 1           # position FW2            uint8_t     --> struct.unpack('B')
nBytesEtalonDN          = 2           # voltage etalon DN       uint16_t    --> struct.unpack('H')
nBytesEtalonSign        = 1           # sign etalon             uint8_t     --> struct.unpack('B')
nBytesRocli1_LCVR       = 2           # Rocli1_LCVR             uint16_t    --> struct.unpack('H')
nBytesRocli2_LCVR       = 2           # Rocli2_LCVR             uint16_t    --> struct.unpack('H')
nBytesEtalonVoltLecture = 2           # Etalon Volts reading    uint16_t    --> struct.unpack('H')
nBytesfw1PosReal        = 1           # Real FW1 Pos            uint8_t     --> struct.unpack('B')
nBytesfw2PosReal        = 1           # Real FW2 Pos            uint8_t     --> struct.unpack('B')
nByteslcvr1DNReal       = 2           # Counts lcvr1            uint16_t    --> struct.unpack('H')
nByteslcvr2DNReal       = 2           # Counts lcvr2            uint16_t    --> struct.unpack('H')

# --------------------------------------------------------------------------- #

defaultImgRowSize = 2048
defaultImgColSize = 2048

# --------------------------------------------------------------------------- #

first_type = 10
Types_of_thumbnail = ['0x0A / 10 no data (only header without data)',
         '0x0B / 11 BinnigLos (cam 1/2: 4x4, cam sj FULL)',             
         '0x0C / 12 BinnigNoLos (cam 1/2: 8x8, cam sh 2x2)',           
         '0x0D / 13 Full image (test only)',
         '0x0E / 14 CropCenteredLos',
         '0x0F / 15 CropCenteredNoLos',
         '0x10 / 16 CropUpLeftLos',
         '0x11 / 17 CropUpLeftNoLos',
         '0x12 / 18 CropUpCenterLos',
         '0x13 / 19 CropUpCenterNoLos',
         '0x14 / 20 CropUpRightLos',
         '0x15 / 21 CropUpRightNoLos',
         '0x16 / 22 CropCenterLeftLos',
         '0x17 / 23 CropCenterLeftNoLos',
         '0x18 / 24 CropCenterRightLos',
         '0x19 / 25 CropCenterRightNoLos',
         '0x1A / 26 CropDownLeftLos',
         '0x1B / 27 CropDownLeftNoLos',
         '0x1C / 28 CropDownCenterLos',
         '0x1D / 29 CropDownCenterNoLos',
         '0x1E / 30 CropDownRightLos',
         '0x1F / 31 CropDownRightNoLos']
Types_of_thumbnail_bn = [1,4,8] 

# Si la imagen original es de 4 bytes por pixel, y es de tipo los (line of signal) al thumbnail se calcula aplicando un binning de  8 x 8
# Si la imagen original es de 2 bytes por pixel, y es de tipo los (line of signal) al thumbnail se calcula aplicando un binning de  4 x 4
# Si la imagen original es de 4 bytes por pixel, y es de tipo NO los (no line of signal) al thumbnail se calcula aplicando un binning de  16 x 16
# Si la imagen original es de 2 bytes por pixel, y es de tipo NO los (no line of signal) al thumbnail se calcula aplicando un binning de  8 x 8
# Si es full test la imagen del thumbnail es del size de la original.

# =========================================================================== #


def GetDatafromHeader(receivedHeader):
    """
    This function reads header information
    
    :param receivedHeader: Array of bytes containing the header info
    :type receivedHeader: bytes array

    :return: Header: Dictionary containing the header data
    :rtype: Dict
        
    """
    
    dataLine = []
    
    Keys = [ 'CameraID', 'TimeStamp_start', 'TimeStamp_end', 'ImageSize', 
            'ObservationMode', 'PipelineConfig', 'nAcc', 'ImageIndex', 'ImageIndex_end',
            'Roi_x_offset', 'Roi_x_size', 'Roi_y_offset', 'Roi_y_size', 'ImageType', 'Observation_Counter', 
            'FW1', 'FW2', 'EtalonDN', 'EtalonSign', 'Rocli1_LCVR', 'Rocli2_LCVR', 'EtalonVoltsReading', 
            'FW1_Real', 'FW2_Real', 'LCVR1_DN_Real', 'LCVR2_DN_Real']
    
    # Updating size of Headers
    positionStartSstLength = nBytesMpsHeader
    nBytesSSt = struct.unpack('H', receivedHeader[positionStartSstLength: positionStartSstLength + nBytesLengthSSt])[0]
    positionStartImageStructLength = positionStartSstLength + nBytesLengthSSt + nBytesSSt
    
    # COMMON positions
    positionStartCameraId        = positionStartImageStructLength + nBytesLengthImageStruct 
    positionStartTimeStamp_start = positionStartCameraId + nBytesCameraId
    positionStartTimeStamp_end   = positionStartTimeStamp_start + nBytesTimeStamp_start
    positionStartImageSize       = positionStartTimeStamp_end + nBytesTimeStamp_end
    positionStartObservationMode = positionStartImageSize + nBytesImageSize
    positionStartPipelineConfig  = positionStartObservationMode + nBytesObservationMode
    positionStartnAcc            = positionStartPipelineConfig + nBytesPipelineConfig
    positionStartImageIndex      = positionStartnAcc + nBytesnAcc
    positionStartImageIndex_end  = positionStartImageIndex + nBytesImageIndex
    positionStartRoi_x_offset    = positionStartImageIndex_end + nBytesImageIndex_end
    positionStartRoi_x_size      = positionStartRoi_x_offset + nBytesRoi_x_offset
    positionStartRoi_y_offset    = positionStartRoi_x_size + nBytesRoi_x_size
    positionStartRoi_y_size      = positionStartRoi_y_offset + nBytesRoi_y_offset
    positionStartImageType       = positionStartRoi_y_size + nBytesRoi_y_size
    positionStartOmCounter       = positionStartImageType  + nBytesImageType
    
    # TuMAG positions
    positionFW1                  = positionStartOmCounter + nBytesOMCounter
    positionFW2                  = positionFW1 + nBytesFW1
    positionEtalonDN             = positionFW2 + nBytesFW2
    positionEtalonSign           = positionEtalonDN + nBytesEtalonDN
    positionRocli1_LCVR          = positionEtalonSign + nBytesEtalonSign
    positionRocli2_LCVR          = positionRocli1_LCVR + nBytesRocli1_LCVR
    positionEtalonVoltReading    = positionRocli2_LCVR + nBytesRocli2_LCVR
    positionFW1Real              = positionEtalonVoltReading + nBytesEtalonVoltLecture
    positionFW2Real              = positionFW1Real + nBytesfw1PosReal
    positionLCVR1_DN_real        = positionFW2Real + nBytesfw2PosReal
    positionLCVR2_DN_real        = positionLCVR1_DN_real + nByteslcvr1DNReal
    
    

    # Extract COMMON Data
    dataLine.append(struct.unpack('B', receivedHeader[positionStartCameraId: positionStartCameraId + nBytesCameraId])[0])
    dataLine.append(struct.unpack('Q', receivedHeader[positionStartTimeStamp_start: positionStartTimeStamp_start + nBytesTimeStamp_start])[0])
    dataLine.append(struct.unpack('Q', receivedHeader[positionStartTimeStamp_end: positionStartTimeStamp_end + nBytesTimeStamp_end])[0])
    dataLine.append(struct.unpack('i', receivedHeader[positionStartImageSize: positionStartImageSize + nBytesImageSize])[0])
    dataLine.append(struct.unpack('B', receivedHeader[positionStartObservationMode: positionStartObservationMode + nBytesObservationMode])[0])
    dataLine.append(struct.unpack('B', receivedHeader[positionStartPipelineConfig: positionStartPipelineConfig + nBytesPipelineConfig])[0])
    dataLine.append(struct.unpack('H', receivedHeader[positionStartnAcc: positionStartnAcc + nBytesnAcc])[0])
    dataLine.append(struct.unpack('i', receivedHeader[positionStartImageIndex: positionStartImageIndex + nBytesImageIndex])[0])
    dataLine.append(struct.unpack('i', receivedHeader[positionStartImageIndex_end: positionStartImageIndex_end + nBytesImageIndex_end])[0])
    dataLine.append(struct.unpack('H', receivedHeader[positionStartRoi_x_offset: positionStartRoi_x_offset + nBytesRoi_x_offset])[0])
    dataLine.append(struct.unpack('H', receivedHeader[positionStartRoi_x_size: positionStartRoi_x_size + nBytesRoi_x_size])[0])
    dataLine.append(struct.unpack('H', receivedHeader[positionStartRoi_y_offset: positionStartRoi_y_offset + nBytesRoi_y_offset])[0])
    dataLine.append(struct.unpack('H', receivedHeader[positionStartRoi_y_size: positionStartRoi_y_size + nBytesRoi_y_size])[0])
    dataLine.append(struct.unpack('B', receivedHeader[positionStartImageType: positionStartImageType + nBytesImageType])[0])
    dataLine.append(struct.unpack('H', receivedHeader[positionStartOmCounter: positionStartOmCounter + nBytesOMCounter])[0]) 
    
	# Extract TuMAG Data
    dataLine.append(struct.unpack('B', receivedHeader[positionFW1: positionFW1 + nBytesFW1])[0])
    dataLine.append(struct.unpack('B', receivedHeader[positionFW2: positionFW2 + nBytesFW2])[0])
    dataLine.append(struct.unpack('H', receivedHeader[positionEtalonDN: positionEtalonDN + nBytesEtalonDN])[0])
    dataLine.append(struct.unpack('B', receivedHeader[positionEtalonSign: positionEtalonSign + nBytesEtalonSign])[0])
    dataLine.append(struct.unpack('H', receivedHeader[positionRocli1_LCVR: positionRocli1_LCVR + nBytesRocli1_LCVR])[0])
    dataLine.append(struct.unpack('H', receivedHeader[positionRocli2_LCVR: positionRocli2_LCVR + nBytesRocli2_LCVR])[0])
    dataLine.append(struct.unpack('H', receivedHeader[positionEtalonVoltReading: positionEtalonVoltReading + nBytesEtalonVoltLecture])[0])
    
    # Temporal Parameters
    dataLine.append(struct.unpack('B', receivedHeader[positionFW1Real: positionFW1Real + nBytesfw1PosReal])[0])
    dataLine.append(struct.unpack('B', receivedHeader[positionFW2Real: positionFW2Real + nBytesfw2PosReal])[0])
    dataLine.append(struct.unpack('H', receivedHeader[positionLCVR1_DN_real: positionLCVR1_DN_real + nByteslcvr1DNReal])[0])
    dataLine.append(struct.unpack('H', receivedHeader[positionLCVR2_DN_real: positionLCVR2_DN_real + nByteslcvr2DNReal])[0])
   
    # Saving data in dictionary
    Header = {}
    for ind, key in enumerate(Keys): 
        Header[key] = dataLine[ind]
    
    return Header

def Read_TuMAG(file):
    '''
    Read img images from TuMAG
    returns image and header
    '''

    nBytesSync = 2              # Identify start of record --> 0x44 0x54
    nBytesSystemId = 1          # System ID of generated data
    nBytesDataId = 1            # Description of data type
    nBytesTotalLength = 4       # Total data packet length (including header)
    nBytesTime = 8              # Data acquisition start time
    nBytesSensorId = 1          # 
    nBytesHeaderVersion = 1     # 
    nBytesHeaderLength = 2      # Total length of header (in bytes). Take it as an offset to Data
    
    positionStartTotalLength = nBytesSync + nBytesSystemId + nBytesDataId
    positionStartHeaderLength = nBytesSync + nBytesSystemId + nBytesDataId + nBytesTotalLength + nBytesTime + nBytesSensorId + nBytesHeaderVersion

    data = np.fromfile(file, dtype='u1')
    mpsHeader = data[0 : nBytesMpsHeader]

    totalLength = struct.unpack('i', mpsHeader[positionStartTotalLength: positionStartTotalLength + nBytesTotalLength])[0]
    headerLength = struct.unpack('H', mpsHeader[positionStartHeaderLength: positionStartHeaderLength + nBytesHeaderLength])[0]

    # print(headerLength)

    receivedHeader = data[0:headerLength]
    receivedImage = data[headerLength:totalLength-nBytesLengthTail]

    Head = GetDatafromHeader(receivedHeader) 

    ImageType = Types_of_thumbnail_bn[Head['ImageType'] - first_type]
    print('ImageType: ',Types_of_thumbnail[Head['ImageType'] - first_type],' ',Head['ImageType'],ImageType)

    PipelineConfig = Head['PipelineConfig']
    PipelineConfig = PipelineConfig * 2
    print('PipelineConfig: ',Head['PipelineConfig'],PipelineConfig)

    width, height = Head['Roi_x_size']//ImageType//PipelineConfig, Head['Roi_y_size']//ImageType//PipelineConfig
    dtypef = '<i2'
    if PipelineConfig == 0:
        dtypef = '<i2'
    if PipelineConfig == 1:
        dtypef = '<i4'
    if PipelineConfig == 2:
        dtypef = '<i4'
    # dtypef = '<i2'

    Image = receivedImage.astype(dtypef).reshape([height, width])

    return Image, Head

def HeadernImageSeparator(Image_path):
    
    """
    Function that given image path separates the header from the image in byte 
    arrays

    Parameters
    ----------
    Image_path : str
        Path to the image.

    Returns
    -------
    receivedHeader : bytes array
        bytes array containing header data.
    headerLength : int
        Length of the header in bytes.
    receivedImage : bytes array
        Bytes array containing image data.

    """
    
    # Process the MPS Science Header (nBytesMpsHeader bytes --> 20 bytes)
    nBytesSync = 2              # Identify start of record --> 0x44 0x54
    nBytesSystemId = 1          # System ID of generated data
    nBytesDataId = 1            # Description of data type
    nBytesTotalLength = 4       # Total data packet length (including header)
    nBytesTime = 8              # Data acquisition start time
    nBytesSensorId = 1          # 
    nBytesHeaderVersion = 1     # 
    nBytesHeaderLength = 2      # Total length of header (in bytes). Take it as an offset to Data
    
    positionStartTotalLength = nBytesSync + nBytesSystemId + nBytesDataId
    positionStartHeaderLength = nBytesSync + nBytesSystemId + nBytesDataId + nBytesTotalLength + nBytesTime + nBytesSensorId + nBytesHeaderVersion
    
    file = open(Image_path,"rb")
    fullReceivedImage = file.read()
    mpsHeader = fullReceivedImage[0 : nBytesMpsHeader]
    
    totalLength = struct.unpack('i', mpsHeader[positionStartTotalLength: positionStartTotalLength + nBytesTotalLength])[0]
    headerLength = struct.unpack('H', mpsHeader[positionStartHeaderLength: positionStartHeaderLength + nBytesHeaderLength])[0]

    # print(headerLength)

    receivedHeader = fullReceivedImage[0:headerLength]
    receivedImage = fullReceivedImage[headerLength:totalLength-nBytesLengthTail]

    file.close()

    return receivedHeader, headerLength, receivedImage

def read_image(file, Head, headerLength):
    
    """
    This function reads a RAW image as a Numpy array
    """
    if Head['ImageType'] != 0:
        ImageType = Types_of_thumbnail_bn[Head['ImageType'] - first_type]
        print(Head['ImageType'],ImageType)

        PipelineConfig = Head['PipelineConfig']
        PipelineConfig = PipelineConfig * 2
        print(Head['PipelineConfig'],PipelineConfig)
    else:
        PipelineConfig = 1
        ImageType = 1

    width, height = Head['Roi_x_size']//ImageType//PipelineConfig, Head['Roi_y_size']//ImageType//PipelineConfig
    dtypef = '<i2'
    if PipelineConfig == 0:
        dtypef = '<i2'
    if PipelineConfig == 1:
        dtypef = '<i4'
    if PipelineConfig == 2:
        dtypef = '<i4'
    im_dummy = np.fromfile(file,dtype=dtypef, count=width*height, offset=headerLength)

    if ( im_dummy.size > 0 ):
        image = im_dummy.reshape([height, width]).astype(np.float)

    return image

def Read_Image(path, PlotFlag = False, printHeader_flag = False, vmin = 0, vmax = 4000):
    """
    This function reads TuMAG image and header information
    
    :param Path: Path to the image file (including image name)
    :type Path: str
    :param PlotFlag: If True the image gets represented. Default value = False
    :type PlotFlag: Boolean
    :param printHeader_flag: If True, Header info will be printed out. Default = False
    :type printHeader_flag: Boolean

    :return: header, image.
    :rtype: Dict, np.ndarray(float)
        
    """
     
    # Separate Header and image data 
    H, hl, I = HeadernImageSeparator(path) # Header, header_lenth and image
     
    # Read header
    Head = GetDatafromHeader(H) 
    # Read Image
    Im  = read_image(path, Head, hl) 
     
    if PlotFlag:
        plt.imshow(Im, cmap = 'Greys', vmin = vmin, vmax = vmax)
        plt.colorbar()
        plt.show()
         
    if printHeader_flag:
        for key in Head:
             print(key, ' : ', Head[key])
     
    return Head, Im




