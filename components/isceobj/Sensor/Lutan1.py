#!/usr/bin/python3

# Reader for Lutan-1 SLC data
# Used Sentinel1.py and ALOS.py as templates
# Author: Bryan Marfito, EOS-RS


import os
import glob
import numpy as np
import xml.etree.ElementTree as ET
import datetime
import isce
import isceobj
from isceobj.Planet.Planet import Planet
from iscesys.Component.Component import Component
from isceobj.Sensor.Sensor import Sensor
from isceobj.Scene.Frame import Frame
from isceobj.Orbit.Orbit import StateVector, Orbit
from isceobj.Planet.AstronomicalHandbook import Const
from iscesys.DateTimeUtil.DateTimeUtil import DateTimeUtil as DTUtil
from isceobj.Orbit.OrbitExtender import OrbitExtender
from osgeo import gdal
import warnings
from scipy.interpolate import UnivariateSpline
from scipy.signal import savgol_filter
import pdb

lookMap = { 'RIGHT' : -1,
            'LEFT' : 1}

# Antenna dimensions 9.8 x 3.4 m
antennaLength = 9.8

XML = Component.Parameter('xml',
        public_name = 'xml',
        default = None,
        type = str,
        doc = 'Input XML file')


TIFF = Component.Parameter('tiff',
                            public_name ='tiff',
                            default = None,
                            type=str,
                            doc = 'Input image file')

ORBIT_FILE = Component.Parameter('orbitFile',
                            public_name ='orbitFile',
                            default = None,
                            type=str,
                            doc = 'Orbit file')


class Lutan1(Sensor):

    "Class for Lutan-1 SLC data"
    
    family = 'l1sm'
    logging_name = 'isce.sensor.Lutan1'

    parameter_list = (TIFF, ORBIT_FILE) + Sensor.parameter_list

    def __init__(self, name = ''):
        super(Lutan1,self).__init__(self.__class__.family, name=name)
        self.frame = Frame()
        self.frame.configure()
        self._xml_root = None
        self.doppler_coeff = None

    def parse(self):
        xmlFileName = self.tiff[:-4] + "meta.xml"
        self.xml = xmlFileName

        with open(self.xml, 'r') as fid:
            xmlstr = fid.read()
        
        self._xml_root = ET.fromstring(xmlstr)
        self.populateMetadata()
        fid.close()

        if self.orbitFile:
            # Check if orbit file exists or not
            if os.path.isfile(self.orbitFile) == True:
                orb = self.extractOrbit()
                self.frame.orbit.setOrbitSource(os.path.basename(self.orbitFile))
            else:
                pass
        else:
            warnings.warn("WARNING! No orbit file found. Orbit information from the annotation file is used for processing.")
            orb = self.extractOrbitFromAnnotation()
            self.frame.orbit.setOrbitSource(os.path.basename(self.xml))
            self.frame.orbit.setOrbitSource('Annotation')

        for sv in orb:
            self.frame.orbit.addStateVector(sv)

    def convertToDateTime(self,string):
        dt = datetime.datetime.strptime(string,"%Y-%m-%dT%H:%M:%S.%f")
        return dt


    def grab_from_xml(self, path):
        try:
            res = self._xml_root.find(path).text
        except:
            raise Exception('Tag= %s not found'%(path))

        if res is None:
            raise Exception('Tag = %s not found'%(path))
        
        return res
    

    def populateMetadata(self):
        mission = self.grab_from_xml('generalHeader/mission')
        polarization = self.grab_from_xml('productInfo/acquisitionInfo/polarisationMode')
        frequency = float(self.grab_from_xml('instrument/radarParameters/centerFrequency'))
        passDirection = self.grab_from_xml('productInfo/missionInfo/orbitDirection')
        rangePixelSize = float(self.grab_from_xml('productInfo/imageDataInfo/imageRaster/columnSpacing'))
        azimuthPixelSize = float(self.grab_from_xml('productInfo/imageDataInfo/imageRaster/rowSpacing'))
        rangeSamplingRate = Const.c/(2.0*rangePixelSize)

        prf = float(self.grab_from_xml('instrument/settings/settingRecord/PRF'))
        lines = int(self.grab_from_xml('productInfo/imageDataInfo/imageRaster/numberOfRows'))
        samples = int(self.grab_from_xml('productInfo/imageDataInfo/imageRaster/numberOfColumns'))

        startingRange = float(self.grab_from_xml('productInfo/sceneInfo/rangeTime/firstPixel'))*Const.c/2.0
        #slantRange = float(self.grab_from_xml('productSpecific/complexImageInfo/'))
        incidenceAngle = float(self.grab_from_xml('productInfo/sceneInfo/sceneCenterCoord/incidenceAngle'))
        dataStartTime = self.convertToDateTime(self.grab_from_xml('productInfo/sceneInfo/start/timeUTC'))
        dataStopTime = self.convertToDateTime(self.grab_from_xml('productInfo/sceneInfo/stop/timeUTC'))
        pulseLength = float(self.grab_from_xml('processing/processingParameter/rangeCompression/chirps/referenceChirp/pulseLength'))
        pulseBandwidth = float(self.grab_from_xml('processing/processingParameter/rangeCompression/chirps/referenceChirp/pulseBandwidth'))
        chirpSlope = pulseBandwidth/pulseLength

        if self.grab_from_xml('processing/processingParameter/rangeCompression/chirps/referenceChirp/chirpSlope') == "DOWN":
            chirpSlope = -1.0 * chirpSlope
        else:
            pass

        # Check for satellite's look direction
        if self.grab_from_xml('productInfo/acquisitionInfo/lookDirection') == "LEFT":
            lookSide = lookMap['LEFT']
            print("Look direction: LEFT")
        else:
            lookSide = lookMap['RIGHT']
            print("Look direction: RIGHT")

        processingFacility = self.grab_from_xml('productInfo/generationInfo/level1ProcessingFacility')

        # Platform parameters
        platform = self.frame.getInstrument().getPlatform()
        platform.setPlanet(Planet(pname='Earth'))
        platform.setMission(mission)
        platform.setPointingDirection(lookSide)
        platform.setAntennaLength(antennaLength)

        # Instrument parameters
        instrument = self.frame.getInstrument()
        instrument.setRadarFrequency(frequency)
        instrument.setPulseRepetitionFrequency(prf)
        instrument.setPulseLength(pulseLength)
        instrument.setChirpSlope(chirpSlope)
        instrument.setIncidenceAngle(incidenceAngle)
        instrument.setRangePixelSize(rangePixelSize)
        instrument.setRangeSamplingRate(rangeSamplingRate)
        instrument.setPulseLength(pulseLength)

        # Frame parameters
        self.frame.setSensingStart(dataStartTime)
        self.frame.setSensingStop(dataStopTime)
        self.frame.setProcessingFacility(processingFacility)

        # Two-way travel time 
        diffTime = DTUtil.timeDeltaToSeconds(dataStopTime - dataStartTime) / 2.0
        sensingMid = dataStartTime + datetime.timedelta(microseconds=int(diffTime*1e6))
        self.frame.setSensingMid(sensingMid)
        self.frame.setPassDirection(passDirection)
        self.frame.setPolarization(polarization)
        self.frame.setStartingRange(startingRange)
        self.frame.setFarRange(startingRange +  (samples - 1) * rangePixelSize)
        self.frame.setNumberOfLines(lines)
        self.frame.setNumberOfSamples(samples)

        return


    def extractOrbit(self):

        '''
        Extract orbit information from the orbit file
        '''
        orb = Orbit()
        orb.configure()

        # I based the margin on the data that I have.
        # Lutan-1 position and velocity sampling frequency is 1 Hz
        margin = datetime.timedelta(minutes=30.0)
        tstart = self.frame.getSensingStart() - margin
        tend = self.frame.getSensingStop() + margin

        file_ext = os.path.splitext(self.orbitFile)[1].lower()

        if file_ext == '.xml':
            try:
                fp = open(self.orbitFile, 'r')
            except IOError as strerr:
                print("IOError: %s" % strerr)
            
            _xml_root = ET.ElementTree(file=fp).getroot()
            node = _xml_root.find('Data_Block/List_of_OSVs')
            
            for child in node:
                timestamp = self.convertToDateTime(child.find('UTC').text)
                if (timestamp >= tstart) and (timestamp <= tend):
                    pos = []
                    vel = []
                    for tag in ['VX', 'VY', 'VZ']:
                        vel.append(float(child.find(tag).text))

                    for tag in ['X', 'Y', 'Z']:
                        pos.append(float(child.find(tag).text))

                    vec = StateVector()
                    vec.setTime(timestamp)
                    vec.setPosition(pos)
                    vec.setVelocity(vel)
                    orb.addStateVector(vec)

            fp.close()

        elif file_ext == '.txt':
            with open(self.orbitFile, 'r') as fid:
                for line in fid:
                    if not line.startswith('#'):
                        break
                
                for line in fid:
                    fields = line.split()
                    if len(fields) >= 13:
                        year = int(fields[0])
                        month = int(fields[1])
                        day = int(fields[2])
                        hour = int(fields[3])
                        minute = int(fields[4])
                        second = float(fields[5])
                        
                        int_second = int(second)
                        microsecond = int((second - int_second) * 1e6)
                        # Convert to datetime   
                        timestamp = datetime.datetime(year, month, day, hour, minute, int_second, microsecond)
                        
                        if (timestamp >= tstart) and (timestamp <= tend):
                            pos = [float(fields[6]), float(fields[7]), float(fields[8])]
                            vel = [float(fields[9]), float(fields[10]), float(fields[11])]
                            vec = StateVector()
                            vec.setTime(timestamp)
                            vec.setPosition(pos)
                            vec.setVelocity(vel)
                            orb.addStateVector(vec)
        else:
            raise Exception("Unsupported orbit file extension: %s" % file_ext)
        return orb

    def orbit_filter_iterative(self, time_seconds, pos_data, vel_data, n_key_points=5):
        '''
        轨道数据迭代滤波函数
        参数:
        time_seconds: 时间序列（相对于第一个点的秒数）
        pos_data: 位置数据 shape=(n,3)
        vel_data: 速度数据 shape=(n,3)
        n_key_points: 关键点数量，默认5个
        '''
        time_array = np.array(time_seconds)
        pos_array = np.array(pos_data)
        vel_array = np.array(vel_data)
        filtered_pos = pos_array.copy()
        filtered_vel = vel_array.copy()
        
        for iteration in range(5):  # 进行5次迭代
            # 选择均匀分布的关键点
            n_points = len(time_array)
            key_indices = np.linspace(0, n_points-1, n_key_points, dtype=int)
            key_times = time_array[key_indices]
            key_pos = filtered_pos[key_indices]
            
            # 对关键点进行多项式拟合
            smoothed_pos = np.zeros_like(filtered_pos)
            for i in range(3):  # x,y,z三个方向
                # 位置使用3次多项式
                poly_deg = 3
                
                # 先用关键点拟合一个初始多项式
                coeffs = np.polyfit(key_times, key_pos[:,i], deg=poly_deg)
                # 用初始多项式拟合所有点
                smoothed_pos[:,i] = np.polyval(coeffs, time_array)
                # 再次对所有点进行拟合
                coeffs = np.polyfit(time_array, smoothed_pos[:,i], deg=poly_deg)
                smoothed_pos[:,i] = np.polyval(coeffs, time_array)
            
            smoothed_vel = np.zeros_like(filtered_vel)
            key_vel = filtered_vel[key_indices]
            for i in range(3):
                # 速度使用2次多项式
                vel_poly_deg = 2
                coeffs = np.polyfit(key_times, key_vel[:,i], deg=vel_poly_deg)
                smoothed_vel[:,i] = np.polyval(coeffs, time_array)
                coeffs = np.polyfit(time_array, smoothed_vel[:,i], deg=vel_poly_deg)
                smoothed_vel[:,i] = np.polyval(coeffs, time_array)
            
            filtered_pos = smoothed_pos
            filtered_vel = smoothed_vel
        
        # 将numpy数组转换为Python列表
        return filtered_pos.tolist(), filtered_vel.tolist()
    
    # def extractOrbitFromAnnotation(self):

    #     '''
    #     Extract orbit information from xml annotation
    #     WARNING! Only use this method if orbit file is not available
    #     '''

    #     try:
    #         fp = open(self.xml, 'r')
    #     except IOError as strerr:
    #         print("IOError: %s" % strerr)

    #     _xml_root = ET.ElementTree(file=fp).getroot()
    #     node = _xml_root.find('platform/orbit')
    #     countNode = len(list(_xml_root.find('platform/orbit')))

    #     frameOrbit = Orbit()
    #     frameOrbit.setOrbitSource('Header')
    #     margin = datetime.timedelta(minutes=30.0)
    #     tstart = self.frame.getSensingStart() - margin
    #     tend = self.frame.getSensingStop() + margin
    #     for k in range(1,countNode):
    #         timestamp = self.convertToDateTime(node.find('stateVec[{}]/timeUTC'.format(k)).text)
    #         if (timestamp >= tstart) and (timestamp <= tend):
    #             pos = [float(node.find('stateVec[{}]/posX'.format(k)).text), float(node.find('stateVec[{}]/posY'.format(k)).text), float(node.find('stateVec[{}]/posZ'.format(k)).text)]
    #             vel = [float(node.find('stateVec[{}]/velX'.format(k)).text), float(node.find('stateVec[{}]/velY'.format(k)).text), float(node.find('stateVec[{}]/velZ'.format(k)).text)]

    #             vec = StateVector()
    #             vec.setTime(timestamp)
    #             vec.setPosition(pos)
    #             vec.setVelocity(vel)
    #             frameOrbit.addStateVector(vec)
        
    #     fp.close()
    #     return frameOrbit
    def extractOrbitFromAnnotation(self):
        '''
        从xml注释中提取轨道信息并进行滤波
        如果没有精轨数据，使用此方法
        '''
        try:
            fp = open(self.xml, 'r')
        except IOError as strerr:
            print("IOError: %s" % strerr)
    
        _xml_root = ET.ElementTree(file=fp).getroot()
        node = _xml_root.find('platform/orbit')
        countNode = len(list(_xml_root.find('platform/orbit')))
    
        frameOrbit = Orbit()
        frameOrbit.setOrbitSource('Header')
        margin = datetime.timedelta(minutes=30.0)
        tstart = self.frame.getSensingStart() - margin
        tend = self.frame.getSensingStop() + margin
        
        # 收集所有轨道数据点
        timestamps = []
        positions = []
        velocities = []
        
        for k in range(1,countNode):
            timestamp = self.convertToDateTime(node.find('stateVec[{}]/timeUTC'.format(k)).text)
            if (timestamp >= tstart) and (timestamp <= tend):
                pos = [float(node.find('stateVec[{}]/posX'.format(k)).text), 
                      float(node.find('stateVec[{}]/posY'.format(k)).text), 
                      float(node.find('stateVec[{}]/posZ'.format(k)).text)]
                vel = [float(node.find('stateVec[{}]/velX'.format(k)).text), 
                      float(node.find('stateVec[{}]/velY'.format(k)).text), 
                      float(node.find('stateVec[{}]/velZ'.format(k)).text)]
                
                timestamps.append(timestamp)
                positions.append(pos)
                velocities.append(vel)
        
        fp.close()
        
        # 计算相对时间（秒）
        time_seconds = [(t - timestamps[0]).total_seconds() for t in timestamps]
        
        # 进行轨道滤波
        filtered_pos, filtered_vel = self.orbit_filter_iterative(time_seconds, positions, velocities)
        
        # 将滤波后的结果转换为轨道状态向量
        for i, timestamp in enumerate(timestamps):
            vec = StateVector()
            vec.setTime(timestamp)
            # 确保使用Python列表
            vec.setPosition(filtered_pos[i])
            vec.setVelocity(filtered_vel[i])
            frameOrbit.addStateVector(vec)
        
        return frameOrbit
    
    def extractImage(self):
        self.parse()
        width = self.frame.getNumberOfSamples()
        lgth = self.frame.getNumberOfLines()
        src = gdal.Open(self.tiff.strip(), gdal.GA_ReadOnly)

        # Band 1 as real and band 2 as imaginary numbers
        # Confirmed by Zhang Yunjun
        band1 = src.GetRasterBand(1)
        band2 = src.GetRasterBand(2)
        cJ = np.complex64(1.0j)

        fid = open(self.output, 'wb')
        for ii in range(lgth):
            # Combine the real and imaginary to make
            # them in to complex numbers
            real = band1.ReadAsArray(0,ii,width,1)
            imag = band2.ReadAsArray(0,ii,width,1)
            # Data becomes np.complex128 after combining them
            data = real + (cJ * imag)
            data.tofile(fid)

        fid.close()
        real = None
        imag = None
        src = None
        band1 = None
        band2 = None

        ####
        slcImage = isceobj.createSlcImage()
        slcImage.setByteOrder('l')
        slcImage.setFilename(self.output)
        slcImage.setAccessMode('read')
        slcImage.setWidth(self.frame.getNumberOfSamples())
        slcImage.setLength(self.frame.getNumberOfLines())
        slcImage.setXmin(0)
        slcImage.setXmax(self.frame.getNumberOfSamples())
        self.frame.setImage(slcImage)

    def extractDoppler(self):
        '''
        Set doppler values to zero since the metadata doppler values are unreliable.
        Also, the SLC images are zero doppler.
        '''
        dop = [0., 0., 0.]

        ####For insarApp
        quadratic = {}
        quadratic['a'] = dop[0] / self.frame.getInstrument().getPulseRepetitionFrequency()
        quadratic['b'] = 0.
        quadratic['c'] = 0.

        print("Average doppler: ", dop)
        self.frame._dopplerVsPixel = dop

        return quadratic
