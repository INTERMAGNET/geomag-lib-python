#!/usr/bin/python
'''Magnetic Data Handling routines

It is important to note that all MagneticData inherited classes should
technically follow an abstract class.  Abstaction has only been included
in Python 2.6 hence we will omit the notion of abstraction and rather concentrate
on respecting commonality between classes. All MagneticData inherited classes must
contain the following routines:

	read(filename)
	write(filename)
	
Additional routines can be added but that would not serve for the basis of commonlity
between all data storage formats.

Author:
Charles Blais, 2013, Geomagnetic Laboratory of Canada
'''
import numpy as np
		
# In newer Python versions, it is recommended to inherit from Object class
class MagneticData:
	'''Magnetic Data storage class
	
	Handles any magnetic data related routines.  It is mostly ment to abstract
	larger set of magnetic data storage formats like IAGA2002 and more.

	All data should be stored by components and datetime.  Supported components
	with datetime are:
		datetime,x,y,z,f,g,h,d,i
		
	Is is expected that all set and get routines utilize numpy arrays.

	Author:
	Charles Blais, 2013, Geomagnetic Laboratory of Canada
	'''
	def __init__(self, **kwargs):
		'''Constructor
		
		Keywords:
		- memorySave = save any component that can be generated from conversation methods
		into memory.  For example, if we request X from (H,D), the Y axis can also be calculated.
		Thus we save Y in memory if Y ever comes up.  Another example, is if someone request X multiple times.
		There is no need to recalculate each time. (default: True)
		- copy = copy magnetic data from this object (inherit)
		'''
		self.datetime_index = 'DATETIME'
		if 'copy' in kwargs: self.__magdata = kwargs['copy'].__magdata
		else: self.__magdata = dict()
		
		self.__memorySave = kwargs['memorySave'] if 'memorySave' in kwargs else True

		# define converting library in dictionary
		# this should reduce processing time during acquisition of multiple components
		self.__convertlib = {
			'X' : self.__getX,
			'Y' : self.__getY,
			'F' : self.__getF,
			'G' : self.__getG,
			'H' : self.__getH,
			'D' : self.__getD,
			'I' : self.__getI
		}
		
	def get(self, comp):
		'''Get magnetic data according to requested component.
		
		Magnetic data is stored in a dict which means that components
		are references by the key in the dict.  All components
		are stored in capital form.
		
		Arguments:
		comp - magnetic field component
		
		Returns:
		Array of magnetic data
		'''
		comp = comp.upper()
		if comp in self.__magdata:
			return self.__magdata[comp]
			
		# if the data is not stored in the format specified, see if we can generate
		# it from the other fields
		if comp in self.__convertlib:
			return self.__convertlib[comp]()
		else:
			raise Exception("Could not return component %s"%comp)
			
	def set(self, comp, data, append=False):
		'''Set magnetic data according to component.
		
		Magnetic data is stored in a dict which means that components
		are references by the key in the dict.  All components
		are stored in capital form.
		
		Arguments:
		comp - magnetic field component
		data - magnetic data
		'''
		if (comp.upper() in self.__magdata) and append:
			self.__magdata[comp.upper()] = np.append(self.__magdata[comp.upper()],data)
		else:
			self.__magdata[comp.upper()] = data
	
	def __getX(self):
		return self.__getXY()['X']
	def __getY(self):
		return self.__getXY()['Y']
	def __getXY(self):
		'''Get XY component field from HD'''
		if not 'H' in self.__magdata and not 'D' in self.__magdata:
			raise Exception("Can not generate X or Y component since H and D component must exist in dataset")
		bth = self.__magdata['D'] * ( np.pi/180/60 )
		data = dict()
		data['X'] = self.__magdata['H'] * np.cos( bth )
		data['Y'] = self.__magdata['H'] * np.sin( bth )
		if self.__memorySave:
			self.__magdata['X'] = data['X']
			self.__magdata['Y'] = data['Y']
		return data
		
	def __getH(self):
		return self.__getHD()['H']
	def __getD(self):
		return self.__getHD()['D']
	def __getHD(self):
		'''Get HD component field from XY'''
		if not 'X' in self.__magdata and not 'Y' in self.__magdata:
			raise Exception("Can not generate H or D component since X and Y component must exist in dataset")
		data = dict()
		data['H'] = np.sqrt(np.power(self.__magdata['X'],2) + np.power(self.__magdata['Y'],2))
		data['D'] = np.arctan2(self.__magdata['Y'], self.__magdata['X'])*(180/np.pi)*60
		if self.__memorySave:
			self.__magdata['H'] = data['H']
			self.__magdata['D'] = data['D']
		return data
		
	def __getF(self):
		'''Get F component from X,Y,Z and G'''
		if not 'G' in self.__magdata and not 'Z' in self.__magdata:
			raise Exception("Can not generate F component since G and Z component must exist in dataset")
		if not 'X' in self.__magdata or not 'Y' in self.__magdata:
			xydata = self.__getXY()
		else:
			xydata = {'X': self.__magdata['X'], 'Y' : self.__magdata['Y']}
		
		f = self.__magdata['G'] + np.sqrt(xydata['X'] + xydata['Y'] + self.__magdata['Z'])
		if self.__memorySave: self.__magdata['F'] = f
		return f
	
	def __getG(self):
		'''Get G component from X,Y,Z and F'''
		if not 'F' in self.__magdata and not 'Z' in self.__magdata:
			raise Exception("Can not generate G component since F and Z component must exist in dataset")
		if not 'X' in self.__magdata or not 'Y' in self.__magdata:
			xydata = self.__getXY()
		else:
			xydata = {'X': self.__magdata['X'], 'Y' : self.__magdata['Y']}
		
		g = self.__magdata['F'] - np.sqrt(xydata['X'] + xydata['Y'] + self.__magdata['Z'])
		if self.__memorySave: self.__magdata['G'] = g
		return g
		
	def __getI(self):
		'''Get G component from X,Y,Z and F'''
		if not 'Z' in self.__magdata:
			raise Exception("Can not generate I component since Z component must exist in dataset")
		if not 'X' in self.__magdata or not 'Y' in self.__magdata:
			xydata = self.__getXY()
		else:
			xydata = {'X': self.__magdata['X'], 'Y' : self.__magdata['Y']}
		
		i = np.arctan2(self.__magdata['Z'], np.sqrt(np.power(xydata['X'],2) + np.power(xydata['Y'],2)))
		if self.__memorySave: self.__magdata['I'] = i
		return i
		
	def cast(self, CastClass, **kwargs):
		'''Cast an object to another inherited class'''
		return CastClass(copy = self, **kwargs)
		
	def append(self, AppendClass):
		'''Append data from one class to another'''
		for idx,values in AppendClass.__magdata.iteritems():
			self.__magdata[idx] = np.append(self.__magdata[idx],values) if idx in self.__magdata else values

class Information:
	'''Information storage for magnetic data or observatory information
	
	The information stored in this object has been done so that all magnetic
	data routine uses identical information
	'''
	def __init__(self, **kwargs):
		'''Constructor
		
		Keywords:
		- copy = copy information from this object (inherit)
		'''
		if 'copy' in kwargs: self.__info = kwargs['copy'].__info
		else: self.__info = dict()
		
	def __setInfo(self, key, value):
		self.__info[key] = value
	def __getInfo(self, key):
		return self.__info[key] if key in self.__info else None
		
	def setIAGA(self, value): self.__setInfo('iaga',value)
	def getIAGA(self): return self.__getInfo('iaga')
	
	def setStationName(self, value): self.__setInfo('stationName', value)
	def getStationName(self): return self.__getInfo('stationName')
		
	def setInstitute(self, value): self.__setInfo('institute',value)
	def getInstitute(self): return self.__getInfo('institute')
	
	def setGIN(self, value): self.__setInfo('GIN',value)
	def getGIN(self): return self.__getInfo('GIN')
		
	def setReportedOrientation(self, value): self.__setInfo('reported_orientation',value)
	def getReportedOrientation(self): return self.__getInfo('reported_orientation')
	
	def setSensorOrientation(self, value): self.__setInfo('sensor_orientation',value)
	def getSensorOrientation(self): return self.__getInfo('sensor_orientation')
	
	# type: float, units: sec
	def setSamplingRate(self, value): self.__setInfo('samplingrate',value)
	def getSamplingRate(self): return self.__getInfo('samplingrate')
	
	# type: string
	def setDigitalSampling(self, value): self.__setInfo('digitalsampling',value)
	def getDigitalSampling(self): return self.__getInfo('digitalsampling')
	def setDataIntervalType(self, value): self.__setInfo('dataintervaltype',value)
	def getDataIntervalType(self): return self.__getInfo('dataintervaltype')
	
	def setDataType(self, value, shortform = False): self.__setInfo('datatype_short' if shortform else 'datatype',value)
	def getDataType(self, shortform = False): return self.__getInfo('datatype_short' if shortform else 'datatype')
			
	def setLongitude(self, value):
		if value < 0 or value > 360:  raise Exception("Longitude East must be betwen 0-360") 
		self.__setInfo('longitude', value)
	def getLongitude(self): return self.__getInfo('longitude')
	def setLatitude(self, value):
		if value < -90 or value > 90: raise Exception("Latitude must be between -90-90")
		self.__setInfo('latitude', value)
	def getLatitude(self): return self.__getInfo('latitude')
	def setAltitude(self, value): self.__setInfo('altitude', value)
	def getAltitude(self): return self.__getInfo('altitude')
	
	def setCoordinates(self, value):
		# Coordinates are stored in the format (lon,lat,alt)
		if len(value) != 3: raise Exception("Invalid coordinates, must be (lon,lat,alt)") 
		self.setLongitude(self, value[0])
		self.setLatitude(self, value[1])
		self.setAltitude(self, value[2])
	def getCoordinates(self):
		return (self.getLongitude(),self.getLatitude(),self.getAltitude())
		
	def setDecbas(self, value): self.__setInfo('decbas', value)
	def getDecbas(self): return self.__getInfo('decbas')

		
class IAGA2002(MagneticData, Information):
	'''IAGA2002 file  handling
	
	For detail description of IAGA2002 format, please visit
	http://www.ngdc.noaa.gov/IAGA/vdat/iagaformat.html
	
	NOTE: IAGA2002 data contains magnetic data and extra header
	information.  This is why it inherits properties from
	two parent classes.
	'''
	def __init__(self, **kwargs):
		'''Constructor'''
		MagneticData.__init__(self, **kwargs)
		Information.__init__(self, **kwargs)
		
		# internal storage
		self.__data = dict()
		
		# used by header reading
		self.__infoConvert = {
			'Source of Data' : self.setInstitute,
			'Station Name' : self.setStationName,
			'IAGA CODE' : self.setIAGA,
			'Geodetic Latitude' : self.setLatitude,
			'Geodetic Longitude' : self.setLongitude,
			'Elevation' : self.setAltitude,
			'Reported' : self.setReportedOrientation,
			'Sensor Orientation' : self.setSensorOrientation,
			'Digital Sampling' : self.setDigitalSampling,
			'Data Interval Type' : self.setDataIntervalType,
			'Data Type' : self.setDataType,
		}
		self.__infoFloatConvert = ('Geodetic Latitude','Geodetic Longitude','Elevation')
		
		# internal referenced variables
		self.__datetime_re = "^[0-9]{4}(-[0-9]{2}){2} [0-9]{2}(:[0-9]{2}){2}.[0-9]{3}"
		self.__datetime_format = "%Y-%m-%d %H:%M:%S.%f"
		
	def read(self, resource, **kwargs):
		'''Read IAGA2002 from file
		
		Keywords:
		- from_datetime = read only data starting at this time (not before)
		- to_datetime = read only data ending at this time (not beyond)
		'''
		import gzip
		
		# Open file for reading if its a string
		if isinstance(resource, basestring):
			file = gzip.open(resource, 'r') if resource.endswith(".gz") else open(resource, 'r')
		else: file = resource
		
		# in IAGA2002 files, the header is considered to stop
		# when the DATE header is reached
		isheader = True
		for line in file.readlines():
			# indicate end of header information
			if line.startswith("DATE"):
				isheader = False
				self.__initData()
			elif isheader: self.__decodeHeaderLine(line)
			else: self.__decodeDataLine(line, **kwargs)
			
		# Set the sampling rate using the data recovered
		# The sampling rate is the median time difference in seconds between samples
		sampling_rate = min([td.total_seconds() for td in np.diff(self.__data[self.datetime_index])])
		self.setSamplingRate(sampling_rate)
		
		# Insert data into parent object and clear memory of old object
		# They must be inserted as numpy arrays into the object
		for comp in self.__data:
			self.set(comp,self.__data[comp])
		self.__data = None # clear memory

		
		
		# Close file resource
		if isinstance(resource, basestring): file.close()
		
	def __initData(self):
		'''Initialize data array according to IAGA2002 format'''
		reported = self.getReportedOrientation()
		if reported is None:
			raise Exception("IAGA2002 file is missing mandatory Reported field")
		# initialize comp array data
		for comp in reported.upper():
			if not comp in self.__data:
				self.__data[comp] = np.array([])
		# datetime field must also exist
		if not self.datetime_index in self.__data:
			self.__data[self.datetime_index] = np.array([])
		
	def __decodeHeaderLine(self, line):
		'''Decode IAGA2002 header line
		
		Decode IAGA2002 header line and store information in memory.  When writting
		information to IAGA2002 files, some header information is required.  We use the 
		setInfo routine in that condition.
		
		By definition, content lavel starts at column 2 and description at column 25.
		
		Arguments:
		line - line to decode
		'''
		# ignore lines that are not header or that are commented
		if len(line) < 70: return
		if not line[69].endswith("|"): return
		if line.strip().startswith("#"): return
		
		header = line[:24].strip()
		value = line[24:line.index("|")].strip()
		if len(value) == 0: return #no value associated with the field
		if header in self.__infoConvert:
			self.__infoConvert[header](float(value) if header in self.__infoFloatConvert else value)
	
	def __decodeDataLine(self, line, **kwargs):
		'''Decode IAGA2002 data line
		
		Decoded IAGA2002 data line.  The data line should be in the format:
		2013-09-01 00:00:00.000 244     19056.50   4862.20 -28999.40  35063.00
		
		By definition, the DATE should always be in the ISO format YYYY-mm-dd
		The TIME should always be in the ISO format hh:mm:ss.sss
		
		D and I are always reprensented as minutes of arc.  All others are nT.
		
		Arguments:
		line - line to decode
		
		Keywords:
		- from_datetime = read only data starting at this datetime (not before)
		- to_datetime = read only data ending at this datetime (not beyond)
		'''
		import re
		from datetime import datetime
		if not re.search(self.__datetime_re,line): return
		
		# We have a valid line, we split the line by spaces (columns) and store their
		# value in object memory temporary until added to parent object.
		line = line.split()
		# the first 3 elements are reserved from time stamping (date,time,doy)
		# the other elements are the reported data
		if len(line) < (3 + len(self.getReportedOrientation())): return  #incomplete line
		
		ts = datetime.strptime("%s %s"%(line[0],line[1]),self.__datetime_format)
		if 'from_datetime' in kwargs:
			if ts < kwargs['from_datetime']: return
		if 'to_datetime' in kwargs:
			if ts > kwargs['to_datetime']: return
		
		self.__data[self.datetime_index] = np.append(self.__data[self.datetime_index],ts)
		reported = self.getReportedOrientation().upper()
		for idx in xrange(len(reported)):
			# validate value, ignore anything above unreported data (88888)
			# we insert np.nan as opposed to None since we are dealing with np arrays
			value = float(line[idx+3]) if float(line[idx + 3]) < 88888 else np.nan
			self.__data[reported[idx]] = np.append(self.__data[reported[idx]],value)
			
	def write(self, resource):
		'''Write IAGA2002 file
		
		IAGA2002 files all start with a header lines of 70 characters long.  The comments
		often found in IAGA2002 are not currently reproduced in the output.  We rather
		put a fix header comment regarding the web service and INTERMAGNET.
		
		NOTE: IAGA2002 files require 4 reported components to be writte
		'''
		if len(self.getReportedOrientation()) != 4:
			raise Exception("Reported orientation must be have 4 components (IAGA2002 rule)")
		
		if isinstance(resource, basestring):
			file = gzip.open(resource,'w') if resource.endswith(".gz") else open(resource,'w')
		else: file = resource
		
		self.__writeHeader(file)
		comps = self.getReportedOrientation()
		for i,value in enumerate(self.get(self.datetime_index)):
			milli = "%03d"%(value.microsecond/1000)
			data = [(99999.00 if np.isnan(self.get(comp)[i]) else self.get(comp)[i]) for comp in comps]
			data.insert(0,value.strftime("%Y-%m-%d %H:%M:%S."+milli+" %j   "))
			file.write("%s %9.2f %9.2f %9.2f %9.2f\r\n"%tuple(data))
		
		# Close file resource
		if isinstance(resource, basestring): file.close()
		
	def __writeHeader(self, file):
		'''Write IAGA2002 header'''
		headers = [["Format", "IAGA2002"],
			["Source of Data", self.getInstitute()],
			["Station Name", self.getStationName()],
			["IAGA CODE", self.getIAGA()],
			["Geodetic Latitude", "%.3f"%self.getLatitude()],
			["Geodetic Longitude", "%.3f"%self.getLongitude()],
			["Elevation", "%d"%self.getAltitude()],
			["Reported", self.getReportedOrientation()],
			["Sensor Orientation", self.getSensorOrientation()],
			["Digital Sampling", self.getDigitalSampling()],
			["Data Interval Type", self.getDataIntervalType()],
			["Data Type", self.getDataType()]]
			
		for header,content in headers:
			if not isinstance(content,basestring): content = ''
			file.write(" %-23s%-45s|\r\n"%(header,content))
		# add fixed IAGA2002 header
		file.write(" # This data file was created for the INTERMAGNET web service.       |\r\n")
		file.write(" # CONDITIONS OF USE: The Conditions of Use for data provided        |\r\n")
		file.write(" # through INTERMAGNET and acknowledgement templates can be found    |\r\n")
		file.write(" # at www.intermagnet.org                                            |\r\n")
		
		iaga = self.getIAGA()
		reported = self.getReportedOrientation()
		
		header = "DATE       TIME         DOY     "
		for comp in reported:
			header += "%3s%1s   "%(iaga,comp)
			header += "|\r\n" if comp == reported[3] else "   "
		file.write(header)

		
class IMF(MagneticData, Information):
	'''IMFv1.22 file  handling
	
	For detail description of IMF format, please visit
	http://www.intermagnet.org/data-donnee/formats/imfv122-eng.php
	
	NOTE: IMFv1.22 data contains magnetic data and extra header
	information.  This is why it inherits properties from
	two parent classes.
	'''
	def __init__(self, **kwargs):
		'''Constructor'''
		MagneticData.__init__(self, **kwargs)
		Information.__init__(self, **kwargs)
		
	def write(self, resource):
		'''Write IMFv1.22 to file
		
		It is important to note that the data stored in memory must be for only
		a days worth and can not exceed.
		
		IMPORTANT! We assume the data is stored chronologically and is minute data.
		'''
		from datetime import datetime
		
		if len(self.getReportedOrientation()) != 4:
			raise Exception("Reported orientation must be have 4 components (IMFv1.22 rule)")
		
		if isinstance(resource, basestring):
			file = gzip.open(resource,'w') if resource.endswith(".gz") else open(resource,'w')
		else: file = resource
		
		# determine the timestamp to print for the header block
		# to do so, we loop over the datetimes for each hour block
		# Python uses pointers to dict elements hence it is not necessary
		# to extract all components before the loop
		dts = self.get(self.datetime_index)
		
		headerDateString = dts[0].strftime("%b%d%y %j").upper()
		
		iaga = self.getIAGA()
		comps = self.getReportedOrientation().upper()
		type = self.getDataType()
		if not type is None: type = type[0].upper() if len(type) > 1 else type.upper()
		else: type = ''
		gin = '' if self.getGIN() is None else self.getGIN().upper()
		colat = (90 - self.getLatitude())*10
		lon = self.getLongitude()*10
		decbas = 0 if self.getDecbas() is None else self.getDecbas()
		
		# constants
		formatHeader = "%3s %s %02d %4s %1s %3s %04d%04d %06d RRRRRRRRRRRRRRRR\r\n"
		formatBody = "%7d %7d %7d %6d"
		nullDataSet = (999999,999999,999999,99999)
		
		# we start the loop at the first hour in the file
		idx = 0
		for hour in xrange(24):
			file.write(formatHeader%(iaga,headerDateString,hour,comps,type,gin,colat,lon,decbas))
			for minute in xrange(60):
				# determine if the next timestamp is hours, if not skip
				ts = datetime(dts[0].year,dts[0].month,dts[0].day,hour,minute)
				if idx >= len(dts): data = nullDataSet #out of range
				elif ts < dts[idx]: data = nullDataSet
				elif ts == dts[idx]:
					data = tuple()
					for i in xrange(len(comps)):
						if np.isnan(self.get(comps[i])[idx]): temp = 999999 if i == 4 else 9999999
						else: temp = self.get(comps[i])[idx]*10 
						data = data + (temp,)
					idx += 1
				else:
					raise Exception("Unexecpted value.  Expected minute resolution data stored chronologically")
				file.write(formatBody%data)
				file.write("\r\n" if minute%2 else "  ")
		
		# Close file resource
		if isinstance(resource, basestring): file.close()
		
class CDF(MagneticData, Information):
	'''NASA CDF file handling
	
	For detail description of CDF format, please visit
	http://cdf.gsfc.nasa.gov/
	
	For details on INTERMAGNET CDF format, please visit:
	
	The following requires the SpacePy package.  Note that it 
	is tricky to install.  Follow instructions on the web and copy
	compiled Python code under python library (it is not copied over).
	
	NOTE: CDF data contains magnetic data and extra header
	information.  This is why it inherits properties from
	two parent classes.
	'''
	def __init__(self, **kwargs):
		'''Constructor'''
		MagneticData.__init__(self, **kwargs)
		Information.__init__(self, **kwargs)
		
	def write(self, resource):
		'''Write CDF to file'''
		from spacepy import pycdf
		from datetime import datetime
		
		nullValue = 99999.0
		
		cdf = pycdf.CDF(resource, '')
		cdf.compress(pycdf.const.GZIP_COMPRESSION, 9)
		self.__writeAttrs(cdf)
		comps = self.getReportedOrientation()
		for i in xrange(len(comps)):
			field = "GeomagneticFieldElement%d"%(i+1)
			data = [ (nullValue if np.isnan(value) else value) for value in self.get(comps[i]).tolist() ]
			zVariable = cdf.new(field, data=data, type=pycdf.const.CDF_DOUBLE)
			cdf[field].attrs['FIELDNAME'] = "Geomagnetic Field Element %d"%(i+1)
			cdf[field].attrs['VALIDMIN'] = -79999.0
			cdf[field].attrs['VALIDMAX'] = 79999.0
			if comps[i] in ('X','Y','Z','F','H','G'): units = 'nT'
			elif comps[i] in ('D','I'): units = 'Minutes of arc'
			else: units = ""
			cdf[field].attrs['UNITS'] = units
			cdf[field].attrs['FILLVAL'] = nullValue
			cdf[field].attrs['StartDate'] = self.get(self.datetime_index)[0].isoformat()
			cdf[field].attrs['StartDateEpoch'] = self.get(self.datetime_index)[0]
			cdf[field].attrs['SampPer'] = self.getSamplingRate()
			cdf[field].attrs['ElemRec'] = comps[i]
			cdf[field].attrs['OrigFreq'] = 99999.0
		
	def __writeAttrs(self, cdf):
		cdf.attrs['Title'] = "Geomagnetic observatory data"
		cdf.attrs['FormatDescription'] = "INTERMAGNET CDF Format"
		cdf.attrs['FormatVersion'] = "1.0"
		cdf.attrs['TermsOfUse'] = '''CONDITIONS OF USE FOR DATA PROVIDED THROUGH INTERMAGNET:
The data made available through INTERMAGNET are provided for
your use and are not for commercial use or sale or distribution
to third parties without the written permission of the institute
(http://www.intermagnet.org/Institutes_e.html) operating
the observatory. Publications making use of the data
should include an acknowledgment statement of the form given below.
A citation reference should be sent to the INTERMAGNET Secretary
(secretary@intermagnet.org) for inclusion in a publications list
on the INTERMAGNET website.

     ACKNOWLEDGEMENT OF DATA FROM OBSERVATORIES
     PARTICIPATING IN INTERMAGNET
We offer two acknowledgement templates. The first is for cases
where data from many observatories have been used and it is not
practical to list them all, or each of their operating institutes.
The second is for cases where research results have been produced
using a smaller set of observatories.

     Suggested Acknowledgement Text (template 1)
The results presented in this paper rely on data collected
at magnetic observatories. We thank the national institutes that
support them and INTERMAGNET for promoting high standards of
magnetic observatory practice (www.intermagnet.org).

     Suggested Acknowledgement Text (template 2)
The results presented in this paper rely on the data
collected at <observatory name>. We thank <institute name>,
for supporting its operation and INTERMAGNET for promoting high
standards of magnetic observatory practice (www.intermagnet.org).
'''
		cdf.attrs['Institution'] = self.getInstitute()
		cdf.attrs['Source'] = "INTERMAGNET (Ottawa GIN)"
		cdf.attrs['History'] = "Unknown"
		cdf.attrs['References'] = "www.intermagnet.org"
		cdf.attrs['ObservatoryName'] = self.getStationName()
		cdf.attrs['IagaCode'] = self.getIAGA()
		cdf.attrs['Latitude'] = self.getLatitude()
		cdf.attrs['Longitude'] = self.getLongitude()
		cdf.attrs['Elevation'] = self.getAltitude()
		cdf.attrs['VectorSensOrient'] = 'Unknown' if self.getSensorOrientation() is None else self.getSensorOrientation()
		dataType = "variation" if self.getDataType() is None else self.getDataType().lower()
		cdf.attrs['BaselineType'] = "None" if dataType == "variation" else dataType
		cdf.attrs['PublicationState'] = "raw"
		cdf.attrs['StandardsConformance'] = "None"
		cdf.attrs['PublicationDate'] = ""
		
		
		
class CSV(MagneticData, Information):
	'''Magnetic CSV
	
	Magnetic CSV is a ASCII text table where each column is divided by commas.  The following
	library allows the delimiter to be altered.  The table columns are as followed:
	
	- datetime
	- [component1]
	- [component2]
	- ...
	
	It is important to note that each component columns are actually titles by X,Y,Z,F,...
	
	Datetime column are specified in the standard ISO format YYYY-mm-dd HH:MM:SS
	
	IMPORTANT! It is designed in a way that the response is readable in its raw form and lightweight.  This way, the CSV
	response from a web service is as small as possible.  Any additional information should be determined
	seperatly (data type, sampling rate,...).
	'''
	def __init__(self, **kwargs):
		'''Constructor
		
		Keywords:
		delimiter - change the delimiter character (default: ,)
		lineterminator - change the end of line character (default: \r\n)
		'''
		MagneticData.__init__(self, **kwargs)
		Information.__init__(self, **kwargs)
		self.__delimiter = kwargs['delimiter'] if 'delimiter' in kwargs else ","
		self.__lineterminator = kwargs['lineterminator'] if 'lineterminator' in kwargs else "\r\n"
	
	def write(self, resource):
		'''Write to CSV file'''
		import csv
		if isinstance(resource, basestring):
			file = gzip.open(resource,'w') if resource.endswith(".gz") else open(resource,'w')
		else: file = resource
		csvwriter = csv.writer(file, delimiter=self.__delimiter, lineterminator=self.__lineterminator)
		
		comps = self.getReportedOrientation().upper()
		header = [comp for comp in comps]
		header.insert(0,'datetime')
		csvwriter.writerow(header)
		for i,value in enumerate(self.get(self.datetime_index)):
			data = [(None if np.isnan(self.get(comp)[i]) else self.get(comp)[i]) for comp in comps]
			data.insert(0,value)
			csvwriter.writerow(data)
		# Close file resource
		if isinstance(resource, basestring): file.close()
		
class JSON(MagneticData, Information):
	'''Magnetic JSON
	
	No standards have yet been formulized for magnetic data.  For the moment, we do
	not concentrate on either JSON-LD or any other format.  We rather simplify, the process
	by returning response in the following format:
	'''
	def __init__(self, **kwargs):
		'''Constructor'''
		MagneticData.__init__(self, **kwargs)
		Information.__init__(self, **kwargs)
		
	def write(self, resource):
		'''Write JSON response
		
		More then often, we expect the resource field to be stdout.  We do allow the JSON
		repsonse to be added to files even though it serves little use.
		'''
		from calendar import timegm
		import json
		if isinstance(resource, basestring):
			file = gzip.open(resource,'w') if resource.endswith(".gz") else open(resource,'w')
		else: file = resource
		
		info = {
			'station':{
				'IAGA' : self.getIAGA(),
				'coordinates' : self.getCoordinates(),
				'institute':self.getInstitute(),
				'name':self.getStationName(),
				'GIN':self.getGIN()
			},
			'reported':self.getReportedOrientation().upper(),
			'dataType':self.getDataType()
		}
		
		comps = self.getReportedOrientation().upper()
		data = {'datetime':{'start':self.get(self.datetime_index)[0].strftime("%Y-%m-%d %H:%M:%S"),'timestamp':timegm(self.get(self.datetime_index)[0].utctimetuple()),'interval':self.getSamplingRate()}}
		# json does not accept Numpy arrays, it must therefore be converted to a list
		for comp in comps:
			data[comp] =  [ (None if np.isnan(value) else value) for value in self.get(comp).tolist() ]
		jsonText = {"@info":info,"data":data}
		# set JSON float encoder to reduce number resolution
		json.encoder.FLOAT_REPR = lambda o: format(o, '.2f')
		file.write(json.dumps(jsonText))
		
		# Close file resource
		if isinstance(resource, basestring): file.close()

		
if __name__ == "__main__":
	ReadObj = IAGA2002()
	ReadObj.read( "ott20130101vmin.min" )
	ReadObj2 = IAGA2002()
	ReadObj2.read( "ott20130102vmin.min" )
	ReadObj.append( ReadObj2 )
	
	WriteObj = ReadObj.cast(MagneticCSV)
	import sys
	WriteObj.write(sys.stdout)
