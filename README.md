geomag-lib-python
=================

General library for magnetic file writing/reading in Python

#magnetic.py

Author: Charles Blais, Natural Resources Canada, Government of Canada

##Description
Python library for read/write routines of geomagnetic data. Supports the following formats:
- IAGA2002  (read/write)
- IMFv122	(write)

It also supports non-standardized formats:
- CSV		(write)
- CDF		(write)
- JSON		(write)

This library was designed to translate IAGA2002 magnetic data to other formats.  Plans are to expend the library for reading/writing of multiple formats.

##Requirements

- Python 2.6+
- NumPy (http://www.numpy.org)
- SpacePy (http://spacepy.lanl.gov/) - for CDF format only

##How to Use

Create the object for reading the data

'''
ReadObj = IAGA2002()
'''

Tell the object which file to read

'''
ReadObj.read('ott20130101vmin.min')
'''

If there more then 1 file to read, create another object and append it to previous

'''
ReadObj2 = IAGA2002()
ReadObj2.read('ott20130102vmin.min')
ReadObj.append(ReadObj2)
'''

You can now get the time and the components X,Y,Z,F,G,H,D,I from the object

'''
time = ReadObj.get(ReadObj.datetime_index)
x = ReadObj.get('x')
'''

The components are numpy arrays and the time is an array of datetime objects.

If the input file contains other components that is being requested, these are calculated.  For example, if the input file contains XYZF and we want H component, this component will be calculated using the proper conversion algorithm as listed at http://www.geomag.nrcan.gc.ca/mag_fld/comp-eng.php.

The data can also be converted to another format using the cast function

'''
WriteObj = ReadObj.cast(IMF)
WriteObj.write(output_file)
'''

The above example will transform are IAGA2002 file to IMFv122 format.

In addition, information stored in the file can also be acquired using any of the following routines.  Some of these function may return empty responses (None) if the information was not stored in the input file.

- getIAGA
- getStationName
- getIntitute
- getGIN
- getReportedOrientation
- getSensorOrientation
- getSamplingRate
- getDigitalSampling
- getDataIntervalType
- getDataType
- getLongitude
- getLatitude
- getAltitude
- getCoordinates
-- returns tuple with (longitude,latitude,altitude)
- getDecbas


