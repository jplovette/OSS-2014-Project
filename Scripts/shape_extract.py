#shape_extract
#author: Tony Chang
#abstract: This script takes the HUC shapefiles and extracts underlying rasters for analysis
#date: 08.05.2014
#organization: OSS 2014 
#location: RENCI, Chapel Hill, NC
#written for Python 3.3

import numpy as np
from matplotlib import pyplot as plt
import gdal, osr
import ogr
import pandas
from PIL import Image, ImageDraw

def getFeatureExtent(pnts):
	xmin = np.min(pnts[0])
	ymin = np.min(pnts[1])
	xmax = np.max(pnts[0])
	ymax = np.max(pnts[1])
	out = [xmin,ymin,xmax,ymax]
	#extent = source_layer.GetExtent() #use this to reference which tile will be used to gather data from
	#print(fextent) 
	return(out)
	
def findSource(extent,var):
	logfile = 'C://Users/tony/OSS/Git/OSS-2014-Project/Scripts/tile_extents_%s.csv' %(var) 
	log = pandas.io.parsers.read_csv(logfile)
	sxmin = extent[0]
	symin = extent[1]
	sxmax = extent[2]
	symax = extent[3]
	#first pass check upper left (UL) corner
	filter_ul = log[(log[' minX']<sxmin) & (log[' maxX']>sxmin) & (log[' minY']<symin) & (log[' maxY']>symin)]
	filter_lr = log[(log[' minX']<sxmax) & (log[' maxX']>sxmax) & (log[' minY']<symax) & (log[' maxY']>symax)]
	#if len(filter_one) = 0, then we have a shape that is between tiles..
	f_ul = str(filter_ul['filename'].iloc[0])
	f_lr = str(filter_lr['filename'].iloc[0])
	if (f_ul == f_lr):
		fname_list = np.array([f_ul])
	else:
		fname_list = np.array([f_ul, f_lr])
	return(fname_list)

def getPoints(poly):
	poly.DumpReadable()
	# dump the data, just to check the first time we look at it
	geom = poly.GetGeometryRef()
	pts = geom.GetGeometryRef(0)
	points = []
	for p in range(pts.GetPointCount()):
		points.append((pts.GetX(p), pts.GetY(p))) #this picks up all the vertices for the n-th feature
	pnts = np.array(points).transpose() #make it a column array
	#print(points[:10]) #print the first ten points to check
	return(pnts)

def world2Pixel(geoMatrix, x, y):
	"""
	Uses a gdal geomatrix (gdal.GetGeoTransform()) to calculate
	the pixel location of a geospatial coordinate
	"""
	ulX = geoMatrix[0]
	ulY = geoMatrix[3]
	xDist = geoMatrix[1]
	yDist = geoMatrix[5]
	rtnX = geoMatrix[2]
	rtnY = geoMatrix[4]
	pixel = np.round((x - ulX) / xDist).astype(np.int)
	line = np.round((ulY - y) / xDist).astype(np.int)
	return (pixel, line)

def imageToArray(i):
# This function will convert the rasterized clipper shapefile
# to a mask for use within GDAL.
    """
    Converts a Python Imaging Library array to a
    numpy array.
    """
    a=np.fromstring(i.tostring(),'b')
    a.shape=i.im.size[1], i.im.size[0]
    return(a)

def transformProj(swkt, rwkt, poly):
	oSRS = osr.SpatialReference ()
	oSRSop = osr.SpatialReference()
	oSRSop.ImportFromWkt(rwkt)
	# wkt from above, is the wicket from the shapefile
	oSRS.ImportFromWkt(swkt)
	# now make sure we have the shapefile geom
	geom = poly.GetGeometryRef()
	pts = geom.GetGeometryRef(0)
	# pts is the polygon of interest
	pts.AssignSpatialReference(oSRS)
	# so transform it to the MODIS geometry
	pts.TransformTo(oSRSop)
	# extract and plot the transformed data
	points = []
	for p in range(pts.GetPointCount()):
		points.append((pts.GetX(p), pts.GetY(p)))
	pnts = np.array(points).transpose()
	return(pnts)

def rasterizer(pnts,raster):
	geo_t = raster.GetGeoTransform()
	pixel, line = world2Pixel(geo_t,pnts[0],pnts[1])
	rasterPoly = Image.new("L", (raster.RasterXSize, raster.RasterYSize),1)
	rasterize = ImageDraw.Draw(rasterPoly)
	listdata = [(pixel[i],line[i]) for i in range(len(pixel))]
	rasterize.polygon(listdata,0)
	mask = 1 - imageToArray(rasterPoly)
	return(mask)

def subsetByMask(mask, fullset):
	#subset the raster 
	#cut the raster to a smaller extent to perform mask analysis...
	aoa_region = np.where(mask == 1)
	minxi = np.min(aoa_region[0])
	minyi = np.min(aoa_region[1])
	maxxi = np.max(aoa_region[0])
	maxyi = np.max(aoa_region[1])
	return(fullset[minxi:maxxi,minyi:maxyi])
		
def forestLossByYear(year_loss, tc_values):
	pixelsize = 30*30
	lossByYearArray = []
	nyears = 13
	min_tc = 90 #minimum amount to consider as forested
	total_tc = len(tc_values[tc_values>=min_tc])
	for i in range(1,nyears):
		losses = (tc_values[year_loss==i])
		#filter only areas where tree cover is greater than 90%
		tc_losses = losses[losses>=min_tc]
		lossByYearArray.append(len(tc_losses)/total_tc)
	return(lossByYearArray)

def culmForestLoss(flby):
	culm_array = []
	for i in range(len(flby)):
		culm_array.append(np.sum(flby[:i]))
	return(np.array(culm_array))
	
#===================================================================================================
#===================================================================================================
#===================================================================================================

if __name__ == "__main__":
	shapeworkspace = 'C://Users/tony/OSS/Project/data/Basin_ref-2014-08-05/Basin_ref/'
	shapefile = 'CONUS_base_ref.shp' #using the huc8 for test file
	sname = '%s%s' %(shapeworkspace,shapefile)
	orig_data_source = ogr.Open(sname)
	source_ds = ogr.GetDriverByName("Memory").CopyDataSource(orig_data_source, "")
	source_layer = source_ds.GetLayer(0)
	source_srs = source_layer.GetSpatialRef()
	wkt = source_srs.ExportToWkt() #this is the full definition of the projection as a string
	total_n = source_layer.GetFeatureCount()
	watershed_summary = []
	for n in range(total_n):
		poly = source_layer.GetFeature(n) #look up record n (based on the id)
		g_id = poly.GetField("GAGE_ID")
		pnts = getPoints(poly)
		fextent = getFeatureExtent(pnts)
		var = 'lossyear' #<======== run this again for another variable (treecover2000)
		fname_lossyear = findSource(fextent,var) #look for the correct raster file to use
		if (len(fname_lossyear) >1):
			continue
			#skip boundary watershed for now.....
			
		#======Open up the raster...
		raster = gdal.Open(fname_lossyear[0])
		fcWKT = raster.GetProjectionRef() # get the raster wicket (the projection information)

		#======This step can be ignored since the shapefile and raster are in the same projection..
		pnts_T = transformProj(wkt, fcWKT, poly) #transforms the shape to match the raster
		mask = rasterizer(pnts, raster) #rasterize the points

		r_array = raster.ReadAsArray()
		raster_sub = subsetByMask(mask,r_array)
		mask_sub = subsetByMask(mask,mask)
		
		#now pull information for that watershed
		#find where mask equals 1 and pull information from the raster
		mask_i = np.where(mask_sub == 1)
		year_loss = raster_sub[mask_i]
		
		#now look at tree cover 2000
		var = 'treecover2000'
		fname_tcover = findSource(fextent,var)
		tc_raster = (gdal.Open(fname_tcover[0])).ReadAsArray()
		tc_sub = subsetByMask(mask,tc_raster)
		tc_values = tc_sub[mask_i]
		
		flby = forestLossByYear(year_loss,tc_values)
		culm_loss = culmForestLoss(flby)
		watershed_summary.append([g_id,flby,culm_loss])
	