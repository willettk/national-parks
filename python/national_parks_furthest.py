'''
What point in the continental US is furthest from a National Park:

Option 1: 

    The correct way to do it, I think. 

    Treat the National Parks each as a single point (probably the centroid of the boundary). 
    Create a Voronoi tessellation of the set of points, which splits it into
        polygons that are closer to a single enclosed point than any other. 
    Intersect the Voronoi diagram with the outer border of the US, which
        provides the boundary constraints of the program.
    The point with the greatest distance from any park must lie along one of the
        vertices, either between individual Voronoi cells or the intersection of
        an edge with the outer border. 

    My problems trying to implement this:

        - I can't get the basemap package to install in Python.
        - R has good facilities for this, but I'm having trouble finding an
            outline of the US that's represented as a single polygon; disconnected
            islands and other landmasses result in self-intersections and errors
            in geometry when I try to compute it.
        - CARTO is good for mapping, but the interface changed recently and I don't
            think it's set up for the distance calculations involved.

Option 2:

    The brute force approach.

    Take a rectangular area that bounds the contiguous United States. 
    Create a rectilinear grid over the area and compute the distance to every one of the
        47 National Parks within the area.
    Evaluate the great circle distances to the parks for each point and take the minimum.
    Plot the results on top of a shapefile for the US.

'''

import numpy as np
import shapefile
from matplotlib import pyplot as plt
import matplotlib.path as mplPath

a = 6378137.0 # m
b = 6356752.3 # m

'''

Extreme points of the contiguous (lower 48) United States:

Location                    Latitude    Longitude
---------------------------------------------------
East: Lubec, ME             44.840833,  -67.015556
West: Cape Alava, WA        48.166667,  -124.733333
North: Northwest Angle, MN  49.266667,  -95.05
South: Key West, FL         24.54409,   -81.804905

'''

lon_min = -124.733333
lon_max =  -67.015556
lat_min =   24.54409 
lat_max =   49.266667

def meanRadius(M,N):

    return 2./(1./M + 1./N)

def meridionalRadius(phi,a,b):

    numerator = (a*b)**2
    denominator = ((a*np.cos(phi))**2 + (b*np.sin(phi))**2)**1.5 

    return numerator/denominator

def primeVerticalRadius(phi,a,b):

    numerator = a**2
    denominator = np.sqrt((a*np.cos(phi))**2 + (b*np.sin(phi))**2)

    return numerator/denominator

def centralAngle(phi1,phi2,lambda1,lambda2):

    dlambda = lambda1 - lambda2

    return np.arccos(np.sin(phi1)*np.sin(phi2) + np.cos(phi1)*np.cos(phi2)*np.cos(dlambda))

def arcLength(r,dsigma):

    return r * dsigma

def haversineDistance(phi1,phi2,lambda1,lambda2,r):

    delta_phi = phi1-phi2
    delta_lambda = lambda1-lambda2

    a = (np.sin(delta_phi/2.))**2 + np.cos(phi1)*np.cos(phi2)*(np.sin(delta_lambda/2.))**2
    c = 2 * np.arctan2(np.sqrt(a),np.sqrt(1. - a))

    return r*c

def loadNPSData():

    reader = shapefile.Reader("/Users/willettk/Desktop/nps_boundary/nps_boundary.shp")

    fields = [field[0] for field in reader.fields[1:]]

    npData = []

    excludeCodes = ('WOTR','HAVO','HALE','VIIS',
                    'DENA','GAAR','GLBA','KATM',
                    'KEFJ','KOVA','LACL','WRST',
                    'NPSA')

    npcount = 0
    for sr in reader.iterShapeRecords():
        attr = dict(zip(fields, sr.record))
        if attr['UNIT_TYPE'] == "National Park" and attr['UNIT_CODE'] not in excludeCodes:
            d = dict()
            d['name'] = attr['UNIT_NAME']
            d['code'] = attr['UNIT_CODE']
            d['lonCenter'] = (sr.shape.bbox[0] + sr.shape.bbox[2])/2.
            d['latCenter'] = (sr.shape.bbox[1] + sr.shape.bbox[3])/2.

            npData.append(d)
            npcount += 1

    #print "Total of {} National Parks found in dataset".format(npcount)
    
    return npData

def loadUSPolygon():

    reader = shapefile.Reader("/Users/willettk/Desktop/cb_2015_us_nation_20m/cb_2015_us_nation_20m.shp")
    s = reader.shape()

    # Find the primary polygon that defines the outline of the contiguous US

    s1 = np.array(s.parts[1:])
    s2 = np.array(s.parts[:-1])

    diffarr = s2 - s1

    startval = 0
    maxind = 0
    for i,x in enumerate(diffarr):
        if abs(x) > startval:
            maxind = i
            startval = abs(x)

    usPoints = s.points[s1[maxind-1]:s1[maxind]]
    
    return usPoints

def loadStates():

    reader = shapefile.Reader("/Users/willettk/Desktop/cb_2015_us_state_20m/cb_2015_us_state_20m.shp")

    statePoints = []

    for sr in reader.iterShapeRecords():
        if sr.record[4] not in ('DC','PR','AK','HI'):
            if len(sr.shape.parts) > 1:
                s1 = sr.shape.parts[:-1]
                s2 = sr.shape.parts[1:]
                for i1,i2 in zip(s1,s2):
                    statePoints.append(sr.shape.points[i1:i2])
                statePoints.append(sr.shape.points[sr.shape.parts[-1]:])
            else:
                statePoints.append(sr.shape.points)

    return statePoints

def containsPoint(polygon,point):

    plist = [[x[0],x[1]] for x in polygon]
    p0 = polygon[0]
    plist.append([p0[0],p0[1]])

    bbPath = mplPath.Path(np.array(plist))
    
    result = bbPath.contains_point(point)

    return result
    
def polygonTest():

    usPoints = loadUSPolygon()

    pointTrue = (-105,40)
    pointFalse = (105,40)

    assert containsPoint(usPoints,pointTrue), \
        "Point {} should be in the US lower 48 polygon".format(pointTrue)
    assert np.logical_not(containsPoint(usPoints,pointFalse)), \
        "Point {} should not be in the US lower 48 polygon".format(pointFalse)

def makeGrid(nlat=10,nlon=10):

    latarr = np.linspace(lat_min,lat_max,nlat)
    lonarr = np.linspace(lon_min,lon_max,nlon)

    return latarr,lonarr

def testIowa():

    latIowa = 41.66
    lonIowa = -91.53

    usPoints = loadUSPolygon()
    npData = loadNPSData()
    
    results = []

    # Need to convert these to radians
    
    deg2rad = np.pi / 180.
    
    M = meridionalRadius(lat * deg2rad,a,b)
    N = primeVerticalRadius(lat * deg2rad,a,b)
    r = meanRadius(M,N)
    
    # Check against all National Parks
    
    minDist = 1e7
    minPark = None
    for park in npData:
        '''
        dsigma = centralAngle(lat*deg2rad,park['latCenter']*deg2rad,
                                lon*deg2rad/15.,park['lonCenter']*deg2rad/15.)
        d = arcLength(r,dsigma)
        '''
        d = haversineDistance(latIowa*deg2rad,park['latCenter']*deg2rad,
                                lonIowa*deg2rad,park['lonCenter']*deg2rad,
                                r)
    
        print "Distance of {} ({:7.2f},{:7.2f}) to Iowa ({:7.2f},{:7.2f}) is {:.1f} km.".format(park['code'],park['lonCenter'],park['latCenter'],lonIowa,latIowa,d*1e-3)
        if d < minDist:
            minDist = d
            minPark = park
    
    print "\nNearest park: {} ({:.1f} km).".format(minPark['name'],minDist*1e-3)

if __name__ == "__main__":

    latarr,lonarr = makeGrid(75,75)
    
    usPoints = loadUSPolygon()
    npData = loadNPSData()
    
    pointsOutsidePolygon = 0
    maxDistGlobal = 0.
    maxParkGlobal = None
    maxParkData = None
    maxLocation = ()

    results = []
    for lat in latarr:
        for lon in lonarr:

            # Need to convert these to radians

            deg2rad = np.pi / 180.
    
            M = meridionalRadius(lat * deg2rad,a,b)
            N = primeVerticalRadius(lat * deg2rad,a,b)
            r = meanRadius(M,N)
    
            # Check to see if it's in the US
            if containsPoint(usPoints,(lon,lat)):
                # Check against all National Parks
    
                minDist = 1e7
                minPark = None

                parkData = []

                for park in npData:
                    '''
                    dsigma = centralAngle(lat*deg2rad,park['latCenter']*deg2rad,
                                            lon*deg2rad/15.,park['lonCenter']*deg2rad/15.)
                    d = arcLength(r,dsigma)
                    '''
                    d = haversineDistance(lat*deg2rad,park['latCenter']*deg2rad,
                                            lon*deg2rad,park['lonCenter']*deg2rad,
                                r)

                    if d < minDist:
                        minDist = d
                        minPark = park['name']

                    parkData.append((park,d*1e-3))
    
                #print "Closest park to ({:7.2f},{:7.2f}) is {} at {:.1f} km.".format(lon,lat,minPark,minDist*1e-3)
                if minDist > maxDistGlobal:
                    maxDistGlobal = minDist
                    maxParkGlobal = minPark
                    maxParkData = parkData
                    maxLocation = (lon,lat)

                results.append((lon,lat,minDist*1e-3))

            else:
                pointsOutsidePolygon += 1
                #print "Point {:.3f},{:.3f} is not found in lower 48 US polygon.".format(lon,lat)

    print "{}/{} points in grid were outside polygon.".format(pointsOutsidePolygon,len(latarr)*len(lonarr))
    print "Largest distance within grid is ({:7.2f},{:7.2f}), \n\twhich is {:.1f} km from {}.".format(maxLocation[0],maxLocation[1],maxDistGlobal*1e-3,maxParkGlobal)

    # List the 5 closest National Parks. I think the bottom three should be approximately the same.

    pdList = [(x[0]['name'],x[1]) for x in maxParkData]
    pdListSorted = sorted(pdList,key = lambda x:x[1])
    print "\nThree closest National Parks to the point of inaccessibility:"
    for p in pdListSorted[:3]:
        print "\t{}, {:.0f} km".format(*p)

    # Make a map!

    fig = plt.figure(figsize=(12,8))
    ax = fig.add_subplot(111)

    # Draw the outline of the US

    xarr = [x[0] for x in usPoints]
    yarr = [x[1] for x in usPoints]
    ax.plot(xarr,yarr,color='k')

    # Draw the outlines of the states

    statePoints = loadStates()

    for state in statePoints:
        xState = [x[0] for x in state]
        yState = [x[1] for x in state]
        ax.plot(xState,yState,color='grey')

    lonplot = [x[0] for x in results]
    latplot = [x[1] for x in results]
    dists = [x[2] for x in results]

    plt.scatter(lonplot,latplot,c=dists,marker='h',cmap='viridis_r',edgecolors='none',s=80)
    cbar = plt.colorbar()
    cbar.ax.set_ylabel("distance [km]",fontsize=16)

    # Plot the nadir
    ax.scatter([maxLocation[0]],[maxLocation[1]],s=200,marker='*',c='r')

    # Plot locations of all NPs
    npLats = [x['latCenter'] for x in npData] 
    npLons = [x['lonCenter'] for x in npData] 
    ax.scatter(npLons,npLats,s=20,marker='D',c='brown')

    ax.set_xlabel("Longitude",fontsize=14)
    ax.set_ylabel("Latitude",fontsize=14)
    ax.set_title("Distance to the nearest U.S. National Park",fontsize=20)

    fig.tight_layout()

    plt.show()

