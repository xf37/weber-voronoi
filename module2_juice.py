#-------------------------------------------------------------------------------
# Name:        module2
# Purpose:
#
# Author:      Xin Feng
#
# Created:     11/12/2019
# Copyright:   (c) Xin Feng 2019
# Licence:     <your licence>
#-------------------------------------------------------------------------------

import arcpy
import numpy
import math
import csv
import time
from sympy.geometry import *
import csv
from datetime import datetime
from dbfread import DBF
import random


##########################################################################
# setting variables
##########################################################################
# number of points
n_pts = 1180
# number of markets
n_mkts = 4
# number of resource types plus one
n_rsc = 9
# number of test
n_test = 0
# tolerence value
diff_value = 0.0000001
# iteration number
diff_inter = 500



##########################################################################
# define path
##########################################################################
# dictionary for imput data
path0 = "W:\\data\\upload\\application 2\\"
# dictionary for output data
path = "W:\\data\\upload\\test\\"



##########################################################################
# define others
##########################################################################
# add coordinates
field_x = "POINT_X"
field_y = "POINT_Y"
boundary = path0 + "boundary.shp"
inFeatures = []
Extent = "6246825.05248557 1998162.29366227 6441995.37663836 2297729.15019917"
resouce_vor =  "voronoi_RM"
market = "markets"
vor_intersect = "voronoi_intersect"
vor_points = "vor_its_points"


##########################################################################
# Process: Create Thiessen Polygons based on layers of resouce points (except markets)
##########################################################################
start2 = time.time()
for i in range(1, n_rsc):
    points_shp = str(i)
    voronoi_shp = path + resouce_vor + str(i) + ".shp"
    tempEnvironment0 = arcpy.env.extent
    arcpy.env.extent = Extent
    arcpy.CreateThiessenPolygons_analysis(path0 + points_shp + ".shp", voronoi_shp , "ONLY_FID")
    arcpy.env.extent = tempEnvironment0
    inFeatures.append(voronoi_shp)
end2 = time.time() - start2
print end2



##########################################################################
# Process: Intersect
##########################################################################
start3 = time.time()
inFeatures.append(boundary)
multi_voronoi = path + vor_intersect + ".shp"
arcpy.Intersect_analysis(inFeatures, multi_voronoi, "ALL", "", "INPUT")
end3 = time.time() - start3
print end3



##########################################################################
# Process: Feature To Point
###########################################################################
multi_vpoints = path + vor_points + ".shp"
arcpy.FeatureToPoint_management(multi_voronoi, multi_vpoints, "INSIDE")


##########################################################################
# creat fields and calculate coordinates
##########################################################################
# for voronoi polygons' centroids   (as the start point for heuristic)
arcpy.AddField_management(multi_vpoints, field_x, "DOUBLE")
arcpy.AddField_management(multi_vpoints, field_y, "DOUBLE")
arcpy.CalculateField_management(multi_vpoints, field_x, "!SHAPE.CENTROID.X!", "PYTHON_9.3")
arcpy.CalculateField_management(multi_vpoints, field_y, "!SHAPE.CENTROID.Y!", "PYTHON_9.3")



##########################################################################
#get coordinate & demand(weight) for RM & market
##########################################################################
start0 = time.time()
coo = {}   # coordinate for RM
for i in range(1, n_rsc):
     coo[i] = []
     table = DBF(path0 + str(i) + '.dbf', load = True)
     length = len(table)
     for j in range(length):
         point = []
         point.append(table.records[j][field_x])
         point.append(table.records[j][field_y])
         point.append(1)  # demand weight, may change later
         coo[i].append(point)

coo_m = []   # coordinate for market
table = DBF(path0 + market + '.dbf', load = True)
length_m = len(table)
for j in range(length_m):
     point = []
     point.append(table.records[j][field_x])
     point.append(table.records[j][field_y])
     point.append(1)    # market weight, may change later
     coo_m.append(point)



##########################################################################
#get info from multi-voronoi-points layer
##########################################################################
id_v = []
coo_vpts = []
table_v = DBF(path + vor_points +  ".dbf", load = True)

length_v = len(table_v)
for k in range(length_v):
    v_polygon = []
    for j in range(1, n_rsc):
        j_1 = j - 1
        if j_1 == 0:
            id_j_1 = "Input_FID"
            v_polygon.append(table_v.records[k][id_j_1])
        elif j_1 < 10:
            id_j_1 = "Input_FI_" + str(j_1)
            v_polygon.append(table_v.records[k][id_j_1])
        elif j_1 >= 10:
            id_j_1 = "Input_F_" + str(j_1)
            v_polygon.append(table_v.records[k][id_j_1])
    id_v.append(v_polygon)
    v_point = []
    v_point.append(table_v.records[k][field_x])
    v_point.append(table_v.records[k][field_y])
    coo_vpts.append(v_point)

end0 = time.time() - start0
print end0


##########################################################################
#Function: get Coordinate based on ID of given RM & market
##########################################################################
def getCoordinate(set_rm):
    X_l = []
    Y_l = []
    w_l = []
    for i_mkt in range(n_mkts):
        X_l.append(coo_m[i_mkt][0])
        Y_l.append(coo_m[i_mkt][1])
        w_l.append(coo_m[i_mkt][2])
    for i_rm in range(1,n_rsc):
        id_rm = set_rm[i_rm - 1]
        X_l.append(coo[i_rm][id_rm][0])
        Y_l.append(coo[i_rm][id_rm][1])
        w_l.append(coo[i_rm][id_rm][2])
    return X_l,Y_l,w_l


##########################################################################
#Function: get weber point based on given RM & market
##########################################################################
def getWeber(xc,yc,demand_X_l,demand_Y_l,demand_a_l):
    xk_n = 0
    xk_den = 0
    yk_n = 0
    yk_den = 0
    total_w_d = 0
    for i_xk in range(len(demand_X_l)):
        aixi = demand_a_l[i_xk] * demand_X_l[i_xk]
        xxiyyi = math.sqrt(pow((xc - demand_X_l[i_xk]),2) + pow((yc - demand_Y_l[i_xk]),2))+0.0000000000001
        xk_n = xk_n + aixi / xxiyyi
        xk_den = xk_den + demand_a_l[i_xk] / xxiyyi
        aiyi = demand_a_l[i_xk] * demand_Y_l[i_xk]
        yk_n = yk_n + aiyi / xxiyyi
        yk_den = yk_den + demand_a_l[i_xk] / xxiyyi
        weighted_dist = demand_a_l[i_xk] * (math.sqrt(pow((xc - demand_X_l[i_xk]),2) + pow((yc - demand_Y_l[i_xk]),2)))
        total_w_d = total_w_d + weighted_dist
##    print 'total_w_d is :' + repr(total_w_d)
    xk = xk_n / xk_den
    yk = yk_n / yk_den
    xyk_l = []
    xyk_l.append(xk)
    xyk_l.append(yk)
    return xyk_l, total_w_d


##########################################################################
#Function: get TWD
##########################################################################
def getTWD(xc,yc,demand_X_l,demand_Y_l,demand_a_l):
    total_w_d = 0
    for i_xk in range(len(demand_X_l)):
        weighted_dist = demand_a_l[i_xk] * (math.sqrt(pow((xc - demand_X_l[i_xk]),2) + pow((yc - demand_Y_l[i_xk]),2)))
        total_w_d = total_w_d + weighted_dist
##    print 'total_w_d is :' + repr(total_w_d)
    return total_w_d




##########################################################################
#Get Weber point based on each voronoi polygon
##########################################################################
start1 = time.time()
#define a list to record weighted distance
wd_l = []    # weighted distance list
coo_l = []    # coordinate list
twd_l = []    # total weighted distance list

for i_v in range(length_v):
    wd_l.append([])
    coo_l.append([])
    twd_l.append([])

for row in arcpy.da.SearchCursor(multi_voronoi, ["FID", "SHAPE@"]):
    print row[0]
    [X_l, Y_l, w_l] = getCoordinate(id_v[row[0]])
    xco = coo_vpts[row[0]][0]
    yco = coo_vpts[row[0]][1]
    coo_l[row[0]].append(xco)
    coo_l[row[0]].append(yco)
    [xyk_list, total_w_d_new] = getWeber(xco, yco, X_l, Y_l, w_l)
    twd_l[row[0]].append(total_w_d_new)
    diff = math.sqrt(pow((xyk_list[0] - xco),2) + pow((xyk_list[1] - yco),2))
    #print 'original diff for ' + repr(i_v) + ' is:' + repr(diff)
    n_inter = 0

    # Iteration process
    while (diff > diff_value and n_inter < diff_inter):
    #while (n_inter < diff_inter):
        pt = arcpy.Point()
        pt.X = xyk_list[0]
        pt.Y = xyk_list[1]
    # Judge whether Weber point is within the polygon
        if row[1].contains(pt):
            n_inter = n_inter + 1
            [xyk_list_new, total_w_d_new] = getWeber(xyk_list[0],xyk_list[1],X_l, Y_l, w_l)
            coo_l[row[0]].append(pt.X)
            coo_l[row[0]].append(pt.Y)
            twd_l[row[0]].append(total_w_d_new)
            diff = math.sqrt(pow((xyk_list_new[0] - xyk_list[0]),2) + pow((xyk_list_new[1] - xyk_list[1]),2))
            xyk_list = xyk_list_new
        else:
            coo_l[row[0]].append(pt.X)
            coo_l[row[0]].append(pt.Y)
            total_w_d_new = getTWD(pt.X, pt.Y, X_l, Y_l, w_l)
            twd_l[row[0]].append(total_w_d_new)
            #total_w_d_new = 0
            break
    wd_l[row[0]].append(total_w_d_new)

end1 = time.time() - start1
print end1


##########################################################################
# save total weighted distance
##########################################################################
with open(path +'twd.csv', 'wb') as myfile:
    wr = csv.writer(myfile, quoting=csv.QUOTE_ALL)
    for i in range(length_v):
        wr.writerow(twd_l[i])

n_total_iteration = 0
for i in range(len(twd_l)):
    n_total_iteration += len(twd_l[i])

min_each_voronoi = list(map(min, twd_l))
min_twd = min(min_each_voronoi)
min_index =  numpy.argmin(min_each_voronoi)
n_iteration =  len(twd_l[min_index])
n_polygons = len(twd_l)



##########################################################################
# write funtion
##########################################################################
def writeResult(n_pts, n_rsc, n_mkts, n_test, end0, end1, end2, end3, n_iteration, min_twd, min_index, n_polygons):
    result = path + "results" + "_" + str(n_pts) + "_" + str(n_rsc) + "_" + str(n_mkts)  +  ".txt"
    file = open(result, "w+")
    file.write('This is the result for: ' + str(n_pts) + "_" + str(n_rsc) + "_" + str(n_mkts) + "_" + str(n_test) + '\n')
    file.write('data imput time is: ' + str(end0) + 's' + '\n')
    file.write('calculation time is: ' + str(end1) + 's' + '\n')
    file.write('voronoi time is: ' + str(end2) + 's' + '\n')
    file.write('intersect time is: ' + str(end3) + 's' + '\n')
    file.write('the number of iterations for the optimal polygon is: ' + str(n_iteration) + '\n')
    file.write('the total number of iterations is: ' + str(n_total_iteration) + '\n')
    file.write('minimum total cost is:' + str(min_twd) + '\n')
    file.write('the total number of faces is: ' + str(n_polygons) + '\n')
    file.write('the voronoi ID for the optimal solution is: ' + str(min_index) + '\n')
    file.close()

writeResult(n_pts, n_rsc, n_mkts, n_test, end0, end1, end2, end3, n_iteration, min_twd, min_index, n_polygons)




"""
##########################################################################
# record shapefile of weber points
##########################################################################
name = "WeberPts_RM"
fields_cursor = ["SHAPE@XY", "faceID", "iteration", "TWD"]
fc = path + name + ".shp"
cursor = arcpy.da.InsertCursor(fc, fields_cursor)

for nn in range(length_v):
    for i in range(len(coo_l[nn])/2):
        point = (float(coo_l[nn][2*i]),float(coo_l[nn][2*i+1]))
        newRow = [point, nn, i, twd_l[nn][i]]
        cursor.insertRow(newRow)
del cursor
"""

"""
##########################################################################
# load twd list and check minimum total distance & the optimal voronoi ID
##########################################################################
path =  "W:\\data\\weber_simulation\\" + str(100000) + "_" + str(5) + "_" + str(3) + "_" + str(0) + "\\"
with open(path + 'twd.csv', 'rb') as f:
     reader = csv.reader(f)
     twd_l = []
     for subset in reader:
        subset_l = [float(i) for i in subset]
        twd_l.append(subset_l)

min_each_voronoi = list(map(min, twd_l))
min_twd = min(min_each_voronoi)
min_index =  numpy.argmin(min_each_voronoi)
print "minimum total distance is: " + repr(min_twd)
print "the optimal voronoi ID is: " + repr(min_index)
"""





