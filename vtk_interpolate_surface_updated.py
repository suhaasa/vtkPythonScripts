import sys
import vtk
import numpy as np
from vtk.util import numpy_support as VN
import argparse

def calcDistance2Points(model, pt1,pt2):
	x1,y1,z1 = model.GetPoint(pt1)
	x2,y2,z2 = model.GetPoint(pt2)
	distance = ((x1-x2)**2 + (y1-y2)**2 + (z1-z2)**2)**(.5)
	return distance

def maxDistanceBetweenPoints(model, seedPt, connectedPts_list):
	max = 0
	for i in xrange(0,connectedPts_list.GetNumberOfIds()):
		distance = calcDistance2Points(model, seedPt,connectedPts_list.GetId(i))
		if(distance > max):
			max = distance
	return max

def minDistanceBetweenPoints(model, seedPt, connectedPts_list):
	min = 100000
	min_pt_index = 0
	for iPt in xrange(0,connectedPts_list.GetNumberOfIds()):
		distance = calcDistance2Points(model, seedPt,connectedPts_list.GetId(iPt))
		if(distance < min):
			min = distance
			min_pt_index = iPt
	return min,min_pt_index

def minDistanceBetweenPointsinSet(model, seedPt, connectedPts_list):
	min = 100000
	for iPt in connectedPts_list:
		distance = calcDistance2Points(model, seedPt,iPt)
		if(distance < min):
			min = distance
	return min

def calcDistanceAlongSurface(mesh,startPt,endPt):
	dijkstra = vtk.vtkDijkstraGraphGeodesicPath()
	dijkstra.SetInput(mesh)
	dijkstra.SetStartVertex(startPt)
	dijkstra.SetEndVertex(endPt)
	dijkstra.Update()
	pts = dijkstra.GetOutput().GetPoints()
	dist = 0.0
	for ptId in range(pts.GetNumberOfPoints()-1):
		pts.GetPoint(ptId, p0)
		pts.GetPoint(ptId+1, p1)
		dist += math.sqrt(vtk.vtkMath.Distance2BetweenPoints(p0, p1))
	return dist


def getConnectedVerticesWithinRadius(model, seedPt, radius):
	connectedPts_list = vtk.vtkIdList()
	x,y,z = model.GetPoint(seedPt)
	sphere = vtk.vtkSphere()
	sphere.SetRadius(radius)
	sphere.SetCenter(x,y,z)
	clipper = vtk.vtkClipPolyData()
	clipper.SetClipFunction(sphere)
	clipper.SetInputData(model)
	clipper.GenerateClippedOutputOn()
	clipper.Update()
	# slice_writer = vtk.vtkXMLPolyDataWriter()
	# slice_writer.SetInputData(clipper.GetClippedOutput())
	# slice_writer.SetFileName('Region_cut' + '.vtp')
	# slice_writer.Write()
	connect = vtk.vtkConnectivityFilter()
	connect.SetInputData(clipper.GetClippedOutput())
	connect.SetExtractionModeToClosestPointRegion()
	connect.SetClosestPoint(x,y,z)
	connect.Update()
	region = connect.GetOutput()
	for ptID in xrange(0,region.GetNumberOfPoints()):
		xp,yp,zp = region.GetPoint(ptID)
		modelpt = model.FindPoint(xp,yp,zp)
		connectedPts_list.InsertUniqueId(modelpt)
	return connectedPts_list

def smoothPointSet(mesh,data,seedPt,radius,p):
	ptList = getConnectedVerticesWithinRadius(mesh,seedPt,radius)
	weights = []
	sum_distances = 0
	for i in xrange(0,ptList.GetNumberOfIds()):
		distance = calcDistance2Points(mesh,seedPt,ptList.GetId(i))
		value = getValue(mesh,ptList.GetId(i),data)
		if(distance==0):
			weighting = 1
		else:
			weighting = 1/distance**p*value
			sum_distances += 1/distance**p
		weights.append(weighting)
		
	interpolated_value = 0
	if(sum_distances==0):
		interpolated_value = 0
	else:
		for i in xrange(0,len(weights)):
			interpolated_value += weights[i]/sum_distances
	return interpolated_value

def getConnectedVerticesNotIncludingSeed(model, seedPt):
	cell_list = vtk.vtkIdList()
	connectedPts_list = vtk.vtkIdList()
	model.GetPointCells(seedPt,cell_list)
	for j in xrange(0,cell_list.GetNumberOfIds()):
		pt_list = vtk.vtkIdList()
		pt_list = model.GetCell(cell_list.GetId(j)).GetPointIds()
		for k in xrange(0,pt_list.GetNumberOfIds()):
			if (pt_list.GetId(k) != seedPt):
				connectedPts_list.InsertUniqueId(pt_list.GetId(k))
	return connectedPts_list

def interpolatePointSet(mesh,data,seedPt,ptList,p,maxRadius):
	weights = []
	sum_distances = 0
	for i in xrange(0,ptList.GetNumberOfIds()):
		distance = calcDistance2Points(mesh,seedPt,ptList.GetId(i))
		value = getValue(mesh,ptList.GetId(i),data)
		if(distance==0):
			weighting = 1
		elif(distance < maxRadius):
			weighting = 1/distance**p*value
			sum_distances += 1/distance**p
			weights.append(weighting)
	interpolated_value = 0
	if(sum_distances==0):
		interpolated_value = 0
	else:
		for i in xrange(0,len(weights)):
			interpolated_value += weights[i]/sum_distances
	return interpolated_value

def interpolatePointSetShepardsMethod(mesh,data,seedPt,ptList,p,maxRadius):
	weights = []
	sum_distances = 0
	for i in xrange(0,ptList.GetNumberOfIds()):
		distance = calcDistance2Points(mesh,seedPt,ptList.GetId(i))
		distance_along_mesh = calcDistanceAlongSurface(mesh,seedPt,ptList.GetId(i))
		print(distance)
		print(distance_along_mesh)
		value = getValue(mesh,ptList.GetId(i),data)
		if(distance==0):
			weighting = 1
		elif(distance < 1/3*maxRadius):
			p=2
			weighting = 1/distance**p*value
			sum_distances += 1/distance**p
			weights.append(weighting)
		else:
			p=4
			weighting = 1/distance**p*value
			sum_distances += 1/distance**p
			weights.append(weighting)
	interpolated_value = 0
	if(sum_distances==0):
		interpolated_value = 0
	else:
		for i in xrange(0,len(weights)):
			interpolated_value += weights[i]/sum_distances
	return interpolated_value

def interpolatePointSetModifiedShepardsMethod(mesh,data,seedPt,ptList,p,maxRadius):
	weights = []
	sum_distances = 0
	for i in xrange(0,ptList.GetNumberOfIds()):
		distance = calcDistance2Points(mesh,seedPt,ptList.GetId(i))
		#distance_along_mesh = calcDistanceAlongSurface(mesh,seedPt,ptList.GetId(i))
		#print(distance)
		#print(distance_along_mesh)
		value = getValue(mesh,ptList.GetId(i),data)
		if(distance==0):
			weighting = 1
		elif(distance > maxRadius):
			weights.append(0)
		else:
			weighting = (maxRadius - distance)/(maxRadius * distance)**p*value
			sum_distances = sum_distances + (maxRadius - distance)/(maxRadius * distance)**p
			weights.append(weighting)
	interpolated_value = 0
	if(sum_distances==0):
		interpolated_value = 0
	else:
		for i in xrange(0,len(weights)):
			interpolated_value += weights[i]/sum_distances
	return interpolated_value

def findLocalMaximas(mesh,data):
	numPts = mesh.GetNumberOfPoints()
	localMaxPtList = vtk.vtkIdList()
	for ptID in xrange(0,numPts):
		pt_list = getConnectedVerticesNotIncludingSeed(mesh,ptID)
		seed_data_value = getValue(mesh,ptID,data)
		max_found = 1
		for i in xrange(0,pt_list.GetNumberOfIds()):
			data_value = getValue(mesh,pt_list.GetId(i),data)
			if(data_value>=seed_data_value):
				max_found = 0
		if max_found == 1:
			localMaxPtList.InsertUniqueId(ptID)
	return localMaxPtList

def findNuclearLocalMaximas(mesh,data,nuclear_data,Threshold):
	numPts = mesh.GetNumberOfPoints()
	localMaxPtList = vtk.vtkIdList()
	for ptID in xrange(0,numPts):
		pt_list = getConnectedVerticesNotIncludingSeed(mesh,ptID)
		seed_data_value = getValue(mesh,ptID,data)

		max_found = 1
		nuclear_found = 1
		nuclear_value = 0
		for i in xrange(0,pt_list.GetNumberOfIds()):
			data_value = getValue(mesh,pt_list.GetId(i),data)
			nuclear_value = nuclear_value + getValue(mesh,pt_list.GetId(i),nuclear_data)
			if(data_value>=seed_data_value):
				max_found = 0
		if(pt_list.GetNumberOfIds()<=0):
			nuclear_value = 0
		else:
			nuclear_value = nuclear_value / pt_list.GetNumberOfIds()
		if(nuclear_value < Threshold):
			nuclear_found = 0
		if max_found == 1 and nuclear_found ==1:
			localMaxPtList.InsertUniqueId(ptID)
	return localMaxPtList

def findDistancelocalMaximas(mesh):
	numPts = mesh.GetNumberOfPoints()
	localMaxPtList = vtk.vtkIdList()
	for ptID in xrange(0,numPts):
		pt_list = getConnectedVerticesNotIncludingSeed(mesh,ptID)
		seed_data_value = getValue(mesh,ptID,'Nuclei_distances')
		max_found = 1
		for i in xrange(0,pt_list.GetNumberOfIds()):
			data_value = getValue(mesh,pt_list.GetId(i),'Nuclei_distances')
			if(data_value>=seed_data_value):
				max_found = 0
		if max_found == 1:
			localMaxPtList.InsertUniqueId(ptID)
	return localMaxPtList

def findNuclearPoints(mesh,nuclear_data,Threshold):
	numPts = mesh.GetNumberOfPoints()
	NuclearPtList = vtk.vtkIdList()
	for ptID in xrange(0,numPts):
		seed_nuclear_value = getValue(mesh,ptID,nuclear_data)
		if seed_nuclear_value  > Threshold:
			NuclearPtList.InsertUniqueId(ptID)
	return NuclearPtList

def averagePointList(mesh,pt_list,data):
	numPts = mesh.GetNumberOfPoints()
	data_array = vtk.vtkDoubleArray()
	for ptID in xrange(0,numPts):
		value = getValue(mesh,ptID,data)
		if (pt_list.IsId(ptID) == -1):
			data_array.InsertNextValue(value)
		else:
			connectedPt_list = getConnectedVerticesNotIncludingSeed(mesh, ptID)
			sum = 0
			for i in xrange(0,connectedPt_list.GetNumberOfIds()):
				sum = sum + getValue(mesh,pt_list.GetId(i),data)
			average = sum / connectedPt_list.GetNumberOfIds()
			data_array.InsertNextValue(average)
	return data_array

def averagePointListWithinRadius(mesh,data,radius):
	numPts = mesh.GetNumberOfPoints()
	data_array = vtk.vtkDoubleArray()
	for ptID in xrange(0,numPts):
		update_progress(ptID,numPts)
		value = getValue(mesh,ptID,data)
		connectedPt_list = getConnectedVerticesWithinRadius(mesh, ptID, radius)
		sum = 0
		for i in xrange(0,connectedPt_list.GetNumberOfIds()):
			sum = sum + getValue(mesh,connectedPt_list.GetId(i),data)
		average = sum / connectedPt_list.GetNumberOfIds()
		data_array.InsertNextValue(average)
	return data_array

def kNN(mesh,data,k,controlPt_list):
	numPts = mesh.GetNumberOfPoints()
	data_array = vtk.vtkDoubleArray()
	#loop through each point in the mesh
	for ptID in xrange(0, numPts):
		update_progress(ptID,numPts)
		#search until k nearest neighbors are found
		found_NN = 0
		edgePt_list = vtk.vtkIdList()
		patchPt_list = vtk.vtkIdList()
		NNPt_list = vtk.vtkIdList()
		edgePt_list.InsertUniqueId(ptID)
		patchPt_list.InsertUniqueId(ptID)
		while found_NN < k :
			found_NN = 0
			temp_list = vtk.vtkIdList()
			numIds = edgePt_list.GetNumberOfIds()
			for i in  xrange(0,numIds):
				connnectedPt_list = getConnectedVerticesNotIncludingSeed(mesh,edgePt_list.GetId(i))
				for j in xrange(0,connnectedPt_list.GetNumberOfIds()):
					if(patchPt_list.IsId(connnectedPt_list.GetId(j)) == -1):
						temp_list.InsertUniqueId(connnectedPt_list.GetId(j))
						patchPt_list.InsertUniqueId(connnectedPt_list.GetId(j))
			edgePt_list.Reset()
			for i in xrange(0,temp_list.GetNumberOfIds()):
				edgePt_list.InsertUniqueId(temp_list.GetId(i))
			for i in xrange(0,patchPt_list.GetNumberOfIds()):
				if(controlPt_list.IsId(patchPt_list.GetId(i)) > -1):
					found_NN = found_NN + 1
					NNPt_list.InsertUniqueId(patchPt_list.GetId(i))
			if(patchPt_list.GetNumberOfIds()>100):
				break
		print(patchPt_list.GetNumberOfIds())
		#average those k nearest neighbors using the array "data"
		sum = 0
		for i in xrange(0,NNPt_list.GetNumberOfIds()):
			sum = sum + getValue(mesh,NNPt_list.GetId(i),data)
		if(NNPt_list.GetNumberOfIds()==0):
			data_array.InsertNextValue(0)
		else:
			average = sum / NNPt_list.GetNumberOfIds()
			#set the interpolated value of that point
			data_array.InsertNextValue(average)
	return data_array

def patchAverage(mesh,data,maxIds):
	numPts = mesh.GetNumberOfPoints()
	data_array = vtk.vtkDoubleArray()
	#loop through each point in the mesh
	for ptID in xrange(0, numPts):
		update_progress(ptID,numPts)
		#search until k nearest neighbors are found
		edgePt_list = vtk.vtkIdList()
		patchPt_list = vtk.vtkIdList()
		edgePt_list.InsertNextId(ptID)
		patchPt_list.InsertNextId(ptID)
		while patchPt_list.GetNumberOfIds() < maxIds:
			temp_list = vtk.vtkIdList()
			numIds = edgePt_list.GetNumberOfIds()
			for i in  xrange(0,numIds):
				connnectedPt_list = getConnectedVerticesNotIncludingSeed(mesh,edgePt_list.GetId(i))
				for j in xrange(0,connnectedPt_list.GetNumberOfIds()):
					if(patchPt_list.IsId(connnectedPt_list.GetId(j)) == -1):
						temp_list.InsertNextId(connnectedPt_list.GetId(j))
						patchPt_list.InsertUniqueId(connnectedPt_list.GetId(j))
			edgePt_list.Reset()
			for i in xrange(0,temp_list.GetNumberOfIds()):
				edgePt_list.InsertNextId(temp_list.GetId(i))
		#average those k nearest neighbors using the array "data"
		sum = 0
		for i in xrange(0,patchPt_list.GetNumberOfIds()):
			sum = sum + getValue(mesh,patchPt_list.GetId(i),data)
		average = sum / patchPt_list.GetNumberOfIds()
		#set the interpolated value of that point
		data_array.InsertNextValue(average)
	return data_array

def findNuclei(mesh,ptList,nuclear_data,Threshold,threshold_gradient,vprint):
	nuclei_list = []
	for i in xrange(0,ptList.GetNumberOfIds()):
		update_progress(i,ptList.GetNumberOfIds(),vprint)
		#search until all nuclear points are found
		edgePt_list = vtk.vtkIdList()
		patchPt_list = vtk.vtkIdList()
		edgePt_list.InsertNextId(ptList.GetId(i))
		patchPt_list.InsertNextId(ptList.GetId(i))
		seed_nuc_val = getValue(mesh,ptList.GetId(i),nuclear_data)
		while(edgePt_list.GetNumberOfIds() >= 6 or edgePt_list.GetNumberOfIds() == 1):
			temp_list = vtk.vtkIdList()
			numIds = edgePt_list.GetNumberOfIds()
			for i in  xrange(0,numIds):
				connnectedPt_list = getConnectedVerticesNotIncludingSeed(mesh,edgePt_list.GetId(i))
				for j in xrange(0,connnectedPt_list.GetNumberOfIds()):
					nuc_val = getValue(mesh,connnectedPt_list.GetId(j),nuclear_data)
					if(patchPt_list.IsId(connnectedPt_list.GetId(j)) == -1 and (nuc_val >= Threshold and nuc_val/getValue(mesh,patchPt_list.GetId(0),nuclear_data) >= threshold_gradient)):
						temp_list.InsertNextId(connnectedPt_list.GetId(j))
						patchPt_list.InsertNextId(connnectedPt_list.GetId(j))
			edgePt_list.Reset()
			for i in xrange(0,temp_list.GetNumberOfIds()):
				edgePt_list.InsertNextId(temp_list.GetId(i))
		nuclei_list.append(patchPt_list)

	return nuclei_list

def findNuclei(mesh,ptList,nuclei_data,threshold_gradient,vprint):
	#make intial data structures
	edgePt_list = []
	patchPt_list = []
	usedPts_set = set()
	for i in xrange(0,ptList.GetNumberOfIds()):
		edgePts = vtk.vtkIdList()
		patchPts = vtk.vtkIdList()
		edgePts.InsertNextId(ptList.GetId(i))
		patchPts.InsertNextId(ptList.GetId(i))
		edgePt_list.append(edgePts)
		patchPt_list.append(patchPts)
		usedPts_set.add(ptList.GetId(i))
	count = 0
	while(count<10):
		update_progress(count,10,vprint)
		count = count + 1
		temp_list = []
		for i in xrange(0,ptList.GetNumberOfIds()):
			tempPts = vtk.vtkIdList()
			temp_list.append(tempPts)
		for i in xrange(0,len(edgePt_list)):
			numIds = edgePt_list[i].GetNumberOfIds()
			for j in xrange(0,numIds):
				connectedPts_list = getConnectedVerticesNotIncludingSeed(mesh,edgePt_list[i].GetId(j))
				for k in xrange(0,connectedPts_list.GetNumberOfIds()):
					nuc_val = getValue(mesh,connectedPts_list.GetId(k),nuclei_data)
					if(patchPt_list[i].IsId(connectedPts_list.GetId(j)) == -1 and connectedPts_list.GetId(j) not in usedPts_set and (nuc_val/getValue(mesh,patchPt_list[i].GetId(0),nuclei_data) >= threshold_gradient)):
						temp_list[i].InsertNextId(connectedPts_list.GetId(j))
						patchPt_list[i].InsertNextId(connectedPts_list.GetId(j))
						usedPts_set.add(connectedPts_list.GetId(j))
			edgePt_list[i].Reset()
			for j in xrange(0,temp_list[i].GetNumberOfIds()):
				edgePt_list[i].InsertNextId(temp_list[i].GetId(j))

	return patchPt_list


	#loop through each edgelist in the List
	#add to a used node list

def findNuclei(mesh,nuclei_data,threshold,threshold_gradient,vprint):
	nuclei_pts = set()
	numPts = mesh.GetNumberOfPoints()
	for ptID in xrange(0,numPts):
		if getValue(mesh,ptID,nuclei_data) > threshold:
			nuclei_pts.add(ptID)
	
	all_nuclei_pts = []
	for i in nuclei_pts:
		all_nuclei_pts.append(i)

	nuclei_list = []
	nuclei_border_list = []
	for ptID in all_nuclei_pts:
		if(ptID in nuclei_pts):
			#search until all nuclear points are found
			edgePt_list = vtk.vtkIdList()
			patchPt_list = vtk.vtkIdList()
			edgePt_list.InsertNextId(ptID)
			patchPt_list.InsertNextId(ptID)
			nuclei_border_pts = set()
			while(edgePt_list.GetNumberOfIds() >= 6 or edgePt_list.GetNumberOfIds() == 1):
				temp_list = vtk.vtkIdList()
				numIds = edgePt_list.GetNumberOfIds()
				for i in  xrange(0,numIds):
					connnectedPt_list = getConnectedVerticesNotIncludingSeed(mesh,edgePt_list.GetId(i))
					for j in xrange(0,connnectedPt_list.GetNumberOfIds()):
						nuc_val = getValue(mesh,connnectedPt_list.GetId(j),nuclei_data)
						if(patchPt_list.IsId(connnectedPt_list.GetId(j)) == -1 and connnectedPt_list.GetId(j) in nuclei_pts and nuc_val/getValue(mesh,patchPt_list.GetId(0),nuclei_data) >= threshold_gradient):
							temp_list.InsertNextId(connnectedPt_list.GetId(j))
							patchPt_list.InsertNextId(connnectedPt_list.GetId(j))
							nuclei_pts.remove(connnectedPt_list.GetId(j))
						elif(patchPt_list.IsId(connnectedPt_list.GetId(j)) == -1):
							nuclei_border_pts.add(connnectedPt_list.GetId(j))

				edgePt_list.Reset()
				for i in xrange(0,temp_list.GetNumberOfIds()):
					edgePt_list.InsertNextId(temp_list.GetId(i))
			nuclei_list.append(patchPt_list)
			nuclei_border_list.append(nuclei_border_pts)
	filtered_nuclei_list = []
	filtered_nuclei_border_list = []
	for i in xrange(0,len(nuclei_list)):
		if(nuclei_list[i].GetNumberOfIds() > 10):
			filtered_nuclei_list.append(nuclei_list[i])
			filtered_nuclei_border_list.append(nuclei_border_list[i])

	return filtered_nuclei_list,filtered_nuclei_border_list

def setDistanceValues(mesh,nuclei_list,nuclei_border_set_list):
	numPts = mesh.GetNumberOfPoints()
	nuclei_dist_data = [0]*numPts
	for i in xrange(0,len(nuclei_list)):
		for j in xrange(0,nuclei_list[i].GetNumberOfIds()):
			nuclei_dist_data[nuclei_list[i].GetId(j)] = minDistanceBetweenPointsinSet(mesh,nuclei_list[i].GetId(j),nuclei_border_set_list[i])

	vtk_nuclei_dist_data = vtk.vtkDoubleArray()
	for i in xrange(0,len(nuclei_dist_data)):
		vtk_nuclei_dist_data.InsertNextValue(nuclei_dist_data[i])
	vtk_nuclei_dist_data.SetName('Nuclei_distances')
	mesh.GetPointData().AddArray(vtk_nuclei_dist_data)

def assignNuclearPoints(mesh,nuclei_list,local_max_distance_pt_list):
	numPts = mesh.GetNumberOfPoints()
	nuclei_data = [0]*numPts

	for i in xrange(0,len(nuclei_list)):
		for j in xrange(0,nuclei_list[i].GetNumberOfIds()):
			min, min_pt_index = minDistanceBetweenPoints(mesh,nuclei_list[i].GetId(j),local_max_distance_pt_list)
			nuclei_data[nuclei_list[i].GetId(j)] = min_pt_index

	vtk_nuclei_data = vtk.vtkDoubleArray()
	for i in xrange(0,len(nuclei_data)):
		vtk_nuclei_data.InsertNextValue(nuclei_data[i])
	vtk_nuclei_data.SetName('Filtered_Nuclei')
	mesh.GetPointData().AddArray(vtk_nuclei_data)
			

def intializeVTP(filename):
	datareader=vtk.vtkXMLPolyDataReader()
	datareader.SetFileName(filename)
	datareader.Update()

	mesh=vtk.vtkDataSetMapper()
	mesh=datareader.GetOutput()
	return mesh

def calcMagnitude(vector):
	mag = 0
	for i in xrange(0,len(vector)):
		mag += vector[i]**2
	return mag ** (.5)

def getValue(model, ptID, factor):
	return calcMagnitude(model.GetPointData().GetArray(factor).GetTuple(ptID))

def update_progress(progress, total, vprint):  
    vprint('\r[{0:10}]{1:>2}'.format('#' * int(progress * 10 /total + 1), progress))

def createParser():
	parser = argparse.ArgumentParser(description='Interpolates the given data across mesh surface using Nuclear points as control points.')
	parser.add_argument('filename', type=str, help='the filename (include file ext) to analyze')
	parser.add_argument('-d','--data_to_interpolate', nargs='?', type=str, default='Nuclear Dach1 (Cx40Nuclei)', help='the name of array to interpolate (default=Nuclear Dach1 (Cx40Nuclei))')
	parser.add_argument('-r','--max_radius', type=int, nargs='?', default=20, help='the max radius to interpolate within (default=20)')
	parser.add_argument('-t','--threshold', type=int, nargs='?', default=20, help='the threshold for Nucleus (deafult=20)')
	parser.add_argument('-tg','--threshold_gradient', type=float, nargs='?', default=.80, help='the threshold for the precent Nucleus falloff (deafult=.80)')
	parser.add_argument('-p', type=int, default=16, nargs='?', help='the value of the non-linearity when interpolating')
	parser.add_argument('-vtp', type=str, default=None, help='the name of a vtp output file (without file ext) (def suff=_interpolated)')
	parser.add_argument('-v', '--verbose', type=int, nargs='?', const=1, default=0, help='turn on verbosity')
	parser.add_argument('-a', '--algorithm', type=int, choices=[1,2], nargs='?', default=2, help='1:findNuclearPoints 2:findNuclearLocalMaximas (def: 2)')
	parser.add_argument('-n', '--nucleus', type=str, default='Nucleus (Cx40Nuclei)', help='the name of the Nucleus array')
	parser.add_argument('-calc_max', action='store_true', help='option to turn on the automatic max_radius determination')
	return parser

def main(args):
	filename = args.filename
	data_to_interpolate = args.data_to_interpolate
	nuclear_data = args.nucleus
	Threshold = args.threshold
	threshold_gradient = args.threshold_gradient
	maxRadius = args.max_radius
	p = args.p

	if(args.vtp is None):
		output_filename = args.filename.split('.')[0] + '_interpolated.vtp'
	else:
		output_filename = args.vtp + '.vtp'

	if args.verbose:
		def vprint(*args):
			# Print each argument separately so caller doesn't need to
			# stuff everything to be printed into a single string
			for arg in args:
				print(arg),
	else:
		vprint = lambda *a: None      # do-nothing function


	mesh = intializeVTP(filename)
	interpolated_value = vtk.vtkDoubleArray()
	control_pts = vtk.vtkDoubleArray()
	numPts = mesh.GetNumberOfPoints()
	
	if(args.algorithm==1):
		ptList = findNuclearPoints(mesh,nuclear_data,Threshold)
	else:
		ptList = findNuclearLocalMaximas(mesh,data_to_interpolate,nuclear_data,Threshold)

	nuclei_list,nuclei_border_list = findNuclei(mesh,nuclear_data,Threshold,threshold_gradient,vprint)
	nuclei_data = [0]*numPts
	for i in xrange(0,len(nuclei_list)):
		for j in xrange(0,nuclei_list[i].GetNumberOfIds()):
			nuclei_data[nuclei_list[i].GetId(j)] = i
	vtk_nuclei_data = vtk.vtkDoubleArray()
	for i in xrange(0,len(nuclei_data)):
		vtk_nuclei_data.InsertNextValue(nuclei_data[i])
	vtk_nuclei_data.SetName('Nuclei')
	mesh.GetPointData().AddArray(vtk_nuclei_data)
	vprint('Nuclei determined.')

	# assigns the distance of each nuclear point to its border
	setDistanceValues(mesh,nuclei_list,nuclei_border_list)

	# finds the local maximas of the nucleus distance array
	local_max_distance_pt_list = findDistancelocalMaximas(mesh)

	# assign nuclear points to a local maxima
	assignNuclearPoints(mesh,nuclei_list,local_max_distance_pt_list)

	w = vtk.vtkXMLPolyDataWriter()
	w.SetInputData(mesh)
	w.SetFileName(output_filename)
	w.Write()

	if(args.calc_max):
		maxRadius = 0 
		for i in xrange(0,ptList.GetNumberOfIds()):
			max_val = maxDistanceBetweenPoints(mesh, ptList.GetId(i), ptList)
			if max_val > maxRadius:
				maxRadius = max_val
		maxRadius = maxRadius / 2
		vprint('Radius= ', maxRadius)

	for ptID in xrange(0,numPts):
		update_progress(ptID,numPts,vprint)
		if(ptList.IsId(ptID) == -1):
			interpolated_value.InsertNextValue(interpolatePointSetModifiedShepardsMethod(mesh,data_to_interpolate,ptID,ptList,p,maxRadius))
			control_pts.InsertNextValue(0)
		else:
			interpolated_value.InsertNextValue(getValue(mesh,ptID,data_to_interpolate))
			control_pts.InsertNextValue(1)

	# TO IMPLEMENT: The seed point nuclei averaging

	interpolated_value.SetName('Interpolated ' + data_to_interpolate)
	mesh.GetPointData().AddArray(interpolated_value)

	control_pts.SetName('Control Points')
	mesh.GetPointData().AddArray(control_pts)

	w = vtk.vtkXMLPolyDataWriter()
	w.SetInputData(mesh)
	w.SetFileName(output_filename)
	w.Write()

if __name__ == '__main__':
	parser = createParser()
	args = parser.parse_args()
	main(args)
