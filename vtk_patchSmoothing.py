import sys
import vtk
import numpy as np
from vtk.util import numpy_support as VN

def calcDistance2Points(model, pt1,pt2):
	x1,y1,z1 = model.GetPoint(pt1)
	x2,y2,z2 = model.GetPoint(pt2)
	distance = ((x1-x2)**2 + (y1-y2)**2 + (z1-z2)**2)**(.5)
	return distance

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

def findNuclearLocalMaximas(mesh,data,nuclear_data,Threshold):
	numPts = mesh.GetNumberOfPoints()
	localMaxPtList = vtk.vtkIdList()
	for ptID in xrange(0,numPts):
		pt_list = getConnectedVerticesNotIncludingSeed(mesh,ptID)
		seed_data_value = getValue(mesh,ptID,data)

		max_found = 1
		nuclear_found = 1
		for i in xrange(0,pt_list.GetNumberOfIds()):
			data_value = getValue(mesh,pt_list.GetId(i),data)
			nuclear_value = getValue(mesh,pt_list.GetId(i),nuclear_data)
			if(data_value>=seed_data_value):
				max_found = 0
			if(nuclear_value < Threshold):
				nuclear_found = 0
		if max_found == 1 and nuclear_found ==1:
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

def patchAverage(mesh,data,maxIds):
	numPts = mesh.GetNumberOfPoints()
	print(numPts)
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
			print maxIds
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

def patchAverageIgnoreNullValues(mesh,data,maxIds,nullValue):
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
			print maxIds
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
		sum = 0
		validPts = 0
		for i in xrange(0,patchPt_list.GetNumberOfIds()):
			val = getValue(mesh,patchPt_list.GetId(i),data)
			if(val != nullValue):
				sum = sum + val
				validPts = validPts + 1
		average = sum / validPts
		#set the interpolated value of that point
		data_array.InsertNextValue(average)
	return data_array
def intializeVTP(filename):
	datareader=vtk.vtkXMLPolyDataReader()
	datareader.SetFileName(filename + '.vtp')
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

def update_progress(progress, total):  
	print('\r[{0:10}]{1:>2}'.format('#' * int(progress * 10 /total), progress))

def main():
	filename = sys.argv[1]
	data_to_interpolate = sys.argv[2]
	if(len(sys.argv)==3):
		N = 100
		nullValue = -1
	elif(len(sys.argv)==4):
		N = int(sys.argv[3])
		nullValue = -1
	elif(len(sys.argv)==5):
		N = int(sys.argv[3])
		nullValue = int(sys.argv[4])
	mesh = intializeVTP(filename)
	smoothed_value = vtk.vtkDoubleArray()
	smoothed_value = patchAverageIgnoreNullValues(mesh,data_to_interpolate,N,nullValue)
	smoothed_value.SetName('Smoothed ' + data_to_interpolate)
	mesh.GetPointData().AddArray(smoothed_value)
	w = vtk.vtkXMLPolyDataWriter()
	w.SetInputData(mesh)
	w.SetFileName('patchSmoothing_testing.vtp')
	w.Write()


if __name__ == '__main__':
	main()