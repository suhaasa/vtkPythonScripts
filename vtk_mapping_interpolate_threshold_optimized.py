# --- VTK-PYTHON SCRIPT FOR READING MESH VTP AND DACH1 VTI TO MAP DACH1 VALUES ONTO MESH SURFACE
# --- WRITTEN BY SUHAAS ANBAZHAKAN
# --- BASED ON A SCRIPT BY AEKAANSH VERMA

import sys
import vtk
import numpy as np
from vtk.util import numpy_support as VN
import argparse

def main(args):
	WSS_CONV_FACTOR = 10000
	P_CONV_FACTOR = .133333333333333
	INTERPOLATING_RADIUS = 5
	WALL_THRESHOLD = 30
	NUCLEAR_THRESHOLD = 30

	# First, read in the .vtp file for mesh and Dach1
	datareader=vtk.vtkXMLPolyDataReader()
	datareader.SetFileName(args.filename)
	datareader.Update()

	ref = vtk.vtkXMLImageDataReader()
	ref.SetFileName('ImageData_Combined.vti')
	ref.Update()

	# Read your data into another polydata variable for reading
	mesh=vtk.vtkPolyData()
	mesh=datareader.GetOutput()

	imageData=vtk.vtkPolyData()
	imageData=ref.GetOutput()

	numPts=mesh.GetNumberOfPoints()

	factors = ['Dach1']
	nuclearFactor = 'Cx40Nuclei'

	for molecularFactor in xrange(0, len(factors)):
		# Create new Dach1 Array to fill
		maxMFArray = vtk.vtkDoubleArray()
		meanMFArray = vtk.vtkDoubleArray()
		thresholdMFArray = vtk.vtkDoubleArray()
		nuclearMFArray = vtk.vtkDoubleArray()
		NuclearArray = vtk.vtkDoubleArray()
		# Loop through mesh points, find closest point in Dach1.vti, and add to Dach1 Array

		for meshPointID in xrange(0, numPts):
			meshPointCoordinate = mesh.GetPoint(meshPointID)
			x,y,z = mesh.GetPoint(meshPointID)
			pointsPool = np.array([[x,y,z]]);
			pointsPool = np.append(pointsPool,[[x+INTERPOLATING_RADIUS,y,z],[x-INTERPOLATING_RADIUS,y,z]],axis=0)
			pointsPool = np.append(pointsPool,[[x,y+INTERPOLATING_RADIUS,z],[x,y-INTERPOLATING_RADIUS,z]],axis=0)
			pointsPool = np.append(pointsPool,[[x,y,z+INTERPOLATING_RADIUS],[x,y,z-INTERPOLATING_RADIUS]],axis=0)
			maxMFValue=0
			meanMFValue=0
			x = pointsPool[0,0]
			y = pointsPool[0,1]
			z = pointsPool[0,2]
			meshPointCoordinate = x,y,z
			thresholdWallValue = imageData.GetPointData().GetArray('ImageFile').GetValue(imageData.FindPoint(meshPointCoordinate))
			for i in xrange(0,len(pointsPool)):
				x = pointsPool[i,0]
				y = pointsPool[i,1]
				z = pointsPool[i,2]
				meshPointCoordinate = x,y,z
				meshPointID = imageData.FindPoint(meshPointCoordinate)
				MFValue = imageData.GetPointData().GetArray(factors[molecularFactor]).GetValue(meshPointID)
				meanMFValue +=MFValue

				WallValue = imageData.GetPointData().GetArray('ImageFile').GetValue(meshPointID)

				if maxMFValue<MFValue:
					maxMFValue = MFValue

			if thresholdWallValue>WALL_THRESHOLD:
				x = pointsPool[0,0]
				y = pointsPool[0,1]
				z = pointsPool[0,2]
				meshPointCoordinate = x,y,z
				thresholdMFArray.InsertNextValue(imageData.GetPointData().GetArray(factors[molecularFactor]).GetValue(imageData.FindPoint(meshPointCoordinate)))
			else:
				thresholdMFArray.InsertNextValue(0)
			if imageData.GetPointData().HasArray(nuclearFactor)==1:
				thresholdNuclearValue = imageData.GetPointData().GetArray(nuclearFactor).GetValue(imageData.FindPoint(meshPointCoordinate))
				NuclearArray.InsertNextValue(thresholdNuclearValue)
				if thresholdNuclearValue>NUCLEAR_THRESHOLD:
					x = pointsPool[0,0]
					y = pointsPool[0,1]
					z = pointsPool[0,2]
					meshPointCoordinate = x,y,z
					if mesh.GetPointData().HasArray('Filtered ' + nuclearFactor)==1:
						NuclearMultiplier = mesh.GetPointData().GetArray('Filtered '+nuclearFactor).GetValue(meshPointCoordinate)/thresholdNuclearValue
					else:
						NuclearMultiplier = 1
					nuclearMFArray.InsertNextValue(NuclearMultiplier * imageData.GetPointData().GetArray(factors[molecularFactor]).GetValue(imageData.FindPoint(meshPointCoordinate)))
				else:
					nuclearMFArray.InsertNextValue(-1)

			meanMFValue = meanMFValue/len(pointsPool)
			maxMFArray.InsertNextValue(maxMFValue)
			meanMFArray.InsertNextValue(meanMFValue)

		# Add Dach1 Array to point data
		maxMFArray.SetName('Max ' + factors[molecularFactor] + ' (' + str(INTERPOLATING_RADIUS) + ')') 
		mesh.GetPointData().AddArray(maxMFArray)
		meanMFArray.SetName('Mean ' + factors[molecularFactor] + ' (' + str(INTERPOLATING_RADIUS)+ ')') 
		mesh.GetPointData().AddArray(meanMFArray)
		thresholdMFArray.SetName('Threshold ' + factors[molecularFactor] + ' (' + str(INTERPOLATING_RADIUS)+ ')') 
		mesh.GetPointData().AddArray(thresholdMFArray)
		if imageData.GetPointData().HasArray(nuclearFactor)==1:
			nuclearMFArray.SetName('Nuclear ' + factors[molecularFactor] + ' (' + nuclearFactor + ')') 
			mesh.GetPointData().AddArray(nuclearMFArray)
			NuclearArray.SetName('Nucleus ' + '(' + nuclearFactor + ')') 
			mesh.GetPointData().AddArray(NuclearArray)
	# Convert TA_WSS and Pressure to micron units and new arrays to point data
	if mesh.GetPointData().HasArray('vTAWSS')==1:
		TA_WSS = VN.vtk_to_numpy(mesh.GetPointData().GetArray('vTAWSS'))
		TA_WSS_numpy = WSS_CONV_FACTOR*TA_WSS
		TA_WSS_vtk = VN.numpy_to_vtk(TA_WSS_numpy)
		TA_WSS_vtk.SetName('vTAWSS (dynes/cm^2)') 
		mesh.GetPointData().AddArray(TA_WSS_vtk)
	else:
		print('vTAWSS not found. No conversion done.')

	if mesh.GetPointData().HasArray('pressure_avg')==1:
		pressure_avg = VN.vtk_to_numpy(mesh.GetPointData().GetArray('pressure_avg'))
		pressure_numpy = pressure_avg/P_CONV_FACTOR
		pressure_vtk = VN.numpy_to_vtk(pressure_numpy)
		pressure_vtk.SetName('pressure_avg (mmHg)') 
		mesh.GetPointData().AddArray(pressure_vtk) 
	else:
		print('pressure_avg not found. No conversion done.')

	# Write a new .vtp file that has all converted values and mapped values.
	w = vtk.vtkXMLPolyDataWriter()
	w.SetInputData(mesh)
	w.SetFileName('all_results_mapped_interpolated_threshold.vtp')
	w.Write()

def createParser():
	parser = argparse.ArgumentParser(description='Master postprocessing script.')
	parser.add_argument('filename', type=str, help='filename of mesh')
	return parser

if __name__ == '__main__':
	parser = createParser()
	args = parser.parse_args()
	main(args)