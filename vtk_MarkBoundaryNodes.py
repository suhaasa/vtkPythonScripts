import sys
import os
import vtk
import numpy as np
import time
import glob
import math
import argparse

def intializeVTU(filename,vprint):
	datareader=vtk.vtkXMLUnstructuredGridReader()
	datareader.SetFileName(filename)
	datareader.Update()

	mesh=datareader.GetOutput()
	vprint('Loaded ' + filename + ' file.')
	return mesh

def intializeVTP(filename,vprint):
	datareader=vtk.vtkXMLPolyDataReader()
	datareader.SetFileName(filename)
	datareader.Update()

	mesh=datareader.GetOutput() 
	vprint('Loaded ' + filename + ' file.')
	return mesh

def writeVTP(mesh,filename,vprint):
	vprint('Writing ' + filename + ' file...')
	w = vtk.vtkXMLPolyDataWriter()
	w.SetInputData(mesh)
	w.SetFileName(filename + '.vtp')
	w.Write()
	vprint('done.')

def writeVTU(mesh,filename,vprint):
	vprint('Writing ' + filename + ' file...')
	w = vtk.vtkXMLUnstructuredGridWriter()
	w.SetInputData(mesh)
	w.SetFileName(filename + '.vtu')
	w.Write()
	vprint('done.')

def calcMagnitude(vector):
	mag = 0
	for i in xrange(0,len(vector)):
		mag += vector[i]**2
	return mag ** (.5)

def getValue(model, ptID, factor):
	return calcMagnitude(model.GetPointData().GetArray(factor).GetTuple(ptID))

def getBoundaryNodes(mesh,radius):
	numPts = mesh.GetNumberOfPoints()
	boundary_nodes = vtk.vtkIdList()
	#First get the points within a sphere on the boundaries using the boundary points as seed points
	for ptID in xrange(0,numPts):
		if(mesh.GetPointData().GetArray('GlobalBoundaryPoints').GetValue(ptID)>0):
			x0,y0,z0 = mesh.GetPoint(ptID)

			sphere = vtk.vtkSphere()
			sphere.SetCenter(x0,y0,z0)
			sphere.SetRadius(float(radius))
			neighbrs=vtk.vtkExtractGeometry()
			neighbrs.SetImplicitFunction(sphere)
			neighbrs.ExtractInsideOn()
			neighbrs.ExtractBoundaryCellsOn()
			neighbrs.SetInputData(mesh)
			neighbrs.Update()
			neighbors = neighbrs.GetOutput()
			writeVTU(neighbors,'test_boundary.vtu')

			for neighborID in xrange(0,neighbors.GetNumberOfPoints()):
				meshID = mesh.FindPoint(neighbors.GetPoint(neighborID))
				boundary_nodes.InsertUniqueId(meshID)
	return boundary_nodes

	parser.add_argument('filename', type=str, help='the filename (include file ext) to analyze')

def extractBoundaryPointsfromCells(mesh,boundary_nodes):
	numPts = mesh.GetNumberOfPoints()
	for ptID in xrange(0,numPts):
		cell_list = vtk.vtkIdList()
		temp_ptList = mesh.GetPointCells(ptID,cell_list)
		unique_values = vtk.vtkIdList()
		for cellID in xrange(0,cell_list.GetNumberOfIds()):
			unique_values.InsertUniqueId(mesh.GetCellData().GetArray('ModelFaceID').GetValue((cell_list.GetId(cellID))))
		if(unique_values.GetNumberOfIds()>1):
			boundary_nodes.InsertUniqueId(ptID)

def createParser():
	parser = argparse.ArgumentParser(description='Generates "GlobalBoundaryPoints" array from "ModelFaceID" cell data. ')
	parser.add_argument('filename', type=str, help='the filename (include file ext) to analyze')
	parser.add_argument('-vtp', type=str, nargs='?', const='marked_model', default='marked_model', help='the name of a vtp output file (without file ext) (default=marked_model)')
	parser.add_argument('-v', '-verbose', type=int, nargs='?', const=1, default=0, help='turn on verbosity')
	return parser

def main(args):
	FILENAME = args.filename
	OUTPUT_FILENAME = args.vtp
	if args.v:
		def vprint(*args):
			# Print each argument separately so caller doesn't need to
			# stuff everything to be printed into a single string
			for arg in args:
				print(arg)
	else:
		vprint = lambda *a: None      # do-nothing function
	mesh = intializeVTP(FILENAME,vprint)
	numPts = mesh.GetNumberOfPoints()
	boundary_nodes = vtk.vtkIdList()

	# Extracts the nodes from each bifurcation and puts them into a vtkIdList
	if(mesh.GetPointData().HasArray('GlobalBoundaryPoints')==1):
		for ptID in xrange(0,numPts):
			if(mesh.GetPointData().GetArray('GlobalBoundaryPoints').GetValue(ptID)>0):
				boundary_nodes.InsertUniqueId(ptID)
	elif(mesh.GetCellData().HasArray('ModelFaceID')==1):
		extractBoundaryPointsfromCells(mesh,boundary_nodes)
	else:
		vprint('No boundary points on mesh able to be extracted.')
		return 0

	# Generates "GlobalBoundaryPoints" array and adds it to the mesh
	boundary_nodes_array = vtk.vtkDoubleArray()
	for ptID in xrange(0,numPts):
		if(boundary_nodes.IsId(ptID) >= 0):
			boundary_nodes_array.InsertNextValue(1)
		else:
			boundary_nodes_array.InsertNextValue(0)
	boundary_nodes_array.SetName('GlobalBoundaryPoints')
	mesh.GetPointData().AddArray(boundary_nodes_array)

	# Write the mesh with boundary points marked
	writeVTP(mesh,OUTPUT_FILENAME,vprint)


if __name__ == '__main__':
	parser = createParser()
	args = parser.parse_args()
	main(args)
