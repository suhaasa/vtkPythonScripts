import sys
import os
import vtk
import numpy as np
import time
import glob
import math
import argparse

def intializeVTU(filename):
	datareader=vtk.vtkXMLUnstructuredGridReader()
	datareader.SetFileName(filename)
	datareader.Update()

	mesh=datareader.GetOutput()
	print('Loaded .vtu file.')
	return mesh

def intializeVTP(filename):
	datareader=vtk.vtkXMLPolyDataReader()
	datareader.SetFileName(filename)
	datareader.Update()

	mesh=datareader.GetOutput() 
	print('Loaded .vtp file.')
	return mesh

def writeVTP(mesh,filename):
	print('Writing .vtp file...')
	w = vtk.vtkXMLPolyDataWriter()
	w.SetInputData(mesh)
	w.SetFileName(filename)
	w.Write()
	print('done.')

def writeVTU(mesh,filename):
	print('Writing .vtu file...')
	w = vtk.vtkXMLUnstructuredGridWriter()
	w.SetInputData(mesh)
	w.SetFileName(filename)
	w.Write()
	print('done.')

def calcMagnitude(vector):
	mag = 0
	for i in xrange(0,len(vector)):
		mag += vector[i]**2
	return mag ** (.5)

def getValue(model, ptID, factor):
	return calcMagnitude(model.GetPointData().GetArray(factor).GetTuple(ptID))

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
	parser = argparse.ArgumentParser(description='Returns .vtp file of mesh with global boundary nodes labeled using ModelFaceID.')
	parser.add_argument('input_filename', type=str, help='the input filename (include file ext)')
	parser.add_argument('-output_filename', type=str, nargs='?', default='output.vtp', help='the output filename (include file ext)')
	return parser

def main(args):
	FILENAME = args.input_filename
	OUTPUT_FILENAME = args.output_filename
	mesh = intializeVTP(FILENAME)
	boundary_nodes = vtk.vtkIdList()
	extractBoundaryPointsfromCells(mesh,boundary_nodes)
	data_array = vtk.vtkDoubleArray()
	data_array.SetName('GlobalBoundaryPoints')
	for ptID in xrange(0,mesh.GetNumberOfPoints()):
		if(boundary_nodes.IsId(ptID)>0):
			data_array.InsertNextValue(1)
		else:
			data_array.InsertNextValue(0)
	mesh.GetPointData().AddArray(data_array)
	writeVTP(mesh,OUTPUT_FILENAME)

if __name__ == '__main__':
	parser = createParser()
	args = parser.parse_args()
	main(args)
