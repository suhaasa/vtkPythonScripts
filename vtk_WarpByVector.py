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
	w = vtk.vtkXMLPolyDataWriter()
	w.SetInputData(mesh)
	w.SetFileName(filename)
	w.Write()

def writeVTU(mesh,filename):
	w = vtk.vtkXMLUnstructuredGridWriter()
	w.SetInputData(mesh)
	w.SetFileName(filename)
	w.Write()

def calcMagnitude(vector):
	mag = 0
	for i in xrange(0,len(vector)):
		mag += vector[i]**2
	return mag ** (.5)

def getValue(model, ptID, factor):
	return calcMagnitude(model.GetPointData().GetArray(factor).GetTuple(ptID))

def createParser():
	parser = argparse.ArgumentParser(description='Warps a mesh based on a vector at each point.')
	parser.add_argument('indir', type=str, help='where to find the input file')
	parser.add_argument('input_filename', type=str, help='the input filename (include file ext)')
	parser.add_argument('outdir', type=str, help='where to put the output file')
	parser.add_argument('-output_filename', type=str, nargs='?', default=None, help='the prefix for the output filename')
	parser.add_argument('warp_vec_name', type=str, nargs='?', default='MS_Displacement', help='name of the array that holds the displacement vector')
	return parser

def main(args):
	
	INDIR = args.indir
	FILENAME = args.input_filename
	OUTDIR = args.outdir
	WARP_VECTOR_NAME = args.warp_vec_name

	mesh = intializeVTP(INDIR + FILENAME)
	mesh.GetPointData().SetActiveVectors(WARP_VECTOR_NAME)
	warpVector = vtk.vtkWarpVector()
	warpVector.SetInputData(mesh)
	warpVector.Update()
	if(args.output_filename is None):
		output_filename = OUTDIR + '/Warped_' + FILENAME
	else:
		output_filename = OUTDIR + args.output_filename
	writeVTP(warpVector.GetOutput(), output_filename)

if __name__ == '__main__':
	parser = createParser()
	args = parser.parse_args()
	main(args)
