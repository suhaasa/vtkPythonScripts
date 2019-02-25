import sys
import vtk
import numpy as np
from vtk.util import numpy_support as VN
import argparse

def intializeVTU(filename):
	datareader=vtk.vtkXMLUnstructuredGridReader()
	datareader.SetFileName(filename)
	datareader.Update()

	mesh=vtk.vtkDataSetMapper()
	mesh=datareader.GetOutput()
	print('Loaded .vtu file.')
	return mesh

def intializeVTP(filename):
	datareader=vtk.vtkXMLPolyDataReader()
	datareader.SetFileName(filename)
	datareader.Update()

	mesh=vtk.vtkDataSetMapper()
	mesh=datareader.GetOutput()
	print('Loaded .vtp file.')
	return mesh

def writeVTP(mesh,filename):
	w = vtk.vtkXMLPolyDataWriter()
	w.SetInputData(mesh)
	w.SetFileName(filename + '.vtp')
	w.Write()

def writeVTU(mesh,filename):
	w = vtk.vtkXMLUnstructuredGridWriter()
	w.SetInputData(mesh)
	w.SetFileName(filename + '.vtu')
	w.Write()

def transferData2(mesh_with_less_data_arrays, mesh_with_more_data_arrays):
	numPts1 = mesh_with_more_data_arrays.GetNumberOfPoints()
	numPts2 = mesh_with_less_data_arrays.GetNumberOfPoints()
	numArrays = mesh_with_less_data_arrays.GetPointData().GetNumberOfArrays()
	for i_array in xrange(0,numArrays):
		print('Working on Array ' + str(i_array+1) + '/' + str(numArrays))
		data_array = vtk.vtkDoubleArray()
		for i_node in xrange(0,numPts1):
			globalID = mesh_with_more_data_arrays.GetPointData().GetArray('GlobalNodeID').GetValue(i_node)
			for g_node in xrange(0,numPts2):
				if(mesh_with_less_data_arrays.GetPointData().GetArray('GlobalNodeID').GetValue(g_node) == globalID):
					data_array.InsertNextValue(mesh_with_less_data_arrays.GetPointData().GetArray(i_array).GetValue(i_node))			
		data_array.SetName(mesh_with_less_data_arrays.GetPointData().GetArrayName(i_array))
		mesh_with_more_data_arrays.GetPointData().AddArray(data_array)
	return mesh_with_more_data_arrays

def transferData(mesh_with_less_data_arrays, mesh_with_more_data_arrays):
	numPts = mesh_with_more_data_arrays.GetNumberOfPoints()
	numArrays = mesh_with_less_data_arrays.GetPointData().GetNumberOfArrays()
	for i_array in xrange(0,numArrays):
		print('Working on Array ' + str(i_array+1) + '/' + str(numArrays))
		data = mesh_with_less_data_arrays.GetPointData().GetArray(i_array).GetTuple(0)
		if(len(data)>1):
			data_array = vtk.vtkFloatArray()
			data_array.SetNumberOfComponents(3)
		else:
			data_array = vtk.vtkDoubleArray()
		for i_node in xrange(0,numPts):
			node_to_copy_from = mesh_with_less_data_arrays.FindPoint(mesh_with_more_data_arrays.GetPoint(i_node))
			data = mesh_with_less_data_arrays.GetPointData().GetArray(i_array).GetTuple(node_to_copy_from)
			if(len(data)>1):
				data_array.InsertNextTuple3(data[0],data[1],data[2])	
			else:
				data_array.InsertNextValue(data[0])
		data_array.SetName(mesh_with_less_data_arrays.GetPointData().GetArrayName(i_array))
		mesh_with_more_data_arrays.GetPointData().AddArray(data_array)
	return mesh_with_more_data_arrays

def createParser():
	parser = argparse.ArgumentParser(description='Merges the data arrays on two similar but different meshes. Can also be used to transfer data from one mesh to another.')
	parser.add_argument('filename1', type=str, help='the first filename (include file ext) to merge')
	parser.add_argument('filename2', type=str, help='the second filename (include file ext) to merge')
	parser.add_argument('output_filename', type=str, help='the output filename (without file ext)')
	parser.add_argument('output', type=int, choices=[1,2], nargs='?', default=2, help='choose which file (1 or 2) that the output data will be mapped on')
	return parser

def main(args):
	FILENAME1 = args.filename1
	FILENAME2 = args.filename2
	output_filename = args.output_filename
	FILE_TO_OUTPUT = args.output
	if(FILENAME1[len(FILENAME1)-3:] == 'vtu' and FILENAME2[len(FILENAME2)-3:] == 'vtu'):
		mesh1 = intializeVTU(FILENAME1)
		mesh2 = intializeVTU(FILENAME2)
		file_type = '.vtu'
	elif(FILENAME1[len(FILENAME1)-3:] == 'vtp' and FILENAME2[len(FILENAME2)-3:] == 'vtp'):
		mesh1 = intializeVTP(FILENAME1)
		mesh2 = intializeVTP(FILENAME2)
		file_type = '.vtp'
	else:
		print('error: filenames in different format or not recongnized.')

	if(FILE_TO_OUTPUT == 1):
		merged_mesh = transferData(mesh2, mesh1)
	else:
		merged_mesh = transferData(mesh1, mesh2)
	if(file_type == '.vtu'):
		writeVTU(merged_mesh, output_filename)
	else:
		writeVTP(merged_mesh, output_filename)

if __name__ == '__main__':
	parser = createParser()
	args = parser.parse_args()
	main(args)
