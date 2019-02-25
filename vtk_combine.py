# --- VTK-PYTHON SCRIPT FOR READING MESH VTP AND DACH1 VTI TO MAP DACH1 VALUES ONTO MESH SURFACE
# --- WRITTEN BY SUHAAS ANBAZHAKAN
# --- BASED ON A SCRIPT BY AEKAANSH VERMA

import sys
import vtk
import numpy as np
from vtk.util import numpy_support as VN

# change this input varible to include the file names of all the .vti files

factors = ['Jag1','Dach1','DAPI']

# -----------------------------------
# DO NOT CHANGE ANTHING BELOW THIS LINE
# -----------------------------------

# First, read in inital image data
datareader=vtk.vtkXMLImageDataReader()
datareader.SetFileName(factors[0] + '.vti')
datareader.Update()

image=vtk.vtkPolyData()
image=datareader.GetOutput()

# Then, loop through each factor and add it to the inital data set
for MF in xrange(0,len(factors)):
	#Read vti file
	ref = vtk.vtkXMLImageDataReader()
	ref.SetFileName(factors[MF] + '.vti')
	ref.Update()

	#Read your data into another polydata variable for reading
	factor=vtk.vtkPolyData()
	factor=ref.GetOutput()

	#Adds image data as new array
	ImageData = VN.vtk_to_numpy(factor.GetPointData().GetArray(0))
	ImageData_vtk = VN.numpy_to_vtk(ImageData)
	ImageData_vtk.SetName(factors[MF]) 
	image.GetPointData().AddArray(ImageData_vtk) 


# Write a new .vti file that has all converted values and mapped values.
w = vtk.vtkXMLImageDataWriter()
w.SetInputData(image)
w.SetFileName('ImageData_Combined.vti')
w.Write()