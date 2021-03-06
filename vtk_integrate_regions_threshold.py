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
	vprint('Loaded .vtu file.')
	return mesh

def intializeVTP(filename,vprint):
	datareader=vtk.vtkXMLPolyDataReader()
	datareader.SetFileName(filename)
	datareader.Update()

	mesh=datareader.GetOutput() 
	vprint('Loaded .vtp file.')
	return mesh

def writeVTP(mesh,filename,vprint):
	vprint('Writing ' + filename + ' file...')
	w = vtk.vtkXMLPolyDataWriter()
	w.SetInputData(mesh)
	w.SetFileName(filename)
	w.Write()
	vprint('done.')

def writeVTU(mesh,filename,vprint):
	vprint('Writing ' + filename + ' file...')
	w = vtk.vtkXMLUnstructuredGridWriter()
	w.SetInputData(mesh)
	w.SetFileName(filename)
	w.Write()
	vprint('done.')

def calcMagnitude(vector):
	mag = 0
	for i in xrange(0,len(vector)):
		mag += vector[i]**2
	return mag ** (.5)

def getValue(model, ptID, factor):
	if(model.GetPointData().HasArray(factor)):
		return calcMagnitude(model.GetPointData().GetArray(factor).GetTuple(ptID))
	else:
		print('Array name not found: ' + factor)
		return 0
		

def getBoundaryNodes(mesh,radius1,radius2,used_nodes):
	numPts = mesh.GetNumberOfPoints()
	boundary_nodes = vtk.vtkIdList()
	#First get the points within a sphere on the boundaries using the boundary points as seed points
	for ptID in xrange(0,numPts):
		if(mesh.GetPointData().GetArray('GlobalBoundaryPoints').GetValue(ptID)>0):
			x0,y0,z0 = mesh.GetPoint(ptID)

			neighbrs=vtk.vtkExtractGeometry()
			sphere1 = vtk.vtkSphere()
			sphere1.SetCenter(x0,y0,z0)
			sphere1.SetRadius(float(radius1))
			sphere2 = vtk.vtkSphere()
			sphere2.SetCenter(x0,y0,z0)
			sphere2.SetRadius(float(radius2))
			diff = vtk.vtkImplicitBoolean()
			diff.AddFunction(sphere2)
			diff.AddFunction(sphere1)
			diff.SetOperationTypeToDifference()
			neighbrs.SetImplicitFunction(sphere2)
			neighbrs.ExtractInsideOn()
			neighbrs.ExtractBoundaryCellsOn()
			neighbrs.SetInputData(mesh)
			neighbrs.Update()

			neighbors= vtk.vtkConnectivityFilter()
			neighbors.SetInputData(neighbrs.GetOutput())
			neighbors.SetExtractionModeToClosestPointRegion()
			neighbors.SetClosestPoint(x0,y0,z0)
			neighbors.Update()
			neighbors = neighbors.GetOutput()
			for neighborID in xrange(0,neighbors.GetNumberOfPoints()):
				meshID = mesh.FindPoint(neighbors.GetPoint(neighborID))
				if(used_nodes.IsId(meshID) < 0):
					boundary_nodes.InsertNextId(meshID)
					used_nodes.InsertNextId(meshID)
	writeVTU(neighbors,'test_extraction.vtu')
	return [boundary_nodes, used_nodes]

def extractRegion(mesh,selection_nodes):
	#Intialize variables
	ids = vtk.vtkIdTypeArray()
	cell_nodes = vtk.vtkIdList()
	cell_vtk_Id_list = vtk.vtkIdList()
	cellIds = vtk.vtkIdTypeArray()
	ids.SetNumberOfComponents(1)
	
	#Determines the cells enclosed by selection_nodes (which are points)
	for i in xrange(0,selection_nodes.GetNumberOfIds()):
		ids.InsertNextValue(selection_nodes.GetId(i))
		mesh.GetPointCells(selection_nodes.GetId(i), cell_nodes)
		for j in xrange(0,cell_nodes.GetNumberOfIds()):
			cell_vtk_Id_list.InsertUniqueId(cell_nodes.GetId(j))

	#Converts the vtkIdList into vtkIdTypeArray
	for i in xrange(0,cell_vtk_Id_list.GetNumberOfIds()):
		cellIds.InsertNextValue(cell_vtk_Id_list.GetId(i))

	#Creates the selection object to extract the subset of cells from the mesh
	region=vtk.vtkExtractSelection()
	region.SetInputData(0,mesh)
	tempCells = vtk.vtkSelectionNode()
	tempCells.SetFieldType(vtk.vtkSelectionNode.CELL)
	tempCells.SetContentType(vtk.vtkSelectionNode.INDICES)
	tempCells.SetSelectionList(cellIds)
	tempSelection = vtk.vtkSelection()
	tempSelection.AddNode(tempCells)
	region.SetInputData(1,tempSelection)
	region.Update()

	#Outputs the mesh as an Polydata object
	output = vtk.vtkPolyData()
	output.ShallowCopy(region.GetOutput())
	return output,cell_vtk_Id_list

def markRegionNodes(mesh,region_nodes):
	region_data = vtk.vtkDoubleArray()
	for ptID in xrange(0,mesh.GetNumberOfPoints()):
		region_val = 0
		for i_region in xrange(0,len(region_nodes)):
			if(region_nodes[i_region].IsId(ptID) >= 0):
				region_val = float(i_region+1)
		region_data.InsertNextValue(region_val)
	region_data.SetName('Regions')
	mesh.GetPointData().AddArray(region_data)
	return mesh

def extractBoundaryPointsfromCells(mesh,boundary_nodes):
	numPts = mesh.GetNumberOfPoints()
	for ptID in xrange(0,numPts):
		temp_ptList = mesh.GetPointCells(mesh.GetPoint(ptID),cell_list)
		unique_values = vtk.vtkIdList()
		for cellID in xrange(0,cell_list.GetNumberOfIds()):
			unique_values.InsertUniqueId(cell_list.GetId(cellID))
		if(unique_values.GetNumberOfIds()>1):
			boundary_nodes.InsertUniqueId(mesh.GetPoint(ptID))

def getIntegratedQOI(model,region_nodes,quantities_to_integrate,output_filename,thresholds,percentiles,vprint):
	# Integrated quantities
	output_collection = []
	# Point Averaged quantities
	output2_collection = []
	# loop through each set of nodes defining a region
	for i_region in xrange(0, len(region_nodes)):
		print('Working on Region ' + str(i_region+1) + '...')
		region, cell_list = extractRegion(model,region_nodes[i_region])
		writeVTP(region,'test.vtp',vprint)
		numPts=region.GetNumberOfPoints()
		numCells=cell_list.GetNumberOfIds()

	  	# Read in your quantities of interest from the .vtp file
		QOI = []
		for i in xrange(0, len(quantities_to_integrate)):
			QOI.append(vtk.vtkDoubleArray())
			QOI[i] = model.GetPointData().GetArray(quantities_to_integrate[i])

	  	# Initialize data structures to store area, point total, and point counts
		total_area = 0.0
		integrated_variables = []
		for i in xrange(0, len(quantities_to_integrate)):
			integrated_variables.append(0.0)
	  
		point_totals = []
		point_counts = []
		point_averages = []
		for iq in xrange(0, len(quantities_to_integrate)):
			point_totals.append(0.0)
			point_counts.append(0)
			point_averages.append(0.0)

		for pointID in xrange(0, numPts):
			# Get the region coordinate
			pointCoordinate = region.GetPoint(pointID) 
			for iq in xrange(0, len(quantities_to_integrate)):
				# Find the region coordinate in the model
				vector = QOI[iq].GetTuple(model.FindPoint(pointCoordinate))
				mag = 0
				for i in xrange(0,len(vector)):
					mag = mag + float(vector[i])**2
					mag = mag**(.5)
					# Ignore 0 value points
					if(mag>0):
						point_totals[iq] = point_totals[iq] + mag
						point_counts[iq] = point_counts[iq] + 1

		# Calculate point averaged QOIs
		for iq in xrange(0, len(quantities_to_integrate)):
			if(point_counts[iq]>0):
				point_averages[iq] = point_totals[iq]/point_counts[iq]

	  # Now loop over the cells and add in contribution to the integrated_variable
	  # that is equal to the average of all the cell nodal values multiplied by
	  # the area of the cell
		for icell in xrange(0,numCells):
	    
			temp_cell = model.GetCell(cell_list.GetId(icell))
			pts_cell = temp_cell.GetPointIds()
	    
			# First, get the area of this cell
			vtkpt = temp_cell.GetPoints()
			p0 = vtkpt.GetPoint(0)
			p1 = vtkpt.GetPoint(1)
			p2 = vtkpt.GetPoint(2)
			temp_area = vtk.vtkTriangle().TriangleArea(p0,p1,p2)
			total_area = total_area + temp_area
	    
			# Now, sum up the values of the quantities of interest at the cell
			# vertices
			averages = []
			for iq in xrange(0, len(quantities_to_integrate)):
				averages.append(0.0)
	    
			for ipt in xrange(0,pts_cell.GetNumberOfIds()):
				for iq in xrange(0, len(quantities_to_integrate)):
					vector = QOI[iq].GetTuple(pts_cell.GetId(ipt))
					mag = 0
					for i in xrange(0,len(vector)):
						mag = mag + float(vector[i])**2
					mag = mag**(.5)
					averages[iq] = averages[iq] + mag
	        
	    # To complete the trapezoidal rule integration, multiply each summed quantity
	    # by the area of the cell, then divide by the number of vertices
			for iq in xrange(0,len(quantities_to_integrate)):
				integrated_variables[iq] = integrated_variables[iq] + \
				averages[iq] * temp_area / float(pts_cell.GetNumberOfIds())
	        
	  # Now that we have integrated the variables of interest, it is time to save
	  # the results to a list for outputting later
		temp_collection = [total_area]
		for iq in xrange(0, len(quantities_to_integrate)):
			temp_collection.append(integrated_variables[iq])

		temp2_collection = []
		for iq in xrange(0, len(quantities_to_integrate)):
			temp2_collection.append(point_averages[iq])
	  # Stores the integrated values
		output_collection.append(temp_collection)
	  # Stores the point averaged values
		output2_collection.append(temp2_collection)
	
	# Now that we have looped over all our .vtp files of interest and integrated
	# the variables, it is time to save them to the output file. We also print
	# out the integrated quantities divided by the surface area for post-processing
	# convenience

	outfile = open(output_filename, 'w')
	## This section prints out the point average quantities in the format:
	# Region | Surface Area | Integrated QOI[0] | Integrated QOI[0] / Area ... | Integrated QOI[n] | Integrated QOI[n] / Area
	# First print a header that tells what each integrated quantity of interest is
	out_string = 'Region, Surface area, '
	for iq in xrange(0, len(quantities_to_integrate)):
		out_string = out_string + quantities_to_integrate[iq] + ', '
		out_string = out_string + quantities_to_integrate[iq] + "/Area, "
	out_string = out_string + '\n'
	outfile.write(out_string)
	# Print data
	for i in xrange(0, len(region_nodes)):
		out_string = str(i) + '_' + str(thresholds[i]) + ',' + str(output_collection[i][0])
		for iq in xrange(1, len(quantities_to_integrate)+1):
			out_string = out_string + ', ' + str(output_collection[i][iq])
			if(output_collection[i][0]!=0):
				out_string = out_string + ', ' + str(output_collection[i][iq]/output_collection[i][0])
			else:
				out_string = out_string + ', ' + str('Nan')
	    
		out_string = out_string + '\n'
		outfile.write(out_string)

	outfile.write('\n')

	## This section prints out the point average quantities in the format:
	# Region | point averaged QOI[0] | point average QOI[1] ... | QOI[n] 
	# Second print a header point averaged quantity of interests
	out_string = 'Region, '
	for iq in xrange(0, len(quantities_to_integrate)):
		out_string = out_string + quantities_to_integrate[iq] + ', '
	out_string = out_string + '\n'
	outfile.write(out_string)
	# Print data
	for i in xrange(0, len(region_nodes)):
		out_string = str(i) + '_' + str(thresholds[i])
		for iq in xrange(0, len(quantities_to_integrate)):
			out_string = out_string + ', ' + str(output2_collection[i][iq])
	    
		out_string = out_string + '\n'
		outfile.write(out_string)

	outfile.close()

# Returns the a vtkIdList of nodes that have a QOI value between the low and high threshold
def getRegionNodes(mesh,lower_threshold,higher_threshold,QOI,nuclear_threshold,nuclear_factor):
	numPts = mesh.GetNumberOfPoints()
	node_list = vtk.vtkIdList()
	for ptID in xrange(0,numPts):
		QOI_value = getValue(mesh, ptID, QOI)
		nuclear_value = getValue(mesh, ptID, nuclear_factor)
		if (QOI_value >= lower_threshold and QOI_value <= higher_threshold and nuclear_value >= nuclear_threshold):
			node_list.InsertNextId(ptID)
	return node_list

def calcQuartiles(mesh,tQOI, percentiles, vprint):
	numPts = mesh.GetNumberOfPoints()
	values = []
	for ptID in xrange(0,numPts):
		values.append(getValue(mesh,ptID,tQOI))
	out = np.percentile(values,percentiles)
	vprint('Percentile values: ', out)
	return out

def createParser():
	parser = argparse.ArgumentParser(description='Determines the area and point averaged value of quantities of interest within a user-defined range of radii away from a bifurcation. NOTE: Model must include either "ModelFaceID" or "GlobalBoundaryPoints" array.')
	parser.add_argument('filename', type=str, help='the filename (include file ext) to analyze')
	parser.add_argument('-data', type=str, default='threshold_integrated_quantities.csv', help='the output data filename (default=integrated_quantities.csv)')
	parser.add_argument('-t', '--threshold', type=float, action='append', default=[], help='the list of thresholds desired')
	parser.add_argument('-p', '--percentiles', type=float, nargs='?', action='append', default=[], help='the list of percentiles desired (overwrites threshold values)')
	parser.add_argument('-vtp', type=str, nargs='?', default='', help='the name of a vtp output file without file ext (default=none)')
	parser.add_argument('-tQOI', type=str, nargs='?', action='append', default=['Gradients','vTAWSS (dynes/cm^2)'],help='name of the array to be thresholded')
	parser.add_argument('-n', '--nucleus', type=str, nargs='?', default='Nucleus (Cx40Nuclei)',help='name of the array of Nucleus data')
	parser.add_argument('-nt', '--nuclear_threshold', type=float, default=20, help='the threshold of nucleus desired')
	parser.add_argument('-iQOI', nargs='?', action='append', default=['vTAWSS (dynes/cm^2)', 'Gradients'],help='name of the array to be integrated')
	parser.add_argument('-v', '--verbose', action='store_true', help='turn on verbosity')	
	return parser

def main(args):
	FILENAME = args.filename
	OUTPUT_FILENAME = args.vtp + '.vtp'
	data_filename = args.data
	THRESHOLD_QOI = args.tQOI
	INTERGRATION_QOI = args.iQOI
	nuclear_factor = args.nucleus
	nuclear_threshold = args.nuclear_threshold

	if args.verbose:
		def vprint(*args):
			# Print each argument separately so caller doesn't need to
			# stuff everything to be printed into a single string
			for arg in args:
				print(arg),
	else:
		vprint = lambda *a: None      # do-nothing function

	mesh = intializeVTP(FILENAME,vprint)
	numPts = mesh.GetNumberOfPoints()

	for QOI_i in xrange(0,len(THRESHOLD_QOI)):
		if(len(args.percentiles)==0):
			THRESHOLD_LEVEL = args.threshold
		else:
			percentiles = args.percentiles
			THRESHOLD_LEVEL = calcQuartiles(mesh,THRESHOLD_QOI[QOI_i], percentiles, vprint)

		region_nodes = []
		first_region_nodes = getRegionNodes(mesh,0,THRESHOLD_LEVEL[0],THRESHOLD_QOI[QOI_i], nuclear_threshold, nuclear_factor)
		vprint(first_region_nodes)
		region_nodes.append(first_region_nodes)
		vprint('first region done.')
		for i in xrange(1,len(THRESHOLD_LEVEL)):
			single_region_nodes = getRegionNodes(mesh,THRESHOLD_LEVEL[i-1],THRESHOLD_LEVEL[i],THRESHOLD_QOI[QOI_i], nuclear_threshold, nuclear_factor)
			vprint(single_region_nodes)
			region_nodes.append(single_region_nodes)
		vprint('all regions found.')

		if(OUTPUT_FILENAME!='.vtp'):
			mesh = markRegionNodes(mesh,region_nodes)
			writeVTP(mesh,OUTPUT_FILENAME,vprint)

		getIntegratedQOI(mesh,region_nodes,INTERGRATION_QOI,THRESHOLD_QOI[QOI_i].split(' ')[0] + '_' + data_filename,THRESHOLD_LEVEL,args.percentiles,vprint)

if __name__ == '__main__':
	parser = createParser()
	args = parser.parse_args()
	main(args)