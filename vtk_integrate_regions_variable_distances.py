import sys
import os
import vtk
import numpy as np
import time
import glob
import math
import argparse

#Reads a vtu file with given filename and outputs the mesh object
def intializeVTU(filename,vprint):
	if '.vtu' not in filename:
		filename = filename + '.vtu'
	datareader=vtk.vtkXMLUnstructuredGridReader()
	datareader.SetFileName(filename)
	datareader.Update()

	mesh=datareader.GetOutput()
	vprint('Loaded .vtu file.')
	return mesh

#Reads a vtp file with given filename and outputs the mesh object
def intializeVTP(filename,vprint):
	if '.vtp' not in filename:
		filename = filename + '.vtp'
	datareader=vtk.vtkXMLPolyDataReader()
	datareader.SetFileName(filename)
	datareader.Update()

	mesh=datareader.GetOutput() 
	vprint('Loaded .vtp file.')
	return mesh

#Writes a vtp with given mesh object and filename
def writeVTP(mesh,filename,vprint):
	if '.vtp' not in filename:
		filename = filename + '.vtp'
	vprint('Writing ' + filename + ' file...')
	w = vtk.vtkXMLPolyDataWriter()
	w.SetInputData(mesh)
	w.SetFileName(filename)
	w.Write()
	vprint('done.')

#Writes a vtu with given mesh object and filename
def writeVTU(mesh,filename,vprint):
	if '.vtu' not in filename:
		filename = filename + '.vtu'
	vprint('Writing ' + filename + ' file...')
	w = vtk.vtkXMLUnstructuredGridWriter()
	w.SetInputData(mesh)
	w.SetFileName(filename)
	w.Write()
	vprint('done.')


# Calculates the magnitude of given vector (any length)
def calcMagnitude(vector):
	mag = 0
	for i in xrange(0,len(vector)):
		mag += vector[i]**2
	return mag ** (.5)

# Gets the value (magnitude) of scalar (tuple) on a given model at a given point from a given array name 
def getValue(model, ptID, factor):
	return calcMagnitude(model.GetPointData().GetArray(factor).GetTuple(ptID))


def getBoundaryNodes(mesh,radius,used_nodes,vprint):
	numPts = mesh.GetNumberOfPoints()
	boundary_nodes = vtk.vtkIdList()
	#First get the points within a sphere on the boundaries using the boundary points as seed points
	for ptID in xrange(0,numPts):
		if(mesh.GetPointData().GetArray('GlobalBoundaryPoints').GetValue(ptID)>0):
			x0,y0,z0 = mesh.GetPoint(ptID)

			neighbrs=vtk.vtkExtractGeometry()
			sphere = vtk.vtkSphere()
			sphere.SetCenter(x0,y0,z0)
			sphere.SetRadius(float(radius))
			neighbrs.SetImplicitFunction(sphere)
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
	writeVTU(neighbors,'test_extraction.vtu',vprint)
	return [boundary_nodes, used_nodes]

def extractRegion(mesh,selection_nodes,vprint):
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
		for i in xrange(0,cell_nodes.GetNumberOfIds()):
			cell_vtk_Id_list.InsertUniqueId(cell_nodes.GetId(i))

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

	#Outputs the mesh as an unstructured grid object
	output = vtk.vtkUnstructuredGrid()
	output.ShallowCopy(region.GetOutput())
	return output

#Adds an data array to mesh to enumerate the region based on the region_nodes list
def markRegionNodes(mesh,region_nodes):
	region_data = [0] * len(region_nodes)
	for i in xrange(0,len(region_data)):
		region_data[i] = []
	for i in xrange(0,len(region_data)):
		for j in xrange(0,len(region_nodes[i])):
			region_data[i].append(vtk.vtkDoubleArray())
	for ptID in xrange(0,mesh.GetNumberOfPoints()):
		counter = 1
		found = 0
		for i_region in xrange(0,len(region_nodes)):
			for j_region in xrange(0,len(region_nodes[i_region])):
				if(region_nodes[i_region][j_region].IsId(ptID) >= 0):
					region_data[i_region][j_region].InsertNextValue(1)
				else:
					region_data[i_region][j_region].InsertNextValue(0)
	for i in xrange(0,len(region_data)):
		for j in xrange(0,len(region_data[0])):
			region_data[i][j].SetName('Region_' + str(i) + '_' + str(j))
			mesh.GetPointData().AddArray(region_data[i][j])
	return mesh

def extractBoundaryPointsfromCells(mesh,boundary_nodes):
	numPts = mesh.GetNumberOfPoints()
	for ptID in xrange(0,numPts):
		cell_list = vtk.vtkIdList()
		temp_ptList = mesh.GetPointCells(ptID,cell_list)
		unique_values = vtk.vtkIdList()
		for cellID in xrange(0,cell_list.GetNumberOfIds()):
			unique_values.InsertUniqueId(cell_list.GetId(cellID))
		if(unique_values.GetNumberOfIds()>1):
			boundary_nodes.InsertUniqueId(ptID)

def getIntegratedQOI(model,sorted_region_nodes,quantities_to_integrate,output_filename,vprint):
	# Integrated quantities
	output_collection = []
	# Point Averaged quantities
	output2_collection = []
	for i_region in xrange(0, len(sorted_region_nodes)):
		for j_region in xrange(0,len(sorted_region_nodes[i_region])):
			print('Working on Region ' + str(i_region+1) + '_' + str(j_region+1) + '...')
			region = extractRegion(model,sorted_region_nodes[i_region][j_region],vprint)
			writeVTU(region,'test.vtu',vprint)
			numPts=region.GetNumberOfPoints()
			numCells=region.GetNumberOfCells()
		  
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
		    
				temp_cell = region.GetCell(icell)
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
					iid = pts_cell.GetId(ipt)
					pointCoordinate = region.GetPoint(iid) 
					for iq in xrange(0, len(quantities_to_integrate)):
						vector = QOI[iq].GetTuple(model.FindPoint(pointCoordinate))
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

	# First print a header that tells what each integrated quantity of interest is
	out_string = 'Region, Surface area, '
	for iq in xrange(0, len(quantities_to_integrate)):
		out_string = out_string + quantities_to_integrate[iq] + ', '
		out_string = out_string + quantities_to_integrate[iq] + "/Area, "
	out_string = out_string + '\n'
	outfile.write(out_string)
	# Print data
	counter = 0
	for i in xrange(0, len(sorted_region_nodes)):
		for j in xrange(0,len(sorted_region_nodes[i])):
			out_string = str(i) + '_' + str(j) + ',' + str(output_collection[counter][0])
			for iq in xrange(1, len(quantities_to_integrate)+1):
				out_string = out_string + ', ' + str(output_collection[counter][iq])
				out_string = out_string + ', ' + str(output_collection[counter][iq]/output_collection[i][0])
			counter += 1
		    
			out_string = out_string + '\n'
			outfile.write(out_string)

	outfile.write('\n')

	# Second print a header point averaged quantity of interests
	out_string = 'Region, '
	for iq in xrange(0, len(quantities_to_integrate)):
		out_string = out_string + quantities_to_integrate[iq] + ', '
	out_string = out_string + '\n'
	outfile.write(out_string)
	# Print data
	counter = 0
	for i in xrange(0, len(sorted_region_nodes)):
		for j in xrange(0,len(sorted_region_nodes[i])):
			out_string = str(i) + '_' + str(j) + ',' + str(output_collection[counter][0])
			for iq in xrange(0, len(quantities_to_integrate)):
				out_string = out_string + ', ' + str(output2_collection[counter][iq])
			counter += 1
		    
			out_string = out_string + '\n'
			outfile.write(out_string)

	outfile.close()

def sortRegionNodes(mesh,region_nodes,branch_groups,nuclear_factor,nuclear_threshold,gradient_factor,gradient_threshold,vprint):
	sorted_region_nodes = [0] * len(region_nodes)
	for i_region in xrange(0,len(sorted_region_nodes)):
		sorted_region_nodes[i_region] = []
	for i_region in xrange(0,len(region_nodes)):
		for i_branch_group in xrange(0,len(branch_groups)):
			sorted_region_nodes[i_region].append(vtk.vtkIdList())
	for i_region in xrange(0,len(region_nodes)):
		selection_nodes = region_nodes[i_region]
		for ptID in xrange(0,selection_nodes.GetNumberOfIds()):
			point = selection_nodes.GetId(ptID)
			for i_branch_group in xrange(0,len(branch_groups)):
				for i_branch in xrange(0,len(branch_groups[i_branch_group])):
					nuc_val = getValue(mesh,point,nuclear_factor)
					grad_val = getValue(mesh,point,gradient_factor)
					branch_val = getValue(mesh,point,branch_groups[i_branch_group][i_branch])
					if(branch_val==1 and nuc_val >= nuclear_threshold and grad_val <= gradient_threshold):
						sorted_region_nodes[i_region][i_branch_group].InsertNextId(point)
	for i_region in sorted_region_nodes[0]:
		vprint(i_region)
	return sorted_region_nodes



def labelBranches(mesh,branch_groups,vprint):
	for ibg in xrange(0,len(branch_groups)):
		for ib in xrange(0,len(branch_groups[ibg])):
			branch = intializeVTP(branch_groups[ibg][ib],vprint)
			branch_label = vtk.vtkDoubleArray()
			branch_pts = vtk.vtkIdList()
			for ip in xrange(0,branch.GetNumberOfPoints()):
				branch_pts.InsertNextId(mesh.FindPoint(branch.GetPoint(ip)))
			for ip in xrange(0,mesh.GetNumberOfPoints()):
				if(branch_pts.IsId(ip)>=0):
					branch_label.InsertNextValue(1)
				else:
					branch_label.InsertNextValue(0)
			branch_label.SetName(branch_groups[ibg][ib])
			mesh.GetPointData().AddArray(branch_label)
	writeVTP(mesh,'result_labeled_branches.vtp',vprint)
	return mesh
	

def createParser():
	parser = argparse.ArgumentParser(description='Determines the area and point averaged value of quantities of interest within a user-defined range of radii away from a bifurcation. NOTE: Model must include either "ModelFaceID" or "GlobalBoundaryPoints" array.')
	parser.add_argument('filename', type=str, help='the filename (include file ext) to analyze')
	parser.add_argument('-data', type=str, default='integrated_quantities.csv', help='the output data filename (default=integrated_quantities.csv)')
	parser.add_argument('radii', type=str, help='the list of radii desired separated by commas')
	parser.add_argument('-group1', type=str, nargs='?', action='append', default=[], help='list of branches in first group')
	parser.add_argument('-group2', type=str, nargs='?', action='append', default=[], help='list of branches in second group')
	parser.add_argument('-n', '--nucleus', type=str, nargs='?', default='Nucleus (Cx40Nuclei)',help='name of the array of Nucleus data')
	parser.add_argument('-nt', '--nuclear_threshold', type=float, default=20, help='the threshold of nucleus desired')
	parser.add_argument('-g', '--gradient', type=str, nargs='?', default='Gradients',help='name of the array of gradient data')
	parser.add_argument('-gt', '--gradient_threshold', type=float, default=10, help='the threshold of gradient desired')
	parser.add_argument('-iQOI', nargs='?', action='append', default=['vTAWSS (dynes/cm^2)', 'Gradients'],help='name of the array to be integrated')
	parser.add_argument('-vtp', type=str, nargs='?', const='marked_model', default='', help='the name of a vtp output file (without file ext) (default=none)')
	parser.add_argument('-v', '-verbose', type=int, nargs='?', const=1, default=0, help='turn on verbosity')
	return parser

def main(args):
	# Initialize variables
	FILENAME = args.filename
	OUTPUT_FILENAME = args.vtp + '.vtp'
	data_out_filename = args.data
	RADII = [int(i) for i in args.radii.decode("UTF-8").split(',')] 
	QOI = args.iQOI
	nuclear_factor = args.nucleus
	nuclear_threshold = args.nuclear_threshold
	gradient_factor = args.gradient
	gradient_threshold = args.gradient_threshold

	boundary_nodes = vtk.vtkIdList()
	if args.v:
		def vprint(*args):
			for arg in args:
				print(arg),
	else:
		vprint = lambda *a: None
	branch_groups = [args.group1, args.group2]
	vprint(branch_groups)
	mesh = intializeVTP(FILENAME,vprint)
	numPts = mesh.GetNumberOfPoints()
	# Extracts the nodes from each bifurcation and puts them into a vtkIdList
	if(mesh.GetPointData().HasArray('GlobalBoundaryPoints')==1):
		for ptID in xrange(0,numPts):
			if(mesh.GetPointData().GetArray('GlobalBoundaryPoints').GetValue(ptID)>0):
				boundary_nodes.InsertUniqueId(ptID)
	elif(mesh.GetCellData().HasArray('ModelFaceID')==1):
		extractBoundaryPointsfromCells(mesh,boundary_nodes)
	else:
		print('No boundary points on mesh or able to be extracted.')
		return 0

	# Gets all the nodes from 0 to the first radius specified.
	region_nodes = []
	used_nodes = vtk.vtkIdList()
	[first_region_nodes, used_nodes] = getBoundaryNodes(mesh,RADII[0],used_nodes,vprint)
	region_nodes.append(first_region_nodes)
	print('first region done.')

	# Gets all the nodes between each subsequent pair of radii.
	if(len(RADII)>1):
		for radius in xrange(1,len(RADII)):
				[single_region_nodes, used_nodes] = getBoundaryNodes(mesh,RADII[radius],used_nodes,vprint)
				print(single_region_nodes)
				region_nodes.append(single_region_nodes)
	furthest_region_nodes = vtk.vtkIdList()
	if(used_nodes.GetNumberOfIds() < numPts):
		for ptID in xrange(0,numPts):
			if(used_nodes.IsId(ptID) < 0):
				furthest_region_nodes.InsertNextId(ptID)
	print('all regions found.')
	region_nodes.append(furthest_region_nodes)
	
	labelBranches(mesh,branch_groups,vprint)
	# Sorts region nodes by branch and filters by Nuclei and WSSG
	sorted_region_nodes = sortRegionNodes(mesh,region_nodes,branch_groups,nuclear_factor,nuclear_threshold,gradient_factor,gradient_threshold,vprint)

	# Outputs a .vtp file with all regions marked visually under the "Regions" array
	if(OUTPUT_FILENAME!='.vtp'):
		mesh = markRegionNodes(mesh,sorted_region_nodes)
		writeVTP(mesh,OUTPUT_FILENAME,vprint)

	# Integrates each QOI for each region defined and generates data output file
	getIntegratedQOI(mesh,sorted_region_nodes,QOI,data_out_filename,vprint)


if __name__ == '__main__':
	parser = createParser()
	args = parser.parse_args()
	main(args)






		