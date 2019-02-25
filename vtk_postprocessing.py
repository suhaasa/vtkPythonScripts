import sys
import os
import vtk
import numpy as np
import time
import glob
import math
import argparse
from argparse import Namespace

def generateMasterFile(sv,fsi,pts,out):
	import vtk_mapping_interpolate_threshold_optimized as vtk_map
	import vtk_merge_mesh_PointData2 as vtk_mmPD2
	import vtk_WarpByVector as vtk_wbv
	parser = vtk_mmPD2.createParser()
	sub_args = parser.parse_args([pts,fsi, 'fsi_pts'])
	vtk_mmPD2.main(sub_args)
	parser = vtk_map.createParser()
	sub_args = parser.parse_args(['fsi_pts.vtp'])
	vtk_map.main(sub_args)
	parser = vtk_wbv.createParser()
	sub_args = parser.parse_args(['./','all_results_mapped_interpolated_threshold.vtp','./','-output_filename','warped.vtp'])
	vtk_wbv.main(sub_args)
	parser = vtk_mmPD2.createParser()
	sub_args = parser.parse_args(['warped.vtp',sv,'master_file'])
	vtk_mmPD2.main(sub_args)
	return

def createArgs():
	parser = argparse.ArgumentParser(description='Master postprocessing script.')
	parser.add_argument('-master','-master_file', type=str, nargs='?', default=None, help='master filename of all data')
	parser.add_argument('-raw','-raw_file', type=str, nargs='?', default=None, help='filename of mesh without WSSG and unit conversion')
	parser.add_argument('-sv','-sv_file', type=str, nargs='?', default=None, help='filename of mesh with rigid simulation data')
	parser.add_argument('-fsi','-fsi_file', type=str, nargs='?', default=None, help='filename of mesh with FSI data')
	parser.add_argument('-pts','-pts_file', type=str, nargs='?', default=None, help='filename of mesh with boundary point data')
	parser.add_argument('-output_filename', type=str, help='the output filename (without file ext)')
	parser.add_argument('-g', '-generate_master_file', type=int, nargs='?', const=1, default=0, help='generates the master vtp file with all data Required:fsi+raw/sv+pts')
	parser.add_argument('-irt', '-integrate_regions_threshold', type=int, nargs='?', const=1, default=0, help='run the vtk_integrate_regions_threshold algorithm Required:master')
	parser.add_argument('-irvd', '-integrate_regions_variable_distances', type=int, nargs='?', const=1, default=0, help='run the vtk_integrate_regions_variable_distances Required:master')
	parser.add_argument('-ib', '-integrate_branches', type=int, nargs='?', const=1, default=0, help='run the vtk_integrate_branches algorithm Required:fsi+raw/sv+branch')
	parser.add_argument('-c1Dq', '-compute1Dquantities', type=int, nargs='?', const=1, default=0, help='run the vtk_compute1Dquantities algorithm Required:fsi+raw/sv+centerline')
	parser.add_argument('-b', '-branches', type=str, nargs='?', const='./Branches/', default='./Branches/', help='location of the branch mesh files (can only contain those files)')
	parser.add_argument('-t', '-threshold', type=str, nargs='?', const='10,50,100', default='10,50,100', help='threshold of values to use')
	parser.add_argument('-r', '-radii', type=str, nargs='?', const='10,25,50', default='10,25,50', help='radii to use')
	parser.add_argument('-c', '-centerlines', type=str, nargs='?', const='./Centerlines/', default='./Centerlines/', help='location of the centerline vtp files (can only contain those files)')
	parser.add_argument('-v', '-verbose', type=int, nargs='?', const=1, default=0, help='turn on verbosity')
	parser.add_argument('-i', '-interactive', type=int, nargs='?', const=1, default=0, help='turn on interactive mode (not implemented yet)')
	
	return parser.parse_args()

def main(args):

	if args.v:
		def vprint(*args):
			# Print each argument separately so caller doesn't need to
			# stuff everything to be printed into a single string
			for arg in args:
				print(arg),
	else:
		vprint = lambda *a: None      # do-nothing function

	if(args.g):
		if(args.sv is None or args.fsi is None or args.pts is None):
			print('Error: need sv, fsi, and pts files (see help for more info)')
		else:
			generateMasterFile(args.sv,args.fsi,args.pts,args.output_filename)

	if(args.irt):
		if(args.master is None):
			print('Error: need master file (see help for more info)')
		else:
			import vtk_integrate_regions_threshold as vtk_irt
			parser=vtk_irt.createParser()
			sub_args = parser.parse_args([args.master,args.t,'-vtp',args.output_filename])
			vtk_irt.main(sub_args)

	if(args.irvd):
		if(args.master is None):
			print('Error: need master file (see help for more info)')
		else:
			import vtk_integrate_regions_variable_distances as vtk_irvd
			parser=vtk_irvd.createParser()
			sub_args = parser.parse_args([args.master,args.r])
			vtk_irvd.main(sub_args)

	if(args.ib):
		if(args.master is None):
			vprint('Error: need master file (see help for more info)')
		else:
			import vtk_integrate_branches as vtk_ib
			parser=vtk_ib.createParser()
			sub_args = parser.parse_args([args.master,args.b])
			vtk_ib.main(sub_args,vprint)

	if(args.c1Dq):
		if(args.master is None):
			vprint('Error: need master file (see help for more info)')
		else:
			import vtk_compute1Dquantities as vtk_c1Dq
			vtk_ib.main(args)




if __name__ == '__main__':
	args = createArgs()
	main(args)