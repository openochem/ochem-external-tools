#!/usr/bin/env python
# encoding: utf-8
"""
cfp.py: Determines circular fingerprints from MOL2 files.

See cfp.py --help for details on usage.

Created by Florian Nigsch on 2007-08-16.
Copyright (c) 2007. All rights reserved.

This script calculates circular fingerprints based on Sybyl atom types.
The resulting fingerprints are similar to the MOLPRINT 2D fingerprints
implemented by Andreas Bender, see http://www.molprint.com/ for
details.

In the current script it is possible to specify the depth in bonds to
which the neighbouring atoms are considered. It is also possible to
name all molecules according to a given label, and they can be
consecutively numbered.

"""

from __future__ import print_function
import sys, os, re
import getopt
#import time
#from timeit import Timer
from time import *

help_message = '''Usage information
Determines circular fingerprints from MOL2 files.
Available options:
	-h, --help		output this help message
	-i, --input=	input file (required)
	-o, --output=	output file (default: cfp.out)
	-d, --depth=	maximum distance in bonds from atoms (default: 3)
	-p, --prefix=	prefix to use as label instead of the one in
			the MOL2 file
	-n, --numbering	turn on consecutive numbering of
			added labels (default: off)
	-v, --verbose	turn on verbose output (not implemented)
	-s, --separate	separate output for all distances (default: off)
'''

re_ATOMBLOCK = re.compile(r'@<TRIPOS>ATOM')
re_BONDBLOCK = re.compile(r'@<TRIPOS>BOND')
re_MOLECULE = re.compile(r'@<TRIPOS>MOLECULE')
re_BONDLINE = re.compile(r'([0-9]+\s*){3,}')
re_ATOMLINE = re.compile(r'[0-9]{1,3}\s*[A-Za-z]{1,2}[0-9]{0,3}\*{0,2}(\s*[-0-9]+\.[0-9]+)+.*$')

SeparateDepths = False

def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

def ParseAtomBlock(atoms, num_atoms):
	atomlist = []
	#print (atoms[1:])
	for line in atoms:
		line = line.strip()
		spline = re.split(r'\s+', line)
		atomlist.append(spline[5])
	if len(atoms) == num_atoms:
		return atomlist
	else:
		print ("Atomblock not correctly parsed. Skipping molecule.")
		return False

	
def ParseBondBlock(bonds, num_bonds):
	bondlist = []
	for line in bonds:
		line = line.strip()
		spline = re.split(r'\s+', line)
		try:
			bondlist.append([int(spline[1]),int(spline[2]),spline[3]])
		except IndexError as e:
			print ("Wrong bondline!")
			print (bonds)
			#print (spline)
			return False
	if len(bondlist) == num_bonds:
		return bondlist
	else:
		print ("Bondblock not correctly parsed. Skipping molecule.")
		return False


def GetNeighbourMap(AtomIDs, Atoms, Bonds):
	
	neighbours = {}
	neighboursID = {}
	
	
	for AtomID in AtomIDs:
		neighbours[AtomID] = []
		neighboursID[AtomID] = []
		for bondnum, bond in enumerate(Bonds):
			if (AtomID) in bond[0:2]:
				#print ("Atom %d in bond %d (%s from %s to %s)" % (AtomID, bondnum+1, bond[2], Atoms[bond[0]-1], Atoms[bond[1]-1]))
				if AtomID == bond[0]:
					neighbours[AtomID].append(Atoms[bond[1]-1])
					neighboursID[AtomID].append(bond[1])
				elif AtomID == bond[1]:
					neighbours[AtomID].append(Atoms[bond[0]-1])
					neighboursID[AtomID].append(bond[0])
			
	return neighbours, neighboursID


def CountOccurencesOfAtoms(atomlist):
	res = {}
	for atom in atomlist:
		if atom not in res.keys():
			res[atom] = int(1)
		else:
			res[atom] += 1
	
	return res
		
	
def GetFingerPrint(Atoms, Bonds, maxdist=3):
	nlevel = {}
	nmap, nlevel[1] = GetNeighbourMap(range(1,len(Atoms)+1), Atoms, Bonds)
	for dist in range(2, maxdist+1):
		nlevel[dist] = {}
		for key in nlevel[dist-1].keys():
			secondneighbours = []
			for elem in [nlevel[1][x] for x in nlevel[dist-1][key]]:
				secondneighbours += elem
			nset = set(secondneighbours)
			for remove in range(1,dist):
				nset = nset.difference(nlevel[remove][key])
			try:
				nset.remove(key)
			except KeyError:
				pass
			nlevel[dist][key] = list(nset)
	
	fp = ''
	Attributes = set()
	for key in nlevel[1].keys():
		if not SeparateDepths:
			fp += str(Atoms[key-1])
		for dist in range(1, maxdist+1):
			counts = CountOccurencesOfAtoms([Atoms[x-1] for x in nlevel[dist][key]])
			if SeparateDepths:
				fp += str(Atoms[key-1])
			for ckey in counts.keys():
				chunk = ";%d-%s-%d" % (dist, ckey, counts[ckey])
				fp += chunk
				Attributes.add(chunk)
			if SeparateDepths:
				fp += "\t"
		if not SeparateDepths:
			fp += "\t"
	return fp[:-1], len(Attributes)



class Usage(Exception):
	def __init__(self, msg):
		self.msg = msg		

class CFPError(Exception):
	pass

class CFPError_Input(CFPError):
	def __init__(self, message):
		self.message = message


class CFPError_Output(CFPError):
	def __init__(self, message):
		self.message = message


def main(argv=None):
	if argv is None:
		argv = sys.argv
	try:
		try:
			opts, args = getopt.getopt(argv[1:], "ho:vi:d:p:ns", ["help", "output=", "input=", "prefix=", "depth=", "numbering", "separate"])
		except getopt.error as msg:
			raise Usage(msg)
	
		# option processing
		for option, value in opts:
			if option == "-v":
				verbose = True
			if option in ("-h", "--help"):
				raise Usage(help_message)
			if option in ("-o", "--output"):
				fn_output = value
			if option in ("-i", "--input"):
				fn_input = value
			if option in ("-d", "--depth"):
				depth = int(value)
			if option in ("-p", "--prefix"):
				prefix = value
			if option in ("-n", "--numbering"):
				UsePrefixNumbering = True
			if option in ("-s", "--separate"):
				SeparateDepths = True
	
	except Usage as err:
		eprint(sys.argv[0].split("/")[-1] + ": " + str(err.msg))
		eprint("\tfor help use --help")
		return 2
	
	try:
		try:
			stat_v = os.stat(fn_input)
		#except OSError as (errno, strerror):
			#print ("%s: Input I/O Error(%s): %s" % (sys.argv[0].split("/")[-1], errno, strerror))
		except OSError:
			print ("%s: Input I/O Error(%s)" % (sys.argv[0].split("/")[-1], errno, str(OSError.msg)))
	except UnboundLocalError:
		print ("%s: No input file specified." % sys.argv[0].split("/")[-1])
		try:
			raise Usage(help_message)
		except Usage as err:
			eprint (sys.argv[0].split("/")[-1] + ": " + str(err.msg))
			eprint  ("\tfor help use --help")
			return 2
		
	try:
		len(fn_output)
	except UnboundLocalError:
		print ("%s: No output file specified. Writing to cfp.out." % sys.argv[0].split("/")[-1])
		fn_output = "cfp.out"
	
	try:
		if (depth > 6) or (depth < 1):
			print ("%s: Valid values for depth are from 1 to 6." % sys.argv[0].split("/")[-1])
			try:
				raise Usage(help_message)
			except Usage as err:
				eprint (sys.argv[0].split("/")[-1] + ": " + str(err.msg))
				eprint ("\tfor help use --help")
				return 2
	except UnboundLocalError:
		print ("%s: No depth specified. Using default value of 3." % sys.argv[0].split("/")[-1])
		depth = 3
	
	try:
		if prefix != "":
			UsePrefix = True
	except UnboundLocalError:
		UsePrefix = False
		pass
	
	try:
		if UsePrefixNumbering:
			pass
	except UnboundLocalError:
		UsePrefixNumbering = False
	
	atoms = []
	bonds = []
	read_ab, read_bb = False, False
	NextLineIsLabel = False
	NextLineIsInfo = False
	SkipThisMolecule = False
	outfile = open(fn_output, 'w')
	CountConverted = 0
	TotalAttributes = 0
	starttime = time()
	MoleculeNumber = 0
	SkippedMolecules = 0
	
	############
	# Database #
	############
	
	for line in open(fn_input, 'r'):
		
		if NextLineIsLabel:
			NextLineIsLabel = False
			label = line[:-1]
			NextLineIsInfo = True
			continue
		
		if NextLineIsInfo:
			NextLineIsInfo = False
			line = line.strip()
			spline = re.split(r' +', line)
			num_atoms = int(spline[0])
			num_bonds = int(spline[1])
			if num_atoms < 2:
				sys.stdout.write('*')
				sys.stdout.flush()
				SkipThisMolecule = True
				SkippedMolecules += 1
				continue
					
		if re_MOLECULE.match(line):
			MoleculeNumber += 1
			atoms = []
			bonds = []
			NextLineIsLabel = True
			SkipThisMolecule = False
			continue
		
		if SkipThisMolecule:
			continue
		
		if re_ATOMBLOCK.match(line):
			read_ab = True
			read_bb = False
			emptyline = False
			continue
			
		if re_BONDBLOCK.match(line):
			read_ab = False
			read_bb = True
			emptyline = False
			continue
			
		if re.match(r'^\r\n$', line):
			emptyline = True
		if re.match(r'^$', line):
			emptyline = True
		
		if read_ab and not read_bb:
			if not emptyline and not SkipThisMolecule:
				if re_ATOMLINE.search(line):
					atoms.append(line)
				else:
					sys.stderr.write("Invalid atom line:\n")
					sys.stderr.write(line)
					SkipThisMolecule = True
					read_ab = False
					continue
		
		if read_bb and not read_ab:
			if not emptyline:
				if re_BONDLINE.search(line):
					bonds.append(line)
				else:
					sys.stderr.write("Invalid bond line!\n")
					sys.stderr.write(line)
					SkipThisMolecule = True
					read_bb = False
					continue
		
		if len(atoms) == num_atoms and len(bonds) == num_bonds:
			Atoms = ParseAtomBlock(atoms, num_atoms)
			if not Atoms:
				sys.stderr.write("ERROR: Not able to parse atom block. Exiting.\n")
				sys.exit(1)
			Bonds = ParseBondBlock(bonds, num_bonds)
			if not Bonds:
				sys.stderr.write("ERROR: Not able to parse bond block. Exiting.\n")
				sys.exit(1)
			Fingerprint, NumberOfAttributes = GetFingerPrint(Atoms, Bonds, depth)
			TotalAttributes += NumberOfAttributes
			CountConverted += 1
			if UsePrefix:
				label = prefix
				if UsePrefixNumbering:
					label += '-' + str(CountConverted)
			outfile.write( "%s\t%s\n" % (label, Fingerprint) )
			if not CountConverted % 25:
				sys.stdout.write('.')
				sys.stdout.flush()
			atoms = []
			bonds = []
			read_ab = False
			read_bb = False
	print ("\n")
	print ("%d molecules converted in %.2f seconds (%5.2f molecules per second)." % (CountConverted, time()-starttime, CountConverted/(time()-starttime)))
	print ("In average there are %.2f unique features per molecule." % (float(TotalAttributes)/CountConverted))
	#print ("CPU time used: %.2f seconds (%.4f per molecule)" % (time.clock(), time.clock()/CountConverted))
	print ("Molecules that have been skipped (<2 atoms): %d" % SkippedMolecules)
	outfile.close()
	return 0

if __name__ == "__main__":
	sys.exit(main())
