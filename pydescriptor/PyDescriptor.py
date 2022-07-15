####################################################################
# Authors: Dr. Vijay Masand (Amravati, India) and Dr. Vesna Rastija (Osijek, Croatia) 
# Funded by: Josip Juraj Strossmayer University of Osijek, Croatia.
# Copyright (c) 2016, Dr. Vijay Masand and Dr. Vesna Rastija
# All rights reserved. 
# The users can modify and distribute it with the condition to keep the names of original developers and funding university as it is. 
# It is unlawful to remove the names of original developers and funding university. 
# This script is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or 
# FITNESS FOR A PARTICULAR PURPOSE. Read the tutorial carefully before using this manuscript. Developed using PyMOl 1.8.2
# This script calculates >16k molecular descriptors of a molecule 
# & saves them into the file "PyDescriptor.csv".
# Modified and adapted to batch analysis by I.V. Tetko, itetko@bigchem.de
####################################################################
import pymol,glob,os,csv
import os.path
from pymol.wizard import Wizard
t = ('C','N','O','S','P','F','Cl','Br','I','acc','don','lipo','ringC','ringN','ringO','ringS','H','all')
tt = ('spC','sp2C','sp3C','aroC','spN','sp2N','sp3N','plaN','amdN','aroN','sp2O','sp3O','acidO','sp2S','sp3S','oxidS','sulfonS','sp3P','F','Cl','Br','I')
nam = [['mol_name'],['nBonds'],['nRings'],['rsa'],['fmsaplus'],['fmsaminus'],['fsasaplus'],['fsasaminus'],['AbSA'],['molweight'],['avg_molweight'],['nN_no_H']]
f = open('PyDescriptor.csv','a+')
for row in nam:
 for columns in row:
  f.write('%s,'%columns)
for m in t:
 xy = 'number of %s'%m
 m_hy = '%s_hy'%m
 m_hy1 = '%s_hy1'%m
 m_hy2 = '%s_hy2'%m
 m_hy3 = '%s_hy3'%m
 m_hy4 = '%s_hy4'%m
 m_hy5 = '%s_hy5'%m
 nmpcminus = 'n%spcminus'%m
 nmpcplus = 'n%spcplus'%m
 m_SASA = '%s_SASA' %m
 m_MSA = '%s_MSA' %m
 m_MSA1 = '%s_MSA1'%m
 m_MSA2 = '%s_MSA2'%m
 m_MSA3 = '%s_MSA3'%m
 m_MSA4 = '%s_MSA4'%m
 m_HASA = '%s_HASA'%m
 m_HASA1 = '%s_HASA1'%m
 m_HASA2 = '%s_HASA2'%m
 m_HASA3 = '%s_HASA3'%m
 m_HASA4 = '%s_HASA4'%m
 m_sumpc = '%s_sumpc'%m
 m_AbSA = '%s_AbSA'%m
 mminus_SASA = '%sminus_SASA'%m
 mminus_MSA = '%sminus_MSA'%m
 mminus_sumpc = '%sminus_sumpc'%m
 mminus_AbSA = '%sminus_AbSA'%m
 mplus_SASA = '%splus_SASA'%m
 mplus_MSA = '%splus_MSA'%m
 mplus_sumpc = '%splus_sumpc'%m
 mplus_AbSA = '%splus_AbSA'%m
 dec = [[xy],[m_hy],[m_hy1],[m_hy2],[m_hy3],[m_hy4],[m_hy5],[nmpcminus],[nmpcplus],[m_SASA],[m_MSA],[m_MSA1],[m_MSA2],[m_MSA3],[m_MSA4],[m_HASA],[m_HASA1],[m_HASA2],[m_HASA3],[m_HASA4],[m_sumpc],[m_AbSA],[mminus_SASA],[mminus_MSA],[mminus_sumpc],[mminus_AbSA],[mplus_SASA],[mplus_MSA],[mplus_sumpc],[mplus_AbSA]]
 for row in dec:
  for columns in row:
   f.write('%s,'%columns)
for m in t[:16]:
 for i in range(1,10):
  mdab = 'da_%s_%dB'%(m,i)
  mdaa = 'da_%s_%dA'%(m,i)
  madb = 'ad_%s_%dB'%(m,i)
  mada = 'da_%s_%dA'%(m,i)
  mplusb = 'plus_%s_%dB'%(m,i)
  mplusa = 'plus_%s_%dA'%(m,i)
  mminusb = 'minus_%s_%dB'%(m,i)
  mminusa = 'minus_%s_%dA'%(m,i)
  mcoma = 'com_%s_%dA'%(m,i)
  mcomaplus = 'com_%splus_%dA'%(m,i)
  mcomaminus = 'com_%sminus_%dA'%(m,i)
  mcomahy = 'com_%shyd_%dA'%(m,i)
  des = [[madb],[mada],[mplusb],[mplusa],[mminusb],[mminusa],[mcoma],[mcomaplus],[mcomaminus],[mcomahy]]
  for row in des:
   for columns in row:
    f.write('%s,'%columns)
for i in range(1,10):
 for m in t[:16]:
  for v in t[:16]:
   mvb = '%s_%s_%dB'%(m,v,i)
   mva = '%s_%s_%dA'%(m,v,i)
   mvbc = '%s_%s_%dBc'%(m,v,i)
   mvac = '%s_%s_%dAc'%(m,v,i)
   fmvb = 'f%s%s%dB'%(v,m,i)
   fmva = 'f%s%s%dA'%(m,v,i)
   data = [[mvb],[mva],[mvbc],[mvac],[fmvb],[fmva]]
   for row in data:
    for columns in row:
     f.write('%s,' %columns)
for m in tt:
 xy = 'number of %s'%m
 m_hy = '%s_hy'%m
 m_hy1 = '%s_hy1'%m
 m_hy2 = '%s_hy2'%m
 m_hy3 = '%s_hy3'%m
 m_hy4 = '%s_hy4'%m
 m_hy5 = '%s_hy5'%m
 dej = [[xy],[m_hy],[m_hy1],[m_hy2],[m_hy3],[m_hy4],[m_hy5]]
 for row in dej:
  for columns in row:
   f.write('%s,'%columns)
for m in tt:
 for i in range(1,6):
  mallb = '%s_all_%dB'%(m,i)
  malla = '%s_all_%dA'%(m,i)
  mcoma = 'com_%s_%dA'%(m,i)
  deg = [[mallb],[malla],[mcoma]]
  for row in deg:
   for columns in row:
    f.write('%s,'%columns)
wr = csv.writer(f)
wr.writerow('')


import timeit
start_time = timeit.default_timer()
from pymol import cmd,stored, util
from chempy import Bond
path = os.path.dirname(pymol.__script__)
cmd.delete('all')
mole = glob.glob(os.path.join(path, '*.mol2'))
for file in mole:
 cmd.load(file)
 ma = cmd.get_model()
 mol_name = cmd.get_names("objects")
 t = ('e. C','e. N','e. O','e. S','e. P','e. F','e. Cl','e. Br','e. I','acc','don','e. C & not (neighbor e. N+O) or e. S & not (neighbor hydro or e. N+O) or e. Cl or e. I or e. Br','e. C & (byring e. C)','e. N & (byring e. N)','e. O & (byring e. O)','e. S & (byring e. S)','e. H','all &!(com)')
 tt = ('tt. C.1','tt. C.2','tt. C.3','tt. C.ar','tt. N.1','tt. N.2','tt. N.3','tt. N.pl3','tt. N.am','tt. N.ar','tt. O.2','tt. O.3','tt. O.co2','tt. S.2','tt. S.3','tt. S.o','tt. S.o2','tt. P.3','e. F','e. Cl','e. Br','e. I')
 cmd.set('dot_density','4')
 cn = cmd.count_atoms
 ca = util.get_area
 cs = util.get_sasa
 cp = util.sum_partial_charges
 MSA = ca("(all)", 1, 1, quiet=1)
 MSAplus = ca('pc.>0', 1, 1, quiet=1)
 MSAminus = ca('pc.<0', 1, 1, quiet=1)
 SASA = cs("(all)", quiet=1, _self=cmd)
 SASAplus = cs('pc.>0', quiet=1, _self=cmd)
 SASAminus = cs('pc.<0', quiet=1, _self=cmd)
 try:rsa = MSA / SASA
 except ZeroDivisionError: rsa = 0
 try:fmsaplus = MSAplus / MSA
 except ZeroDivisionError: fmsaplus = 0
 try:fmsaminus = MSAminus / MSA
 except ZeroDivisionError: fmsaminus = 0
 try:fsasaplus = SASAplus / SASA
 except ZeroDivisionError: fsasaplus = 0
 try:fsasaminus = SASAminus / SASA
 except ZeroDivisionError: fsasaminus = 0
 AbSA = abs(SASA - MSA)
 molweight = util.compute_mass("(all)", implicit=True,quiet=1,_self=cmd)
 avg_molweight = float(molweight) / len(ma.atom)
 nAtoms = cn('(all)')
 nBonds = len(ma.bond)
 nRings = nBonds - nAtoms + 1
 nN_no_H = cn('e. N & !(neighbor hydro)')
 com = cmd.pseudoatom("com", pos=cmd.centerofmass('all'))
 nam = [[mol_name],[nBonds],[nRings],[rsa],[fmsaplus],[fmsaminus],[fsasaplus],[fsasaminus],[AbSA],[molweight],[avg_molweight],[nN_no_H]]
 for row in nam:
  for columns in row:
   f.write('%s,'%columns)
 for m in t:
  xy = cn('%s'%m)
  m_hy = cn('%s & (pc.>-0.2 & pc.<0.2)'%m)
  m_hy1 = cn('%s & (pc.>-0.199999 & pc.<-0.100000)'%m)
  m_hy2 = cn('%s & (pc.>-0.099999 & pc.<0)'%m)
  m_hy3 = cn('%s & (pc.>0 & pc.<0.099999)'%m)
  m_hy4 = cn('%s & (pc.>0.100000 & pc.<0.199999)'%m)
  m_hy5 = cn ('%s & (pc.>0.3 or pc.<-0.3)'%m)
  nmpcminus = cn ('%s & (pc.<0)'%m)
  nmpcplus = cn('%s and (pc.>0)'%m)
  m_SASA = cs('%s'%m, quiet=1, _self=cmd)
  m_MSA = ca('%s'%m, quiet=1, _self=cmd)
  m_MSA1 = ca('%s &(pc.>-0.199999 & pc.<-0.100000)'%m, quiet=1, _self=cmd)
  m_MSA2 = ca('%s &(pc.>-0.099999 & pc.<0)'%m, quiet=1, _self=cmd)
  m_MSA3 = ca('%s &(pc.>0 & pc.<0.099999)'%m, quiet=1, _self=cmd)
  m_MSA4 = ca('%s &(pc.>0.100000 & pc.<0.199999)'%m, quiet=1, _self=cmd)
  m_HASA = cs('%s &(pc.>-0.2 & pc.<0.2)'%m, quiet=1, _self=cmd)
  m_HASA1 = cs('%s &(pc.>0.000000 & pc.<0.099999)'%m, quiet=1, _self=cmd)
  m_HASA2 = cs('%s &(pc.>0.100000 & pc.<0.200000)'%m, quiet=1, _self=cmd)
  m_HASA3 = cs('%s &(pc.<0.000000 & pc.>-0.099999)'%m, quiet=1, _self=cmd)
  m_HASA4 = cs('%s &(pc.<-0.100000 & pc.>-0.200000)'%m, quiet=1, _self=cmd)
  m_sumpc = cp('%s'%m, 1, _self=cmd)
  m_AbSA = abs(m_SASA - m_MSA)
  mminus_SASA = cs('%s &(pc.<0)'%m, quiet=1, _self=cmd)
  mminus_MSA = ca('%s &(pc.<0)'%m, quiet=1, _self=cmd)
  mminus_sumpc = cp('%s &(pc.<0)'%m, 1, _self=cmd)
  mminus_AbSA = abs(mminus_SASA - mminus_MSA)
  mplus_SASA = cs('%s &(pc.>0)'%m, quiet=1, _self=cmd)
  mplus_MSA = ca('%s &(pc.>0)'%m, quiet=1, _self=cmd)
  mplus_sumpc = cp('%s &(pc.>0)'%m, 1, _self=cmd)
  mplus_AbSA = abs(mplus_SASA - mplus_MSA)
  dec = [[xy],[m_hy],[m_hy1],[m_hy2],[m_hy3],[m_hy4],[m_hy5],[nmpcminus],[nmpcplus],[m_SASA],[m_MSA],[m_MSA1],[m_MSA2],[m_MSA3],[m_MSA4],[m_HASA],[m_HASA1],[m_HASA2],[m_HASA3],[m_HASA4],[m_sumpc],[m_AbSA],[mminus_SASA],[mminus_MSA],[mminus_sumpc],[mminus_AbSA],[mplus_SASA],[mplus_MSA],[mplus_sumpc],[mplus_AbSA]]
  for row in dec:
   for columns in row:
    f.write('%s,'%columns) #line number 167
 for m in t[:16]:
  for i in range(1,10):
   mdab = cn("acc & don xt. %d & %s"%(i,m))
   mdaa = cn("%s w. %d of acc & don"%(m,i))
   madb = cn("acc or don xt. %d & %s"%(i,m))
   mada = cn("%s w. %d of acc or don"%(m,i))
   mplusb = cn("%s xt. %d & pc.>0"%(m,i))
   mplusa = cn("%s w. %d of pc.>0"%(m,i))
   mminusb = cn("%s xt. %d & pc.<0"%(m,i))
   mminusa = cn("%s w. %d of pc.<0"%(m,i))
   mcoma = cn("%s w. %d of com"%(m,i))
   mcomaplus = cn("%s & pc.>0 w. %d of com"%(m,i))
   mcomaminus = cn("%s & pc.<0 w. %d of com"%(m,i))
   mcomahy = cn("%s & (pc.<0.2 & pc.>-0.2) w. %d of com"%(m,i))
   des = [[madb],[mada],[mplusb],[mplusa],[mminusb],[mminusa],[mcoma],[mcomaplus],[mcomaminus],[mcomahy]]
   for row in des:
    for columns in row:
     f.write('%s,'%columns)
 for i in range(1,10):
  for m in t[:16]:
   for v in t[:16]:
    mvb = cn("%s xt. %d & %s"%(v,i,m))
    mva = cn("%s w. %d of %s"%(v,i,m))
    mvbc = cp("%s xt. %d & %s"%(v,i,m), 1, _self=cmd)
    mvac = cp("%s w. %d of %s"%(v,i,m), 1, _self=cmd)
    fmvb = cn("%s xt. %d & %s &!(%s xt. %d & %s)"%(v,i,m,v,i-1,m))
    fmva = cn("%s w. %d & %s &!(%s w. %d & %s)"%(v,i,m,v,i-1,m))
    data = [[mvb],[mva],[mvbc],[mvac],[fmvb],[fmva]]
    for row in data:
     for columns in row:
      f.write('%s,'%columns)
 for m in tt:
  py = cn('%s'%m)                                 # line number 220
  m_hy = cn('%s & (pc.>-0.2 & pc.<0.2)'%m)
  m_hy1 = cn('%s & (pc.>-0.199999 & pc.<-0.100000)'%m)
  m_hy2 = cn('%s & (pc.>-0.099999 & pc.<0)'%m)
  m_hy3 = cn('%s & (pc.>0 & pc.<0.099999)'%m)
  m_hy4 = cn('%s & (pc.>0.100000 & pc.<0.199999)'%m)
  m_hy5 = cn ('%s & (pc.>0.3 or pc.<-0.3)'%m)
  dej = [[py],[m_hy],[m_hy1],[m_hy2],[m_hy3],[m_hy4],[m_hy5]]
  for row in dej:
   for columns in row:
    f.write('%s,'%columns)
 for m in tt:
  for i in range(1,6):
   mallb = cn("%s xt. %d & (all &! e. H)"%(m,i))
   malla = cn("(all &! e. H &! com) w. %d of %s"%(i,m))
   mcoma = cn("%s w. %d of com"%(m,i))
   deg = [[mallb],[malla],[mcoma]]
   for row in deg:
    for columns in row:
     f.write('%s,'%columns)
 wr = csv.writer(f)
 wr.writerow('')
 print ("The data for molecule have been appended to PyDescriptors.csv", mol_name)
 cmd.delete('all')
f.close()
elapsed = timeit.default_timer() - start_time
print ('Total time elapsed for molecular descriptor calculation', elapsed)
print ('The data for all molecules have been appended to PyDescriptors.csv ')
