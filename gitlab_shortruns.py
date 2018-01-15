#!/usr/bin/python
# -*- coding: utf-8 -*-

import os
import sys
import math
import argparse

########################################################################################################

# MAIN PROGRAM

########################################################################################################
parser = argparse.ArgumentParser(description='# Tool to submit short runs of GVEC\n'
                                             '# example execution command:\n'
                                             '# python gitlab_shortruns.py builddir',\
                                 formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument('builddir',type=str, 
                                 help='Path to builddir from .gitlab.yml folder, mandatory last argument')

args = parser.parse_args()

builddir=args.builddir
os.chdir("ini")
cmd="../"+builddir+"build/bin/gvec gitlab_shortruns.ini"
os.system(cmd) #+" 2>std_1st_sim.err 1>std_1st_sim.out")
restartFile="GITLAB_RUN_State_0000_99999999.dat"
if ( os.path.isfile(restartFile)):
  success= 1;
  cmd="../"+builddir+"/bin/gvec gitlab_shortruns_restart.ini " + restartFile
  os.system(cmd) #+" 2>std_restart.err 1>std_restart.out")
  
  if ( os.path.isfile("GITLAB_RESTART_State_0001_99999999.dat")) :
    print( " \n \n \n")
    print("---> first simulation & restart test did execute successfully!" )
    sys.exit(0);
  else :
    print( " \n \n \n")
    print( "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" )
    print( "!!!! ---> restart simulation test did not work        !!!" )
    print( "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" )
    sys.exit(100);
else :
  print( " \n \n \n")
  print( "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" )
  print( "!!!! ---> first simulation test did not work          !!!" )
  print( "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" )
  sys.exit(100);

  
