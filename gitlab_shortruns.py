#!/usr/bin/python
# -*- coding: utf-8 -*-

import os
import sys
import math

os.chdir("ini")
cmd="../build/bin/gvec gitlab_shortruns.ini"
os.system(cmd) #+" 2>std_1st_sim.err 1>std_1st_sim.out")
restartFile="GITLAB_RUN_State_0000_99999999.dat"
if ( os.path.isfile(restartFile)):
  success= 1;
  cmd="../build/bin/gvec gitlab_shortruns_restart.ini " + restartFile
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

  
