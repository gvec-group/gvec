#!/usr/bin/python
# -*- coding: utf-8 -*-

import os
import sys
import math
import argparse

########################################################################################################
def check_stderr( stderr_file=""): 
  stderr=open(stderr_file,'r').readlines()
  no_error= True
  for line in stderr :
     if "ERROR" in line :
        no_error= False
  return no_error
########################################################################################################

########################################################################################################
def check_stdout( stderr_file="",finishmsg=""): 
  stderr=open(stderr_file,'r').readlines()
  finished = False
  no_error= True
  for line in stderr :
     if finishmsg in line :
        finished= True
  return finished
########################################################################################################

########################################################################################################
# parse a str with comma separated ranges: 1,5,10-12
def parse_range(astr):
    result = set()
    for part in astr.split(','):
       x = part.split('-')
       result.update(range(int(x[0]), int(x[-1]) + 1))
    return sorted(result)

########################################################################################################

# MAIN PROGRAM

########################################################################################################
parser = argparse.ArgumentParser(description='# Tool to submit short runs of GVEC\n'
                                             '# example execution command:\n'
                                             '# python gitlab_shortruns.py builddir',\
                                 formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument('builddir',type=str, 
                                 help='Path to builddir from .gitlab.yml folder, mandatory last argument')
parser.add_argument('-execdir',type=str, default='ini',
                                 help='Path to execdir from .gitlab.yml folder, default is ini')
parser.add_argument('-param',type=str, default='../ini/gitlab_shortruns.ini',
                                 help='path to parameterfile to be used for the run, default is ../ini/gitlab_shortruns.ini')
parser.add_argument('-restart',type=str, default='../ini/gitlab_shortruns_restart.ini',
                                 help='path to  parameterfile to be used for the restart run, default is ../ini/gitlab_shortruns_restart.ini')
parser.add_argument('-case', type=str, default='0',    help="0 : DEFAULT, run all cases,\n" 
                                                            "1,2-4 : list of specific cases to run (without spaces!) ")

args = parser.parse_args()

builddir=args.builddir
execdir=args.execdir
param=args.param
restartparam=args.restart

cases = parse_range(args.case)

failed= []
success= []


os.chdir(execdir)
######### SIMULATION ############
caseID=0 

casename="1st_simulation"
print("running caseID %d ,%s ... " % (caseID,casename))
cmd="../"+builddir+"/bin/gvec "+param
print(cmd)
os.system(cmd+" 2>std_"+casename+"_err.txt 1>std_"+casename+"_out.txt")
checkerr = check_stderr("std_"+casename+"_err.txt")
checkout = check_stdout("std_"+casename+"_out.txt","GVEC SUCESSFULLY FINISHED!")


restartFile="GITLAB_RUN_State_0000_99999999.dat"
if ( (not os.path.isfile(restartFile)) or (not checkerr) or (not checkout) ):
  msg=("caseID: %d, %s test did not work!!!" % (caseID,casename))
  print(msg)
  failed.extend([msg])
else:
  msg=("caseID: %d, %s test did execute successfully!" % (caseID,casename))
  print(msg)
  success.extend([msg])

  ######### RESTART ############
  caseID=caseID+1
  if(cases[0]==0 or (caseID in cases)) :
    casename="restart_simulation"
    print("running caseID %d ,%s ... " % (caseID,casename))
    cmd="../"+builddir+"/bin/gvec "+ restartparam + " " + restartFile
    print(cmd)
    os.system(cmd+" 2>std_"+casename+"_err.txt 1>std_"+casename+"_out.txt")
    checkerr = check_stderr("std_"+casename+"_err.txt")
    checkout = check_stdout("std_"+casename+"_out.txt","GVEC SUCESSFULLY FINISHED!")
    
    if ((not os.path.isfile("GITLAB_RESTART_State_0001_99999999.dat")) or (not checkerr) or (not checkout)) :
      msg=("caseID: %d, %s test did not work!!!" % (caseID,casename))
      print(msg)
      failed.extend([msg])
    else :
      msg=("caseID: %d, %s test did execute successfully!" % (caseID,casename) )
      print(msg)
      success.extend([msg])

  ######### GVEC TO HOPR ############
  caseID=caseID+1
  if(cases[0]==0 or (caseID in cases)) :
    casename="gvec_to_hopr"
    print("running caseID %d ,%s ... " % (caseID,casename))
    cmd="../"+builddir+"/bin/test_gvec_to_hopr " + restartFile
    print(cmd)
    os.system(cmd+" 2>std_"+casename+"_err.txt 1>std_"+casename+"_out.txt")
    checkerr = check_stderr("std_"+casename+"_err.txt")
    checkout = check_stdout("std_"+casename+"_out.txt","TEST GVEC TO HOPR FINISHED!")
    if((not checkerr) or (not checkout)) :
      msg=("caseID: %d, %s simulation test did not work!!!" %(caseID,casename))
      print(msg)
      failed.extend([msg])
    else :
      msg=("caseID: %d, %s did execute successfully!"  %(caseID,casename))
      print(msg)
      success.extend([msg])
  
  ######### GVEC TO GENE ############
  caseID=caseID+1
  if(cases[0]==0 or (caseID in cases)) :
    casename="gvec_to_gene"
    print("running caseID %d ,%s ... " % (caseID,casename))
    cmd="../"+builddir+"/bin/test_gvec_to_gene " + restartFile
    print(cmd)
    os.system(cmd+" 2>std_"+casename+"_err.txt 1>std_"+casename+"_out.txt")
    checkerr = check_stderr("std_"+casename+"_err.txt")
    checkout = check_stdout("std_"+casename+"_out.txt","GVEC_TO_GENE FINISHED!")
    if((not checkerr) or (not checkout)) :
      msg=("caseID: %d, %s simulation test did not work!!!" %(caseID,casename))
      print(msg)
      failed.extend([msg])
    else :
      msg=("caseID: %d, %s did execute successfully!"  %(caseID,casename))
      print(msg)
      success.extend([msg])

  ######### GVEC TO CASTOR3D ############
  caseID=caseID+1
  if(cases[0]==0 or (caseID in cases)) :
    casename="convert_gvec_to_castor3d"
    print("running caseID %d ,%s ... " % (caseID,casename))
    cmd="../"+builddir+"/bin/convert_gvec_to_castor3d -r 100 -s 2 " + restartFile + " gvec2castor3d_boozer_output.dat"
    print(cmd)
    os.system(cmd+" 2>std_"+casename+"_err.txt 1>std_"+casename+"_out.txt")
    checkerr = check_stderr("std_"+casename+"_err.txt")
    checkout = check_stdout("std_"+casename+"_out.txt","CONVERT GVEC TO CASTOR3D FINISHED!")
    if ((not os.path.isfile("gvec2castor3d_boozer_output.dat")) or (not checkerr) or (not checkout)) :
      msg=("caseID: %d, %s simulation test did not work!!!" %(caseID,casename))
      print(msg)
      failed.extend([msg])
    else :
      msg=("caseID: %d, %s did execute successfully!"  %(caseID,casename))
      print(msg)
      success.extend([msg])

print( "successful tests:")
for line in success :
   print( "---> "+line )
print( " " )
if(len(failed) > 0 ) :
   print( "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" )
   print( "!!!!!!!    WARNING, following tests  failed:     !!!!!!!!" )
   print( " " )
   for line in failed :
      print( "!!!! ---> "+line )
   print( " " )
   print( "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" )
   sys.exit(100)
else :
   print( "/////////////////////////////////////////////////////////" )
   print( " " )
   print( " ==> ALL TESTS SUCCESSFULL!" )
   print( " " )
   print( "/////////////////////////////////////////////////////////" )
   sys.exit(0)

  
