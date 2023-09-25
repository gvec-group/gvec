#!/usr/bin/python
# -*- coding: utf-8 -*-

import os
import sys
import math
import argparse
from gitlab_helpers import *

########################################################################################################
# parse a str with comma separated ranges: 1,5,10-12
def parse_range(astr):
    result = set()
    for part in astr.split(','):
       x = part.split('-')
       result.update(range(int(x[0]), int(x[-1]) + 1))
    return sorted(result)
########################################################################################################


########################################################################################################

# MAIN PROGRAM

########################################################################################################
parser = argparse.ArgumentParser(description='# Tool to submit short runs of GVEC\n'
                                             '# example execution command:\n'
                                             '# python gitlab_shortruns.py builddir',\
                                 formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument('builddir',type=str, 
                                 help='Name of builddir from .gitlab.yml folder, mandatory last argument')
parser.add_argument('-execdir',type=str, default='test-CI',
                                 help='Name of execdir from .gitlab.yml folder, default is test-CI')
parser.add_argument('-execpre',type=str, default='',
                                 help='command prefix to execute (like "mpirun -np 2", or "export OMP_NUM_THREADS=10;")')
parser.add_argument('-case', type=str, default='0',    help="0 : DEFAULT, run all cases,\n" 
                                                            "1,2-4 : list of specific cases to run (without spaces!) ")

args = parser.parse_args()

builddir=args.builddir
execdir=args.execdir
execpre=args.execpre

cases = parse_range(args.case)

failed= []
success= []


print( '='*132 )
print("RUNNING GITLAB TESTS, builddir: %s, execdir: %s"%(builddir,execdir))
print( '='*132 )

os.chdir(execdir)

refdir="../test-CI"
referenceStateFile="REF_GITLAB_RUN_State_0000_00001000.dat"
restartFile=os.path.join(refdir,referenceStateFile)

#########################################################################
### GVEC SIMULATION 
#########################################################################
caseID=1 
if(cases[0]==0 or (caseID in cases)) :
  casename="1st_simulation"
  print("running caseID %d ,%s ... " % (caseID,casename))
  execbin="../"+builddir+"/bin/gvec"
  param="../test-CI/gvec_shortrun.ini"
  finalStateFile="GITLAB_RUN_State_0000_00001000.dat"
  logFile='logMinimizer_GITLAB_RUN_0000.csv'

  if(not os.path.isfile(execbin)):
    msg=("caseID: %d, %s test failed, executable does not exist" % (caseID,casename))
    failed.extend([msg])
  else:
    runfailed=[]
    cmd=execpre+" "+execbin+" "+param
    print(cmd)
    stdout="std_out_"+casename+".txt"
    stderr="std_err_"+casename+".txt"
    os.system(cmd+" 2>"+stderr+" 1>"+stdout)
    checkerr = check_stderr(stderr)
    if(not checkerr):
      msg=("caseID: %d, %s test failed, problem in stderr !!!" % (caseID,casename))
      runfailed.extend([msg])
    checkout = check_stdout(stdout,"GVEC SUCESSFULLY FINISHED!")
    if(not checkout):
      msg=("caseID: %d, %s test failed, problem in stdout !!!" % (caseID,casename))
      runfailed.extend([msg])
    nodiff,msg = compare_by_numdiff(refdir,referenceStateFile,'',finalStateFile,ignore_strings=['#'])
    msg=("caseID: %d, %s , %s" %(caseID,casename,msg))
    if (nodiff):
      success.extend([msg])
    else:
      runfailed.extend([msg])
    nodiff,msg = compare_by_numdiff(refdir,'REF_'+logFile,'',logFile,colcsv='1,3-')
    msg=("caseID: %d, %s , %s" %(caseID,casename,msg))
    if (nodiff):
      success.extend([msg])
    else:
      runfailed.extend([msg])
    if(len(runfailed) > 0) :
      for line in runfailed :
         print( "!!!! ---> "+line )
      msg=("caseID: %d, %s test failed!"  %(caseID,casename))
      failed.extend(runfailed)
    else:
      msg=("caseID: %d, %s did execute successfully!"  %(caseID,casename))
      success.extend([msg])
  print(msg)

#########################################################################
### RESTART
#########################################################################
caseID=caseID+1
if(cases[0]==0 or (caseID in cases)) :
  casename="restart_simulation"
  print("running caseID %d ,%s ... " % (caseID,casename))
  execbin="../"+builddir+"/bin/gvec"
  param="../test-CI/gvec_shortrun_restart.ini"
  finalStateFile="GITLAB_RESTART_State_0001_00000010.dat"
  logFile='logMinimizer_GITLAB_RESTART_0001.csv'
  if(not os.path.isfile(execbin)):
    msg=("caseID: %d, %s test failed, executable does not exist" % (caseID,casename))
    failed.extend([msg])
  else:
    runfailed=[]
    cmd=execpre+" ""../"+builddir+"/bin/gvec "+ param + " " + restartFile
    print(cmd)
    stdout="std_out_"+casename+".txt"
    stderr="std_err_"+casename+".txt"
    os.system(cmd+" 2>"+stderr+" 1>"+stdout)
    checkerr = check_stderr(stderr)
    if(not checkerr):
      msg=("caseID: %d, %s test failed, problem in stderr !!!" % (caseID,casename))
      runfailed.extend([msg])
    checkout = check_stdout(stdout,"GVEC SUCESSFULLY FINISHED!")
    if(not checkout):
      msg=("caseID: %d, %s test failed, problem in stdout !!!" % (caseID,casename))
      runfailed.extend([msg])
    nodiff,msg = compare_by_numdiff(refdir,'REF_'+finalStateFile,'',finalStateFile,ignore_strings=['#'])
    msg=("caseID: %d, %s , %s" %(caseID,casename,msg))
    if (nodiff):
      success.extend([msg])
    else:
      runfailed.extend([msg])
    nodiff,msg = compare_by_numdiff(refdir,'REF_'+logFile,'',logFile,colcsv='1,3-')
    msg=("caseID: %d, %s , %s" %(caseID,casename,msg))
    if (nodiff):
      success.extend([msg])
    else:
      runfailed.extend([msg])
    if(len(runfailed) > 0) :
      for line in runfailed :
         print( "!!!! ---> "+line )
      msg=("caseID: %d, %s test failed!"  %(caseID,casename))
      failed.extend(runfailed)
    else:
      msg=("caseID: %d, %s did execute successfully!"  %(caseID,casename))
      success.extend([msg])
  print(msg)

#########################################################################
### GVEC TO HOPR
#########################################################################
caseID=caseID+1
if(cases[0]==0 or (caseID in cases)) :
  casename="gvec_to_hopr"
  print("running caseID %d ,%s ... " % (caseID,casename))
  execbin="../"+builddir+"/bin/test_gvec_to_hopr"
  if(not os.path.isfile(execbin)):
    msg=("caseID: %d, %s test failed, executable does not exist" % (caseID,casename))
    failed.extend([msg])
  else:
    runfailed=[]
    cmd=execpre+" "+execbin+" " + restartFile
    print(cmd)
    stdout="std_out_"+casename+".txt"
    stderr="std_err_"+casename+".txt"
    os.system(cmd+" 2>"+stderr+" 1>"+stdout)
    checkerr = check_stderr(stderr)
    checkout = check_stdout(stdout,"TEST GVEC TO HOPR FINISHED!")
    if(not checkerr):
      msg=("caseID: %d, %s test failed, problem in stderr !!!" % (caseID,casename))
      runfailed.extend([msg])
    if(not checkout):
      msg=("caseID: %d, %s test failed, problem in stdout !!!" % (caseID,casename))
      runfailed.extend([msg])
    nodiff,msg = compare_by_numdiff(refdir,'REF_'+stdout,'',stdout,ignore_strings=[' sec'],reltol="2e-08")
    msg=("caseID: %d, %s , %s" %(caseID,casename,msg))
    if (nodiff):
      success.extend([msg])
    else:
      runfailed.extend([msg])
    if(len(runfailed) > 0) :
      for line in runfailed :
         print( "!!!! ---> "+line )
      msg=("caseID: %d, %s test failed!"  %(caseID,casename))
      failed.extend(runfailed)
    else:
      msg=("caseID: %d, %s did execute successfully!"  %(caseID,casename))
      success.extend([msg])
  print(msg)

#########################################################################
### GVEC TO GENE
#########################################################################
caseID=caseID+1
if(cases[0]==0 or (caseID in cases)) :
  casename="gvec_to_gene"
  print("running caseID %d ,%s ... " % (caseID,casename))
  execbin="../"+builddir+"/bin/test_gvec_to_gene"
  if(not os.path.isfile(execbin)):
    msg=("caseID: %d, %s test failed, executable does not exist" % (caseID,casename))
    failed.extend([msg])
  else:
    runfailed=[]
    cmd=execpre+" "+execbin+" " + restartFile
    print(cmd)
    stdout="std_out_"+casename+".txt"
    stderr="std_err_"+casename+".txt"
    os.system(cmd+" 2>"+stderr+" 1>"+stdout)
    checkerr = check_stderr(stderr)
    checkout = check_stdout(stdout,"GVEC_TO_GENE FINISHED!")
    if(not checkerr):
      msg=("caseID: %d, %s test failed, problem in stderr !!!" % (caseID,casename))
      runfailed.extend([msg])
    if(not checkout):
      msg=("caseID: %d, %s test failed, problem in stdout !!!" % (caseID,casename))
      runfailed.extend([msg])
    if(len(runfailed) > 0) :
      for line in runfailed :
         print( "!!!! ---> "+line )
      msg=("caseID: %d, %s test failed!"  %(caseID,casename))
      failed.extend(runfailed)
    else:
      msg=("caseID: %d, %s did execute successfully!"  %(caseID,casename))
      success.extend([msg])
  print(msg)

#########################################################################
### GVEC TO CASTOR3D 
#########################################################################
for sflcoord in ["0","1","2"]:
  caseID=caseID+1
  if(cases[0]==0 or (caseID in cases)) :
    casename="convert_gvec_to_castor3d_sflcoord_"+sflcoord
    print("running caseID %d ,%s ... " % (caseID,casename))
    execbin="../"+builddir+"/bin/convert_gvec_to_castor3d"
    outfile="gvec2castor3d_sfl_"+sflcoord+"_output.dat"
    if(not os.path.isfile(execbin)):
      msg=("caseID: %d, %s test failed, executable does not exist" % (caseID,casename))
      failed.extend([msg])
    else:
      runfailed=[]
      cmd=execpre+" "+execbin+" -r 5 -p 8 -t 4 -s "+sflcoord+" " + restartFile + " " + outfile
      print(cmd)
      stdout="std_out_"+casename+".txt"
      stderr="std_err_"+casename+".txt"
      os.system(cmd+" 2>"+stderr+" 1>"+stdout)
      checkerr = check_stderr(stderr)
      if(not checkerr):
        msg=("caseID: %d, %s test failed, problem in stderr !!!" % (caseID,casename))
        runfailed.extend([msg])
      checkout = check_stdout(stdout,"CONVERT GVEC TO CASTOR3D FINISHED!")
      if(not checkout):
        msg=("caseID: %d, %s test failed, problem in stdout !!!" % (caseID,casename))
        runfailed.extend([msg])
      nodiff,msg = compare_by_numdiff(refdir,'REF_'+outfile,'',outfile,ignore_strings=['#'],reltol="5e-8")
      msg=("caseID: %d, %s , %s" %(caseID,casename,msg))
      if (nodiff):
        success.extend([msg])
      else:
        runfailed.extend([msg])
      if(len(runfailed) > 0) :
        for line in runfailed :
           print( "!!!! ---> "+line )
        msg=("caseID: %d, %s test failed!"  %(caseID,casename))
        failed.extend(runfailed)
      else:
        msg=("caseID: %d, %s did execute successfully!"  %(caseID,casename))
        success.extend([msg])
    print(msg)

#########################################################################
### GVEC TO CASTOR3D 
#########################################################################
caseID=caseID+1
if(cases[0]==0 or (caseID in cases)) :
  casename="convert_gvec_to_jorek"
  print("running caseID %d ,%s ... " % (caseID,casename))
  execbin="../"+builddir+"/bin/convert_gvec_to_jorek"
  outfile="gvec2jorek_output.dat"
  if(not os.path.isfile(execbin)):
    msg=("caseID: %d, %s test failed, executable does not exist" % (caseID,casename))
    failed.extend([msg])
  else:
    runfailed=[]
    cmd=execpre+" "+execbin+" -r 4 -p 8 " + restartFile + " "+outfile
    print(cmd)
    stdout="std_out_"+casename+".txt"
    stderr="std_err_"+casename+".txt"
    os.system(cmd+" 2>"+stderr+" 1>"+stdout)
    checkerr = check_stderr(stderr)
    if(not checkerr):
      msg=("caseID: %d, %s test failed, problem in stderr !!!" % (caseID,casename))
      runfailed.extend([msg])
    checkout = check_stdout(stdout,"CONVERT GVEC TO JOREK FINISHED!")
    if(not checkout):
      msg=("caseID: %d, %s test failed, problem in stdout !!!" % (caseID,casename))
      runfailed.extend([msg])
    nodiff,msg = compare_by_numdiff(refdir,'REF_'+outfile,'',outfile,ignore_strings=['#']
                                    ,abstol="1e-8"  # for large numbers,  difference debug/release
                                    ,ignore_line_ranges=["1145,$"]) # ignore the current, computed with FD 
    msg=("caseID: %d, %s , %s" %(caseID,casename,msg))
    if (nodiff):
      success.extend([msg])
    else:
      runfailed.extend([msg])
    if(len(runfailed) > 0) :
      for line in runfailed :
         print( "!!!! ---> "+line )
      msg=("caseID: %d, %s test failed!"  %(caseID,casename))
      failed.extend(runfailed)
    else:
      msg=("caseID: %d, %s did execute successfully!"  %(caseID,casename))
      success.extend([msg])
  print(msg)

#########################################################################
# SUMMARY
#########################################################################


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

  
