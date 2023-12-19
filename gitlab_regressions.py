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
parser = argparse.ArgumentParser(description='# Tool to compare regression results of GVEC\n'
                                             '# example execution command:\n'
                                             '# python gitlab_regression.py execdir1 execdir2',\
                                 formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument('execdir1',type=str,
                    help='Name of execdir1 from .gitlab.yml folder, mandatory last argument')
parser.add_argument('execdir2',type=str,
                    help='Name of execdir2 from .gitlab.yml folder, mandatory last argument')
parser.add_argument('-case', type=str, default='0',    help="0 : DEFAULT, run all cases,\n" 
                                                            "1,2-4 : list of specific cases to run (without spaces!) ")

args = parser.parse_args()

execdir1=args.execdir1
execdir2=args.execdir2

cases = parse_range(args.case)

failed= []
success= []


print( '='*132 )
print("RUNNING GITLAB TESTS, execdir1: %s, execdir2: %s"%(execdir1,execdir2))
print( '='*132 )

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
    finalStateFile="GITLAB_RUN_State_0000_00001000.dat"
    finalStateFile1=execdir1+"/"+finalStateFile
    finalStateFile2=execdir2+"/"+finalStateFile
    print("Final state file 1 is: %s" % finalStateFile1)
    print("Final state file 2 is: %s" % finalStateFile2)
    logFile='logMinimizer_GITLAB_RUN_0000.csv'
    logFile1=execdir1+"/"+logFile
    logFile2=execdir2+"/"+logFile
    print("Log file 1 is: %s" % logFile1)
    print("Log file 2 is: %s" % logFile2)
    runfailed=[]
    if(not os.path.isfile(finalStateFile1) or not os.path.isfile(finalStateFile2)):
        msg=("caseID: %d, %s test failed, final state file does not exist" % (caseID,casename))
        failed.extend([msg])
    else:
        #stdout="std_out_"+casename+execdir1+"_vs_"+execdir2+"_statefile.txt"
        #stderr="std_err_"+casename+execdir1+"_vs_"+execdir2+"_statefile.txt"
        #os.system("touch " + stdout)
        #os.system("touch " + stderr)
        #checkerr = check_stderr(stderr)
        nodiff,msg = compare_by_numdiff(execdir1,finalStateFile,execdir2,finalStateFile,ignore_strings=['#'])
        msg=("caseID: %d, %s , %s" %(caseID,casename,msg))
        if (nodiff):
            success.extend([msg])
        else:
            runfailed.extend([msg])
    if(not os.path.isfile(logFile1) or not os.path.isfile(logFile2)):
        msg=("caseID: %d, %s test failed, log file does not exist" % (caseID,casename))
        failed.extend([msg])
    else:
        #stdout="std_out_"+casename+execdir1+execdir2+"_logfile.txt"
        #stderr="std_err_"+casename+execdir1+execdir2+"_logfile.txt"
        nodiff,msg = compare_by_numdiff(execdir1,logFile,execdir2,logFile,colcsv='1,3-')
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

  
