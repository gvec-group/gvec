#!/usr/bin/python
# -*- coding: utf-8 -*-

import os
import sys
import subprocess
import math

########################################################################################################
def check_stderr( stderr_file=""): 
  stderr=open(stderr_file,'r').readlines()
  no_error= True
  for line in stderr :
     if "error" in line.lower() :
        no_error= False
  return no_error
########################################################################################################

########################################################################################################
def check_stdout( stdout_file="",finishmsg=""): 
  stdout=open(stdout_file,'r').readlines()
  finished = False
  for line in stdout :
     if finishmsg in line :
        finished= True
  return finished
########################################################################################################

########################################################################################################
# use the numdiff tool to compare two files,   using commandline tools "sed" and "numdiff"
#  - ignore_line_ranges: first line ranges can be ignored, in the format "start,end" (end=$ is EOF)
#  - ignore_strings: then lines containing strings can be ignored
#  - parameter 'colcsv': files must be comma-seperated csv files(!) commandline tool "cut" is used
#    one must provide the 'fields' for the cut command. for example, to skip column 2: colcsv='1,3-'
########################################################################################################
def compare_by_numdiff(refpath,reffilename,filepath,filename,reltol="1e-8",abstol="1e-12",ignore_line_ranges=[],ignore_strings=[],colcsv=''):
    msg=[]
    comp=filename.strip()
    p_comp=os.path.join(filepath,comp)
    ref=reffilename.strip()
    p_ref=os.path.join(refpath,ref)
    #check if files exist
    if(not os.path.isfile(p_ref)):
      msg = ('nothing to compare, reference file: "%s" does not exist' % (p_ref))
      return  False , msg
    else:
      #with open(ref,'r') as fp:
      #  reflines = fp.readlines()
      tmpref='tmpref.'+ref
      sbp=subprocess.run(['cp ',p_ref, tmpref],capture_output=True)

        
    if(not os.path.isfile(p_comp)):
      msg = ('nothing to compare, file: "%s" does not exist' % (comp))
      return  False , msg
    else:
      tmpcomp='tmpcomp.'+comp
      sbp=subprocess.run(['cp ',p_comp,tmpcomp],capture_output=True)

    #delete lines  in ranges of ingore_line_ranges (using sed command, format "start,end")
    for ign_lr in ignore_line_ranges:
      sbp=subprocess.run(['cp ',tmpcomp,' tmp'],capture_output=True)
      sbp=subprocess.run(["sed -e '",ign_lr,"d' tmp > ",tmpcomp],capture_output=True) 
      sbp=subprocess.run(['cp ',tmpref,' tmp'],capture_output=True)
      sbp=subprocess.run(["sed -e '",ign_lr,"d' tmp > ",tmpref],capture_output=True)

    
    #delete lines where ingore_strings are found (using sed command)
    for ign_str in ignore_strings:
      sbp=subprocess.run(['cp ',tmpcomp,' tmp'],capture_output=True)
      sbp=subprocess.run(["sed -e '",ign_str,"d' tmp > ",tmpcomp],capture_output=True) 
      sbp=subprocess.run(['cp ',tmpref,' tmp'],capture_output=True)
      sbp=subprocess.run(["sed -e '",ign_str,"d' tmp > ",tmpref],capture_output=True)


    #only compare given columns (works for csv files only, using cut command!)
    if(len(colcsv)>0):
      sbp=subprocess.run(['cp ',tmpcomp,' tmp'],capture_output=True)
      sbp=subprocess.run(['cut -d "," -f',colcsv,' tmp > ',tmpcomp ],capture_output=True)
      sbp=subprocess.run(['cp ',tmpref,' tmp'],capture_output=True)
      sbp=subprocess.run(['cut -d "," -f',colcsv,' tmp > ',tmpref],capture_output=True)

    #compare numerical values with numdiff 
    # first absolute error is checked, 
    #    if abserr < abstol then nodiff
    #    else, relative error is checked (is Inf if one of the values is zero). 
    #        if  relerr < reltol then nodiff
    #        else diff
    sbp_numdiff=subprocess.run(['numdiff -S --separator=" \t\n,;" --relative-tolerance=',reltol,' --absolute-tolerance=',abstol,tmpref,tmpcomp],capture_output=True,text=True)
    sbp= subprocess.run(['rm -f tmp '],capture_output=True)
    
    stdout=sbp_numdiff.stdout
    stderr=sbp_numdiff.stderr
    with open('log_compare_'+ref+'_to_'+comp+'.txt', 'w') as logf:
      logf.write("FILE COMPARISON, INPUT PARAMETERS:\n - abstol=%s \n - reltol=%s\n"%(abstol,reltol))
      for ign_lr in ignore_line_ranges:
        logf.write(' - ignore_line_range="%s"\n' %(ign_lr))
      for ign_str in ignore_strings:
        logf.write(' - ignore_string="%s"\n' %(ign_str))
      if(len(colcsv)>0):
        logf.write(' - colcsv="%s"' %(colcsv))

      logf.write("FILE COMPARISON, EXECUTE COMMAND:\n")
      logf.write(cmd+"\n")
      logf.write("FILE COMPARISON, STDERR FILE:\n")
      for line in stderr :
        logf.write(line)
      logf.write("FILE COMPARISON, STDOUT FILE:\n")
      for line in stdout :
        logf.write(line)
    if(len(stderr)>0):
      nodiff=False
      msg = ('file comparison failed for "%s", check log_compare_* file' %(comp)) 
    else:
      if(len(stdout)>0):
        lastline = stdout[-1]
        nodiff= ('are equal' in lastline.lower())
        if (not nodiff):
          msg= ('file comparison failed for "%s", check log_compare_* file : %s' % (comp,lastline))
      else:
        nodiff=False 
        msg = ('file comparison failed for "%s", check log_compare_* file' %(comp)) 
    if(nodiff):
      msg=('file comparison for "%s" successful! (abstol=%s, reltol=%s)' % (comp,abstol,reltol))
     
    return nodiff, msg
