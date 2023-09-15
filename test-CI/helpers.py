#!/usr/bin/python
# -*- coding: utf-8 -*-

import os
import sys
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
      cmd='cp '+p_ref+' '+ tmpref
      os.system(cmd)
        
    if(not os.path.isfile(p_comp)):
      msg = ('nothing to compare, file: "%s" does not exist' % (comp))
      return  False , msg
    else:
      tmpcomp='tmpcomp.'+comp
      cmd='cp '+p_comp+' '+tmpcomp
      os.system(cmd)

    #delete lines  in ranges of ingore_line_ranges (using sed command, format "start,end")
    for ign_lr in ignore_line_ranges:
      cmd='cp '+tmpcomp+' tmp'
      os.system(cmd)
      cmd="sed -e '"+ign_lr+"d' tmp > "+tmpcomp  
      os.system(cmd)
      cmd='cp '+tmpref+' tmp'
      os.system(cmd)
      cmd="sed -e '"+ign_lr+"d' tmp > "+tmpref
      os.system(cmd)
    
    #delete lines where ingore_strings are found (using sed command)
    for ign_str in ignore_strings:
      cmd='cp '+tmpcomp+' tmp'
      os.system(cmd)
      cmd="sed '/"+ign_str+"/d' tmp > "+tmpcomp  
      os.system(cmd)
      cmd='cp '+tmpref+' tmp'
      os.system(cmd)
      cmd="sed '/"+ign_str+"/d' tmp > "+tmpref  
      os.system(cmd)

    #only compare given columns (works for csv files only, using cut command!)
    if(len(colcsv)>0):
      cmd='cp '+tmpcomp+' tmp'
      os.system(cmd)
      cmd='cut -d "," -f'+colcsv+' tmp > '+tmpcomp  
      os.system(cmd)
      cmd='cp '+tmpref+' tmp'
      os.system(cmd)
      cmd='cut -d "," -f'+colcsv+' tmp > '+tmpref 
      os.system(cmd)

    #compare numerical values with numdiff 
    # first absolute error is checked, 
    #    if abserr < abstol then nodiff
    #    else, relative error is checked (is Inf if one of the values is zero). 
    #        if  relerr < reltol then nodiff
    #        else diff
    cmd='numdiff -S --separator=" \t\n,;" --relative-tolerance='+reltol+' --absolute-tolerance='+abstol+' '+tmpref+' '+tmpcomp
    os.system(cmd+ ' 2>std.err 1>std.out')
    os.system('rm -f tmp ') #+tmpcomp+' '+tmpref )
    
    stdout=open('std.out','r').readlines()
    stderr=open('std.err','r').readlines()
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
