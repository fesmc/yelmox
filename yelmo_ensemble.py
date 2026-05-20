#!/usr/bin/python
# -*- coding: utf-8 -*-
'''
Script to read and modify param namelist file.
'''
import subprocess as subp
import sys, getopt, os, argparse, shutil, glob, datetime, json

# Define parameter classes
class parameter:
    '''
    class to handle parameter input/output/modification
    (only valid for namelist parameter files now!)
    '''

    def __init__(self,group,name="",units="",value="",line="",comment="",module=""):

        # Get parameter info from line, then from arguments
        if not line == "":
            tmp = self.parseline_nml(line)
            name    = tmp['name']
            value   = tmp['value']
            units   = tmp['units']
            comment = tmp['comment']

        # Assign the local values into the class member
        self.name    = name
        self.units   = units
        self.value   = value
        self.line    = line
        self.comment = comment
        self.group   = group
        self.module  = module
        self.mod     = False

        return


    def parseline_nml(self,line):
        '''Parse line from a namelist file.'''

        tmp     = line.split("!")
        par     = tmp[0]

        comment = ""
        if len(tmp)>1: comment = tmp[1].strip()

        tmp     = par.split("=")
        name    = tmp[0].strip()
        value   = tmp[1].strip()

        # Check if comment contains units
        units = ""
        if "[" in comment:
            tmp   = comment.split("[")[1]
            units = tmp.split("]")[0]

        if not comment == "": comment = "! " + comment

        return dict(name=name,value=value,units=units,comment=comment)


    def __repr__(self):
        return repr((self.line,self.name,self.value,self.units,self.group,
                     self.comment,self.mod))

    def short(self):
        '''Output short string representation of parameter and value.
           Used for automatic folder name generation.'''

        # Store the param value as a string
        # Remove the plus sign in front of exponent
        # Remove directory slashes, periods and trailing .nc from string values
        value = "%s" % (self.value)
        if "+" in value: value = value.replace('+','')

        if "/" in value: value = value.replace('/','')
        if ".." in value: value = value.replace('..','')
        if ".nc" in value: value = value.replace('.nc','')

        # Remove all vowels and underscores from parameter name
        name = self.name
        for letter in ['a','e','i','o','u','A','E','I','O','U','_']:
            name = name[0] + name[1:].replace(letter, '')

        return ".".join([name,value])

    def __str__(self):
        '''Output a string suitable for a namelist parameter file'''

        if type(self.value is str):
            return "{} = {} {}".format(self.name,self.value,self.comment)
        else:
            return "{} = {:<9} {}".format(self.name,self.value,self.comment)

class parameters(parameter):
    '''
    This class encapsulates all other classes, and handles
    input and output of entire groups of parameters.
    parameters has parameter
    '''

    def __init__(self,file="",comment="!"):
        '''Initialize instance by loading file containing parameters.'''

        # Store the file name
        self.file = file

        # Initialization
        self.all = []

        if self.file == "":
            # Generate an empty parameter set of one empty parameter
            X = parameter(group="")
            self.all.append(X)

        else:
            # Load all parameters from the input file

            # Make sure file exists, otherwise generate one
            try:
                self.lines = open(self.file,'r').readlines()
            except:
                print "File could not be opened: "+self.file+'\n'
                raise

            # Loop to find parameters and load them into class
            if ".nml" in self.file:

                # Loop through lines and determine which parts correspond to
                # parameters, store these parts in self.all
                inGroup = False

                self.groups = []

                for line in self.lines:
                    line1 = line.strip()
                    if not len(line1)==0 and not line1[0]=="!":

                        if  inGroup and line1[0] == "/":
                            inGroup = False
                            pass

                        if not inGroup and line1[0] == "&":
                            group   = line1.split()[0].strip("&")
                            if not group in self.groups: self.groups.append(group)
                            inGroup = True
                            pass

                        if inGroup and "=" in line1:
                            X = parameter(line=line1,group=group)
                            self.all.append(X)

            else:
                print("Filetype not handled: "+self.file+'\n')
                sys.exit()


        return

    def __str__(self):
        '''Output list of parameters'''

        ll = []
        for p in self.all: ll.append(p.__str__())

        return "\n\n".join(ll)

    def printer(self,mod=False):

        for p in self.all:
            if not mod:
                print p.__str__()
            elif p.mod:
                print p.__str__()

        return

    def write(self,file):
        '''Write a file of all lines in parameter set'''

        try:
            newfile = open(file,'w')

            for group in self.groups:
                newfile.write("&"+group+"\n")
                for p in self.all:
                    if p.group == group:
                        newfile.write(p.__str__()+"\n")
                newfile.write("/\n\n")

            newfile.close()
            print "Parameter file written: %s" % (file)

        except:
            print "Error::p.write: %s not written" % (file)
            sys.exit(2)

        return

    def err(self,exit=2,message="unknown"):
        '''Print an error'''

        print "Error in %s: %s" % (self.file,message)
        print ""
        if exit >= 0: sys.exit(exit)

        return

    def get(self,name):
        '''Return a specific parameter based on the name'''

        for p in self.all:
            if p.name == name:
                return p

        self.err(2,message=name+" not found.")
        return

    def set(self,name="",value="",group="",param=None):
        '''Define the value of a given parameter'''

        # If a parameter is given as argument, take name and value from it
        if not param is None:
            name  = param.name
            value = param.value
            group = param.group

        # Loop over all parameters in set and find the correct one,
        # store new value in its place
        for p in self.all:
            if p.name == name:
                p.value = value; p.group = group; p.mod = True
                return p

        self.err(2,message=name+" not found.")
        return

    def exchange(self,pset=1,file="param_exchange.txt",comment="#",sep=":"):
        '''Exchange parameter values between modules (converting them as needed)'''

        # Open the parameter exchange file
        try:
            flines = open(file,'r').readlines()
        except:
            print "Error: unable to open parameter exchange file: "+file+"\n"
            sys.exit(2)

        #Fill in lists for conversion
        p1 = []; p2 = []; convert = []

        for line in flines:
            if len(line.strip()) > 0:
                if not line.strip()[0] == comment and sep in line:
                    l0 = line.split("=")
                    convert.append(l0[1].strip())
                    l1 = l0[0].split(sep)
                    p1.append(l1[0].strip())
                    p2.append(l1[1].strip())

        # Loop over all exchanges to be made, exchange values!
        for k in range(len(p1)):
            #print "%s : %s  %s" % (p1[k],p2[k],convert[k])

            valnow = pset.get(p1[k]).value

            if not convert[k] == "1":
                con = float(convert[k])
                val = float(valnow)
                valnow = "%g" % (con*val)
                #valnow = valnow.replace("+","")

            pnow = self.set(name=p2[k],value=valnow)

            print pnow

        return

def autofolder(params,outfldr0):
    '''Given a list of parameters,
       generate an appropriate folder name.
    '''

    parts = []

    for p in params:
        parts.append( p.short() )

    # Join the parts together, combine with the base output dir
    autofldr = '.'.join(parts)
    outfldr  = outfldr0 + autofldr + '/'

    return outfldr

def makedirs(dirname):
    '''
    Make a directory (including sub-directories),
    but first ensuring that path doesn't already exist
    or some other error prevents the creation.
    '''

    try:
        os.makedirs(dirname)
        print     'Directory created: ', dirname
    except OSError:
        if os.path.isdir(dirname):
            print 'Directory already exists: ', dirname
            pass
        else:
            # There was an error on creation, so make sure we know about it
            raise

    return

def combiner(a):
    '''a = [[1,2],[3,4,5],[6],[7,8,9,10]]
    '''
    r=[[]]
    for x in a:
        t = []
        for y in x:
            for i in r:
                t.append(i+[y])
        r = t

    return(r)

def load_lhs(filename):
    '''Load Latin-Hypercube sampled table from file.'''

    with open(filename, 'r') as f:
        data = f.readlines()

    values = []
    for line in data:
        values.append(line.split())

    return(values)

def parse_args(args=[],force=False,lhs=None):
    '''Loop over provided arguments and separate them into parameter names and values.
       Default is to assume they are 'rembo' parameters.

       eg, "melt_choice=1 pdd_factor=1"
       eg, "rembo="melt_choice=1 pdd_factor=1" sico="dtime_ser=50"
    '''

    params = []; group = []

    # Sort the arguments
    args.sort()

    # Loop over the arguments and parse them into separate parameters
    if 'rembo=' in " ".join(args) or 'sico=' in " ".join(args) \
      or 'climber=' in " ".join(args) or '&' in " ".join(args):

        for arg in args:
            now = "none"
            if 'rembo='   in arg: now = "rembo"
            if 'sico='    in arg: now = "sico"
            if 'climber=' in arg: now = "climber"
            if 'grisli='  in arg: now = "grisli"
            if '&'        in arg: now = "&"

            if not now == "none":
                tmp = arg.partition("=")
                now = tmp[0].strip("&")
                tmp = tmp[2].split()
                tmp.sort()

                for p in tmp:
                    params.append(p)
                    group.append(now)
            else:
                print "\nArguments should be encapsulated by the group name,"
                print 'eg, rembo="pname1=5" sico="pname2=10"\n'
                sys.exit(2)

    else: # only the default module parameters will be loaded
        for arg in args:
            params.append(arg)
            group.append("none")

    # Check values of arguments
    print ""
    k = 0
    for p in params:
        m = group[k]; k = k + 1
        print m + ":" + p

    if not lhs is None:
        print "**LHS sampling will be loaded from: "+lhs
        print "  to produce the above parameter ranges."

    if not force:
        try:
            response = raw_input("\n[Enter to proceed] or [ctl-c to exit]")
            print "\n"
        except:
            print "\n"
            sys.exit()


    names = []; values = []

    for p in params:
        if "=" not in p:
            print 'Error::parse_args: parameter names and values must be separated by "="'
            print ''
            sys.exit(2)

        # Separate term into name and value(s)
        tmp = p.split("=")

        name = tmp[0]
        vals = [tmp[1]]
        if "," in tmp[1]: vals = tmp[1].split(",")

        names.append(name)
        values.append(vals)

    # Set up an initial list to store all parameter sets (1 set per run)
    batch = []

    if lhs is None:
        # Make all permutations of parameter options
        allvalues = combiner(values)
    else:
        # Load parameter combinations from a saved table
        allvalues = load_lhs(lhs)

        # Make sure number of parameters == number of values loaded per combination
        if not len(names) == len(allvalues[0]):
            print 'Error::parse_args: number of LHS parameters loaded does not match command line options.'
            print ''
            sys.exit(2)

        # Multiply lhs values with parameter ranges
        ns = len(allvalues)
        np = len(allvalues[0])
        for s in range(ns):
            for k in range(np):
                allvalues[s][k] = '%.3f' % ( (float(values[k][1])-float(values[k][0]))*float(allvalues[s][k])+float(values[k][0]) )

    # Loop over each set of parameter values
    for values in allvalues:

        # Initialize the current new set of parameters
        set = []

        for k in range(len(names)):

            # Generate a new parameter and store in the set
            X = parameter(name=names[k],value=values[k],group=group[k])
            set.append(X)

        # Save the set inside the batch list
        batch.append(set)

    ## Check what has been created ##
##    for b in batch:
##        print "Batch set ====="
##        for p in b:
##            print p

    return batch

def makejob(params,out,wtime,executable,suffix="",auto=False,force=False,edit=False,submit=False,case="none"):
    '''Given a set of parameters, generate output folder and
       set up a job, then submit it.
    '''
    
    infldr = "./" # Path to original parameter files

    # Determine the base name of the parameter files (no folders)
    #o1_in  = "par/yelmo_ismip6_Antarctica_spinup.nml"
    #o1_out = "yelmo_ismip6_Antarctica_spinup.nml"
    o1_in  = "par/yelmo_ismip6_Antarctica.nml"
    o1_out = "yelmo_ismip6_Antarctica.nml"

    # Generate a name for the output folder, if desired
    if auto:
        outfldr = autofolder(params,out)
    else:
        outfldr = out
        
    # First generate some empty data sets for each module
    p_yelmo   = parameters()

    # Now load the relevant default parameter files
    p_yelmo = parameters(file=os.path.join(infldr,o1_in))

    # Store parameter lists as new objects for modification and output
    p_yelmo1  = p_yelmo

    # Loop over command line parameters and modify parameter lists
    for p in params:
        p_yelmo1.set(param=p)

    # Output a separator
    print "============================================================"

    # Print info to screen
    print "\nOutput directory: " + outfldr
    print "\nModified parameters: Yelmo"
    p_yelmo1.printer(mod=True)

    # Make output directory and subdirectory
    makedirs(outfldr)

    # Write new parameter file to output directory
    p_yelmo1.write(os.path.join(outfldr,o1_out))

    ## Copy files of interest

    # Generate the appropriate executable command to run job
    # Assume executable is running from outlfdr
    #executable = "./{}".format(exe_fname)

    # Get path of constants parameter file and copy it 
    const_path = "par/yelmo_const_Earth.nml"
    shutil.copy(const_path,outfldr)
    # Get path of ismip6 parameter file and copy it
    ismip6_path = "par/ismip6.nml"
    shutil.copy(ismip6_path,outfldr)

    # Also copy exe file to rundir
    shutil.copy("libyelmox/bin/"+executable,outfldr)
        
    ## Control parameter files 
    # shutil.copy(os.path.join(infldr,o2),os.path.join(outfldr,o2))

    ## Generate symbolic links to input data folders
    srcname = "input"
    dstname = os.path.join(outfldr,srcname)
    srcpath = os.path.abspath(srcname)
    if os.path.islink(dstname): os.unlink(dstname)
    os.symlink(srcpath,dstname)

    srcname = "ice_data"
    dstname = os.path.join(outfldr,srcname)
    if os.path.islink(dstname): os.unlink(dstname)
    if os.path.islink(srcname):
        linkto = os.readlink(srcname)
        os.symlink(linkto, dstname)
    elif os.path.isdir(srcname):
        srcpath = os.path.abspath(srcname)
        os.symlink(srcpath,dstname)
    else:
        print("Warning: path does not exist {}".format(srcname))

    srcname = "ice_data_javi"
    dstname = os.path.join(outfldr,srcname)
    if os.path.islink(dstname): os.unlink(dstname)
    if os.path.islink(srcname):
        linkto = os.readlink(srcname)
        os.symlink(linkto, dstname)
    elif os.path.isdir(srcname):
        srcpath = os.path.abspath(srcname)
        os.symlink(srcpath,dstname)
    else:
        print("Warning: path does not exist {}".format(srcname))

    # Copy the submit.sh file
    submit_file = "submit.sh"
    shutil.copy(submit_file,outfldr)

    # 2. Run the job

    # Generate the appropriate executable command to run job
    # Assume executable is running from rundir
    #executable = "./{}".format(exe_fname)

    # Also copy exe file to rundir 
    #shutil.copy(exe_path,outfldr)

    # Submit job to queue 
    pid = submitjob(outfldr,executable)

    return outfldr


def submitjob(outfldr,executable):
    '''Submit a job to a HPC queue (qsub,sbatch)'''

    # Get info about current system
    # username  = os.environ.get('USER')
    # hostname  = os.environ.get('HOSTNAME')

    # Command to be called
    #cmd = "{} {}".format(executable,par_path)

    # Create the jobscript using current info
    #nm_jobscript   = 'job.submit'
    #path_jobscript = "{}/{}".format(rundir,nm_jobscript)

    # Jobscript for qsub
    #script = jobscript_qsub(cmd,rundir,username,usergroup,wtime,useremail)
    #jobfile = open(path_jobscript,'w').write(script)
    cmd_job = "cd {} && sbatch submit.sh".format(outfldr)

    # Run the command (ie, change to output directory and submit job)
    # Note: the argument `shell=True` can be a security hazard, but should
    # be ok in this context, see https://docs.python.org/2/library/subprocess.html#frequently-used-arguments
    jobstatus = subp.check_call(cmd_job,shell=True)

    return jobstatus

def main():

    # Default values of options #
    executable = 'yelmox_ismip6.x'        # Exectutable program (default: Ant-40)
    submit     = False           # Submit the job to loadleveler (default: no)
    edit       = False           # Interactive editing of options
    outfldr    = 'output/test/'  # Default output folder is just the outbase
    auto       = False
    force      = False
    wtime      = "24"            # Default wall clock time is 24 hours
    case       = "none"          # No special simulation case is being run
    lhs        = None            # No Latin-Hypercube sampling file given

    pname_suffix = ""            # Suffix to append to default parameter filename

    # Get a list of options and arguments
    try:
        opts, args = getopt.getopt(sys.argv[1:], "hlep:o:a:fw:t:s:", ["help", "program=","edit=","out=","auto=","wall=","lhs=","suffix="])
    except getopt.GetoptError, err:
        # print help information and exit:
        usage()
        print "\n    ",str(err) # will print something like "option -a not recognized"
        sys.exit(2)

    # Cycle through options
    for o, a in opts:
        if o in ("-h", "--help"):
            usage()
            sys.exit()
        elif o in ("-p", "--program"):
            executable = a              # Executable program is input argument
        elif o in ("-l"):
            submit = True               # Job will be submitted to loadleveler
        elif o in ("-e", "--edit"):
            edit   = True               # User will interactively edit options
        elif o in ("-f"):
            force   = True              # Jobs will start without user confirmation
        elif o in ("-o", "--out"):
            outfldr = a + '/'
        elif o in ("-a", "--auto"):
            auto = True                 # No automatic folder generation, use user choice
            outfldr = a + '/'
        elif o in ("-w", "--wall"):
            wtime = a
        elif o in ("-s", "--suffix"):
            pname_suffix = a
        elif o in ("-t"):
            case = a
        elif o in ("--lhs"):
            lhs = a
        else:
            assert False, "unhandled option"

    # Get the batch parameter sets from the arguments
    # (returns an empty set if no parameters should be changed)
    batch = parse_args(args,force=force,lhs=lhs)

    # Make sure that if generating multiple runs
    # that the --auto option has been used
    if len(batch) > 1 and not auto:
        print "\nError: automatic folder generation must be used for batch processing!\n"
        sys.exit(2)

    # Loop over the parameter sets and make jobs
    joblist = []
    for params in batch:

        # Go through job setup and submit process
        fldr = makejob(params,outfldr,wtime,executable,pname_suffix,auto,force,edit,submit,case)

        joblist.append(fldr)

    # Write the job list to a file
    # (make the output folder relative to the output/ directory)
    joblist  = "\n".join(joblist)
    joblist1 = joblist.replace("output/","")
    joblist1 = joblist1.replace("outtmp/","")

    try:
        if os.path.isfile(outfldr+"batch"):
            open(outfldr+"batch","a").write("\n"+joblist1)
        else:
            open(outfldr+"batch","w").write(joblist1)

        print "Output folder(s):\n"
        print joblist
        print "\n"

    except:
        print "Unable to write batch list to " + outfldr

    return

main()
