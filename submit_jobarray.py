#!/usr/bin/python
import os,sys,getopt,random,math#,fileinput

#def main(argv):
#    '''Main function that takes in all commandline arguments and parses them to a correct bsub, with output and error files named aptly'''
argv = sys.argv[1:]
mem = '2G'
cores = '1'
name = 'tmp'
wait = 'NULL'
printflag = False
suffix = int(math.floor(random.random()*100000))
suffixflag = True
queue = 'normal'
absPath=os.path.abspath('.')

#Error catcher and definder of option argument
try:
    opts, args = getopt.getopt(argv,"hspr:m:n:w:q:",["room=","mem=","name=","wait=","queue="])
except getopt.GetoptError:
    print 'submit_jobarray.py command-to-be-submittet-from-stdin -[psh] -[rmnwq] <values>'
    sys.exit(2)
#Actually 
for opt, arg in opts:
    if opt == '-h':
        print 'submit_jobarray.py command-to-be-submittet -[psh] -[rmnwq] <values>'
        print '<values>'
        print '-m/--mem: memory required. Default 1g (accepts g as a short for giga, ie. 1g=1000)'
        print '-n/--name: name of job. Default tmp'
        print '-c/--cores: number of cores. Default 1'
        print '-w/--wait: wait on another job, ie dependicies. Default is NULL (add the ID or jobname you want to be dependent on'
        print 'Flags:'
        print '-p: Do you want to bsub command to be printet? Just add the -f flag. Default False, ie no flag)'
        print '-s: Do you want to remove the random number as a suffix to your filenames to increase the probablity of overwriting, appending and all that jazz? Add the flag'
        print '-h prints this help screen and quits'
        sys.exit()
    elif opt in ("-m", "--mem"):
        mem = arg
    elif opt in ("-n", "--name"):
        name = arg
    elif opt in ("-w", "--wait"):
        wait = arg
    elif opt in ("-c", "--cores"):
        cores = arg
    elif opt in ("-p"):
        printflag = True
    elif opt in ("-s"):
        suffixflag = False

jobname = name #maintaining that the jobname is always without the random number

if(suffixflag):
    #adding the random number suffix to the name
    name = name + str(suffix)

#if __name__ == "__main__":
#main(sys.argv[1:])

#print(len(sys.stdin))


   #Writes the command file, adding a line for each command in the case of an array
script = open(name+'.command','a')
maxi=0
for command in sys.stdin:
    script.write(command)
    maxi+=1
maxi = str(maxi)

script.close()
os.system('chmod +x ' + name + '.command')
    
   #Writes the exe file with the simple loop that goes though the command file and sends the command to the grid enginge
exe = open(name+'.exe','w')
exe.write('$(head -n $SGE_TASK_ID '+name+'.command'+' | tail -n 1)')
exe.close()
os.system('chmod +x ' + name + '.exe')

#The bsubcommand with flags
bsubcommand = 'qsub' + \
    ' -S /bin/bash' + \
    ' -cwd' + \
    ' -l h_vmem=' + mem + \
    ' -pe smp ' + cores + \
    ' -N "' + jobname + '"' + \
    ' -o ' + name + '.o ' + \
    ' -e ' + name + '.e ' + \
    ' -t 1-' + maxi
   #Is there any dependices?
if(wait!='NULL'):
    bsubcommand = bsubcommand + ' -hold_jid '+wait
        
   #Finally adds the exe command to the bsubcommand and runs if, execpt if the -p/--print flag is sat    
bsubcommand = bsubcommand + ' "./' + name + '.exe"'
        
if(printflag):
    print bsubcommand
else:
    os.system(bsubcommand)
    print('\n with the following command:')
    print bsubcommand
