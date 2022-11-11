#!/usr/bin/env python3
import os,sys
import glob
import htcondor

# Job flavours:
# espresso     = 20 minutes
# microcentury = 1 hour
# longlunch    = 2 hours
# workday      = 8 hours
# tomorrow     = 1 day
# testmatch    = 3 days
# nextweek     = 1 week

directory_to_run = sys.argv[1]
tb_type = sys.argv[2] if len(sys.argv) > 2 else 'October'
if directory_to_run[:4] != 'RUN_':
  directory_to_run = 'RUN_' + directory_to_run
eos_path = f'/eos/project-r/rd51-straw/TestBeam-2022-{tb_type}/data-tiger'

userlogin = os.getlogin()
logdir = os.path.join(f'/afs/cern.ch/user/{userlogin[0]}/{userlogin}/', 'LSFJOBs/tiger/')
if not os.path.exists(logdir):
  os.makedirs(logdir)

col = htcondor.Collector()
credd = htcondor.Credd()
credd.add_user_cred(htcondor.CredTypes.Kerberos, None)

sub = htcondor.Submit()
sub['executable'] = 'convert.sh'
sub['error'] = os.path.join(logdir, f'convert-{tb_type}-{directory_to_run}-$(ClusterId).$(ProcId).err')
sub['output'] = os.path.join(logdir, f'convert-{tb_type}-{directory_to_run}-$(ClusterId).$(ProcId).out')
sub['log'] = os.path.join(logdir, f'convert-{tb_type}-{directory_to_run}-$(ClusterId).log')
sub['MY.SendCredential'] = True
sub['+JobFlavour'] = '"microcentury"'
sub['request_cpus'] = '1'
sub['arguments'] = f'$(inputfile) {directory_to_run} {eos_path}'
sub['+JobBatchName'] = f'"[{tb_type}] Converting tiger folder {directory_to_run}"'

files = [os.path.basename(f) for f in glob.glob(os.path.join(eos_path, f'{directory_to_run}/SubRUN*.dat'))]
itemdata = [{"inputfile": f} for f in files]

print(sub)
print(itemdata)

# from https://htcondor.readthedocs.io/en/latest/apis/python-bindings/tutorials/Submitting-and-Managing-Jobs.html
schedd = htcondor.Schedd()
submit_result = schedd.submit(sub, itemdata = iter(itemdata)) # submit one job for each item in the itemdata

