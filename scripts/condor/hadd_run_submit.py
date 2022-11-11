#!/usr/bin/env python3
import os,sys
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

userlogin = os.getlogin()
logdir = os.path.join(f'/afs/cern.ch/user/{userlogin[0]}/{userlogin}/', 'LSFJOBs/tiger/')
if not os.path.exists(logdir):
  os.makedirs(logdir)

outdir = f'/eos/home-{userlogin[0]}/{userlogin}/straw/TestBeam-2022-{tb_type}/'

col = htcondor.Collector()
credd = htcondor.Credd()
credd.add_user_cred(htcondor.CredTypes.Kerberos, None)

sub = htcondor.Submit()
sub['executable'] = 'hadd_run.sh'
sub['error'] = os.path.join(logdir, f'hadd-{tb_type}-{directory_to_run}-$(ClusterId).$(ProcId).err')
sub['output'] = os.path.join(logdir, f'hadd-{tb_type}-{directory_to_run}-$(ClusterId).$(ProcId).out')
sub['log'] = os.path.join(logdir, f'hadd-{tb_type}-{directory_to_run}-$(ClusterId).log')
sub['MY.SendCredential'] = True
sub['+JobFlavour'] = '"microcentury"'
sub['request_cpus'] = '1'
sub['arguments'] = f'{directory_to_run} {outdir}'
sub['+JobBatchName'] = f'"[{tb_type}] hadd tiger run folder {directory_to_run}"'

print(sub)

schedd = htcondor.Schedd()
with schedd.transaction() as txn:
  cluster_id = sub.queue(txn)
