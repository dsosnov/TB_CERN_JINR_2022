#!/usr/bin/env python3
import os,sys
import glob
import subprocess
import shutil
import htcondor

# https://stackoverflow.com/a/64049576
def shell_source(script):
  """
  Sometime you want to emulate the action of "source" in bash,
  settings some environment variables. Here is a way to do it.
  """

  pipe = subprocess.Popen(". %s && env -0" % script, stdout=subprocess.PIPE, shell=True)
  output = pipe.communicate()[0].decode('utf-8')
  output = output[:-1] # fix for index out for range in 'env[ line[0] ] = line[1]'

  env = {}
  # split using null char
  for line in output.split('\x00'):
    line = line.split( '=', 1)
    # print(line)
    env[ line[0] ] = line[1]

  os.environ.update(env)


mappings = {
  'October': {
    0: 'empty',
    2: '20221020',
    6: '20221024',
    13: '20221025',
    24: '20221026',
    28: '20221027',
    36: '20221028',
    51: '20221031',
    54: '20221031-s',
    57: '20221031-s2',
    59: '20221031-s3',
  },
  'November': {
    0: '20221031', # 'empty',
  }
}  

# ./run_condor.py analyse "RUN_31" "RD.root" "/eos/home-d/dsosnov/SPD/TestBeam-2022-October"
def analyse():
  directory = sys.argv[2]
  filename = sys.argv[3][:-5]
  outdir = sys.argv[4]
  eospath = sys.argv[5]
  fullpath_tb_analysis = os.path.abspath(sys.argv[6])
  tb_type = sys.argv[7]

  directories_to_link = ['configs','scripts']

  for d in directories_to_link:
    if not os.path.exists(d):
      os.symlink(os.path.join(fullpath_tb_analysis, d), d)
  if not os.path.exists('out'):
    os.mkdir('out')
  if not os.path.exists('code'):
    os.mkdir('code')
  os.chdir('code')
  for src in glob.glob(os.path.join(fullpath_tb_analysis, 'code/*.*')):
    shutil.copy(src, './')

  shutil.copy('../process_tiger.cxx', './')
  shell_source('/cvmfs/sft.cern.ch/lcg/views/LCG_102/x86_64-centos7-gcc11-opt/setup.sh')

  os.makedirs('../data/tiger')
  os.symlink(f'{eospath}/{directory}', f'../data/tiger/{directory}')

  run_num = int(directory.split('_')[1])
  mapping = 'empty'
  for k in sorted(mappings[tb_type].keys()):
    if run_num >= k:
      mapping = mappings[tb_type][k]
    else:
      break;

  # root_command = ['root', '-b', '-q', ] #'-e', f'gROOT->ProcessLine(".L tiger.C"); gROOT->ProcessLine("(new tiger(\\"{filename}\\", \\"{directory}\\", \\"{mapping}\\", \\"20221025\\"))->Loop()")']
  root_command = ['root', '-b', '-q', f'process_tiger.cxx("{filename}", "{directory}", "{mapping}", "20221025")']
  print(root_command)
  subprocess.run(root_command);

  shutil.copy(f'../out/out_tiger_{directory}-{filename}.root', outdir)

def submit_dir(d, tb_type):
  userlogin = os.getlogin()
  script_basename = os.path.basename(sys.argv[0])
  fullpath_tb_analysis = os.path.abspath('../../../')
  eos_path = f'/eos/project-r/rd51-straw/TestBeam-2022-{tb_type}/data-tiger/'
  ddir = 'RUN_'+str(d)

  outdir_base = f'/eos/home-{userlogin[0]}/{userlogin}/straw/TestBeam-2022-{tb_type}/'
  if not os.path.exists(outdir_base):
    os.makedirs(outdir_base)
  outdir_abs = os.path.abspath(os.path.join(outdir_base, ddir))
  if not os.path.exists(outdir_abs):
    os.makedirs(outdir_abs)


  logdir = os.path.join(f'/afs/cern.ch/user/{userlogin[0]}/{userlogin}/', 'LSFJOBs/tiger/')
  if not os.path.exists(logdir):
    os.makedirs(logdir)

  col = htcondor.Collector()
  credd = htcondor.Credd()
  credd.add_user_cred(htcondor.CredTypes.Kerberos, None)

  sub = htcondor.Submit()
  sub['executable'] = sys.argv[0]
  sub['error'] = os.path.join(logdir, f'{script_basename}-{tb_type}-{ddir}-$(ClusterId).$(ProcId).err')
  sub['output'] = os.path.join(logdir, f'{script_basename}-{tb_type}-{ddir}-$(ClusterId).$(ProcId).out')
  sub['log'] = os.path.join(logdir, f'{script_basename}-{tb_type}-{ddir}-$(ClusterId).log')
  sub['MY.SendCredential'] = True
  sub['+JobFlavour'] = '"longlunch"'
  sub['request_cpus'] = '1'
  sub['transfer_input_files'] = 'process_tiger.cxx'
  sub['should_transfer_files'] = True
  sub['arguments'] = f'analyse {ddir} $(inputfile) {outdir_abs} {eos_path} {fullpath_tb_analysis} {tb_type}'
  sub['+JobBatchName'] = f'"[{tb_type}] {script_basename}: analysing tiger folder {ddir}"'

  files = [os.path.basename(f) for f in glob.glob(os.path.join(eos_path, f'{ddir}/SubRUN*.root'))]
  itemdata = [{"inputfile": f} for f in files]

  print(sub)
  print(itemdata)

  # from https://htcondor.readthedocs.io/en/latest/apis/python-bindings/tutorials/Submitting-and-Managing-Jobs.html
  schedd = htcondor.Schedd()
  submit_result = schedd.submit(sub, itemdata = iter(itemdata)) # submit one job for each item in the itemdata

if __name__ == '__main__':
  if len(sys.argv) > 1 and sys.argv[1] == 'analyse':
    analyse()
  elif len(sys.argv) > 1 and sys.argv[1] == 'submit':
    tb_type = sys.argv[3] if len(sys.argv) > 3 else 'October'
    submit_dir(sys.argv[2], tb_type)
  else:
      pass
