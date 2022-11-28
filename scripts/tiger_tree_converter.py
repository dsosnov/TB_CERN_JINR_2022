import ROOT
ROOT.gROOT.SetBatch(True)

import os,sys

folder = sys.argv[1]

py_path = os.path.dirname(__file__)
cxx_path = os.path.join(py_path, 'tiger_tree_converter.cxx')

# ROOT.gROOT.LoadMacro('tiger_tree_converter.cxx');
ROOT.gROOT.LoadMacro(os.path.realpath(cxx_path));
# ROOT.gInterpreter
ROOT.gROOT.ProcessLine('tiger_tree_converter("{}")'.format(folder));
