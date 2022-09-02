import sys, re, os
infile = sys.argv[1]

detectors_map = [{} for i in range(8)]
with open(infile, 'r') as f:
  for line in f.readlines():
    conf_text = line.split('#', 1)[0].strip()
    if not conf_text:
      continue
    conf_splitted = [int(i) for i in re.split('\t| ', conf_text) if i.strip()]
    # print(conf_splitted)
    detectors_map[conf_splitted[3]][conf_splitted[4]] = tuple(conf_splitted[:3])
# print(detectors_map)

detector_map_filename = 'tiger_' + os.path.basename(infile).replace('.txt', '.cxx')
with open(detector_map_filename, 'w') as f:
  for i in range(len(detectors_map)):
    print(f'void tiger_get_det_{i}()' '{' f'printf("tiger_get_det_{i}(gemrocID, tigerID, channelID)\\n");' '}', file=f)
    print(f'int tiger_get_det_{i}(Char_t gemrocID, Short_t tigerID, Char_t channelID)''{', file=f)
    for k, v in detectors_map[i].items():
      print(f'  if(gemrocID=={v[0]} && tigerID=={v[1]} && channelID=={v[2]})''', file=f)
      print(f'    return {k};', file=f)
      print('  else', file=f)
    print('    return -(channelID + tigerID*64 + gemrocID * 64*8);', file=f)
    print('}', file=f)
    print('', file=f)
    pass
  pass
