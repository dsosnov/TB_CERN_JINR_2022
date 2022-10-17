import sys

# connector pin: (tiger, channel), from FEB3_pinout.pdf
tigChMap = {
  1: (1,3),
  2: (1,45),
  3: (1,11),
  4: (1,47),
  5: (1,22),
  6: (1,54),
  7: (1,4),
  8: (1,41),
  9: (1,12),
  10: (1,43),
  11: (1,23),
  12: (1,33),
  13: (1,19),
  14: (1,39),
  15: (1,15),
  16: (1,37),
  17: (1,14),
  18: (1,6),
  19: (1,60),
  20: (1,9),
  21: (1,0),
  22: (1,25),
  23: (1,26),
  24: (1,29),
  25: (1,49),
  26: (1,7),
  27: (1,55),
  28: (1,1),
  29: (1,51),
  30: (1,27),
  31: (1,20),
  32: (1,13),
  33: (1,16),
  34: (1,31),
  35: (1,28),
  36: (1,35),
  37: (1,59),
  38: (1,57),
  39: (1,56),
  40: (1,53),
  41: (1,50),
  42: (1,2),
  43: (1,52),
  44: (1,17),
  45: (1,58),
  46: (1,21),
  47: (1,46),
  48: (1,8),
  49: (1,40),
  50: (1,10),
  51: (1,48),
  52: (1,5),
  53: (1,44),
  54: (1,24),
  55: (1,42),
  56: (1,18),
  57: (1,36),
  58: (1,30),
  59: (1,34),
  60: (1,31),
  61: (1,38),
  62: (0,51),
  63: (0,60),
  64: (0,57),
  65: (0,30),
  66: (0,53),
  67: (0,31),
  68: (0,16),
  69: (0,22),
  70: (0,55),
  71: (0,54),
  72: (0,21),
  73: (0,50),
  74: (0,19),
  75: (0,59),
  76: (0,24),
  77: (0,18),
  78: (0,26),
  79: (0,58),
  80: (0,27),
  81: (0,20),
  82: (0,25),
  83: (0,29),
  84: (0,14),
  85: (0,56),
  86: (0,15),
  87: (0,28),
  88: (0,13),
  89: (0,52),
  90: (0,0),
  91: (1,61),
  92: (0,2),
  93: (-1,0),
  94: (0,3),
  95: (-1,0),
  96: (0,1),
  97: (-1,0),
  98: (0,23),
  99: (-1,0),
  100: (0,49),
  101: (-1,0),
  102: (0,10),
  103: (-1,0),
  104: (0,48),
  105: (-1,0),
  106: (0,42),
  107: (-1,0),
  108: (0,32),
  109: (-1,0),
  110: (0,34),
  111: (-1,0),
  112: (0,35),
  113: (-1,0),
  114: (0,33),
  115: (-1,0),
  116: (0,17),
  117: (-1,0),
  118: (0,8),
  119: (-1,0),
  120: (0,40),
  121: (-1,0),
  122: (0,41),
  123: (-1,0),
  124: (0,6),
  125: (-1,0),
  126: (0,4),
  127: (-1,0),
  128: (0,38),
  129: (-1,0),
  130: (0,7),
  131: (-1,0),
  132: (0,43),
  133: (0,61),
  134: (0,11),
  135: (0,12),
  136: (0,37),
  137: (0,44),
  138: (0,9),
  139: (0,46),
  140: (0,36),
  141: (0,47),
  142: (0,5),
  143: (0,45),
  144: (0,39),
  # 145: (-1,0),
  # 146: (-1,0)
}

detector_names = {
  -1: 'NC',
  0: 'Scint',
  1: 'Straw',
  2: 'MM0',
  3: 'MM1',
  4: 'MM2',
  5: 'MM3',
  6: 'AddStraw',
  7: 'lemo'
}

# pins_to_mapping
def print_connector(connector_num = 0):
  tiger_connector_map = {k: (f'{i+connector_num*2}_{j}' if i >= 0 else '') for k, (i,j) in tigChMap.items()}
  print_connector_map(connector_num, tiger_connector_map)

def print_connector_map(connector_num = 0, tiger_connector_map = {}):
  max_length_keys = max(map(lambda x: len(f'P{x}'), tiger_connector_map.keys()))
  max_length_vals = max(map(len, tiger_connector_map.values()))
  print()
  print(' '*(max_length_vals+3), f'connector{connector_num}')
  print(' '*(max_length_vals+2), '-'*(max_length_keys*2+4))
  for i in range(72):
    k1 = i*2+1
    k2 = i*2+2
    v1 = tiger_connector_map[k1] if k1 in tiger_connector_map.keys() else ''
    v2 = tiger_connector_map[k2] if k2 in tiger_connector_map.keys() else ''
    print(' '*(max_length_vals-len(v1)), f'{v1} | P{k1}',
          ' ' * (max_length_keys * 2 - len(str(k1)) - len(str(k2)) - 2),
          f'P{k2} | {v2}')
  print(' '*(max_length_vals+2), '-'*(max_length_keys*2+4))
  print()

def parse_connector(infile):
  with open(infile, 'r') as f:
    ls = [list(map(lambda x: x.strip(), l.split('#', 1)))
          for l in f.readlines() if l.strip()]
    parsed = [{}]
    for l in ls:
      if (not l[0]) and l[1]:
        # type: comment
        if 'comments' not in parsed[-1].keys():
          parsed[-1]['comments'] = []
        parsed[-1]['comments'].append(l[1])
    print(ls)

def convert_connectors_to_tiger_channels(infile):
  current_connector_map_p = [{} for p in range(8)]
  tiger_map = {}
  with open(infile, 'r') as f:
    for line in f.readlines():
      conf_text = line.split('#', 1)[0].strip()
      if not conf_text:
        continue
      # print(conf_text)
      conn, chann = tuple(map(lambda x: x.strip(), conf_text.split(':', 1)))
      connector, pin = tuple(map(lambda x: int(x.strip()), conn.split(',')))
      detector, strip = tuple(map(lambda x: x.strip(), chann.split(',')))
      detector_names_reverted = {v.lower(): k for k, v in detector_names.items()}
      detector = detector_names_reverted[detector.lower()] if detector.lower() in detector_names_reverted.keys() else int(detector)
      strip = int(strip)

      # print(connector, '-', pin, '--', detector, '-', strip)
      current_connector_map_p[connector][pin] = chann
      tiger, tpos = tigChMap[pin]
      tiger_map[(tiger + 2 * connector, tpos)] = (detector, strip)
  for c in range(len(current_connector_map_p)):
    tmap = current_connector_map_p[c]
    if len(tmap):
      print_connector_map(c, tmap)
  # sorting
  tiger_map = {k: v for k, v in sorted(tiger_map.items(), key=lambda item: item[0][1])} # by pin
  tiger_map = {k: v for k, v in sorted(tiger_map.items(), key=lambda item: item[0][0])} # by tiger
  with open(infile + '.out', 'w') as f:
    comment = '''#
# Map file for Tiger setup. 
# Structure: gemroc \t tiger # \t tiger channel \t detector \t channel
# Detectors:
# - -1: not connected
# - 0: scintillator:
#   - channel 0 - triple scintillator coinsidence [non-scaled]
#   - channel 1 - sci1
#   - channel 2 - sci2
#   - channel 3 - triple scintillator coinsidence prescaled
#   - channel 4 - master clock (10ms)
#   - channel 5 - master clock (50us)
# - 1: straw
# - 2: MM layer 0
# - 3: MM layer 1
# - 4: MM layer 2
# - 5: MM layer 3
# - 6: Additional straws:
#   - channel 0 - SHiP straw
#   - channel 1 - neutron straw
# - 7: lemo connector without additional description'''
    print(comment, file = f)
    for k, v in tiger_map.items():
      print(f'0\t{k[0]}\t{k[1]}\t{v[0]}\t{v[1]}', file = f)

# parse_connector('test.txt')

# print_connector(0)
infile = sys.argv[1] if len(sys.argv) > 1 else 'MM_TIGER_Cross_map.txt'
convert_connectors_to_tiger_channels(infile);
