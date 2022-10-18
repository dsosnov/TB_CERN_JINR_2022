# Since MM started from 1,
# First connector: MM strips 1 -- 122
# Second connector: MM strips 123 -- 238
# third connector: MM strips 239 -- 360
# Also, on channel MM_0 in MM_Mu2E_Cross.pdf set the ground, channel removed

mm_connector_strips={
  # type: N strips in previous connectors, N strips before GND, N GND strips, N strips after GND 
  'l': (0, 32, 6, 26),
  'c': (128-6, 27, 12, 25),
  'r': (128*2-18, 26, 6, 32)
  }
connector_desc={
  'l': 'left',
  'c': 'central',
  'r': 'right'
  }

feb=1
mm=0
for connector_type in ["l", "c", "r"]:
  desc = connector_desc[connector_type]
  with open(f'mm_tiger_board_bautin_template_{connector_type}.txt', 'w') as f:
    print(f'# tigers connected to {desc} connector of MM{mm} with Bautin\'s board', file=f)
    print(f'# connector (FEB), pin (on FEB3_pinout scheme) : detector, channel', file=f)
    print(f'{feb}, 10 :  NC, 0', file=f)
    params = mm_connector_strips[connector_type]
    for i in(range(12,140,2)):
      # channel = 1 + (128-6) + 27 + int((i-12)/2)
      channel = 1 + params[0] + params[1] + int((i-12)/2)
      print(f'{feb}, {i} : MM{mm}, {channel}', file=f)
    lemo = [5,11,17,23,29,35,41,47]
    for ch in lemo:
      print(f'{feb}, {ch} : lemo, {lemo.index(ch)} # lemo {lemo.index(ch)}', file=f)
