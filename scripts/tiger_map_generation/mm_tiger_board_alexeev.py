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
  with open(f'mm_tiger_board_alexeev_template_{connector_type}.txt', 'w') as f:
    print(f'# tigers connected to {desc} connector of MM{mm} with gem-tiger board', file=f)
    print(f'# connector (FEB), pin (on FEB3_pinout scheme) : detector, channel', file=f)
    params = mm_connector_strips[connector_type]
    with open('mm_tiger_board_alexeev_pins.csv', 'r') as f_pins:
      for line in f_pins.readlines():
        line_splitted = list(map(lambda x: x.strip(), line.split(';')))
        if line_splitted[0].startswith('#'):
          continue;
        if line_splitted[1].startswith('GND'):
          continue;
        pin_num = int(line_splitted[0])
        pin_num_feb3 = int(line_splitted[3])
        channel = -1
        if not (pin_num % 2):
          channel = params[0] + params[1] + (pin_num - 2) / 2
        elif pin_num < params[1]*2:
          channel = params[0] + params[1] - (pin_num-1)/2
        elif pin_num < (params[1]+params[2])*2:
          channel = -1
        else:
          channel = params[0] + 64 + params[1] + 1 + (127-pin_num)/2
        channel = int(channel)
        if channel >= 0:
          print(f'{feb}, {pin_num_feb3} :  MM{mm}, {channel}', file=f)
        else:
          print(f'{feb}, {pin_num_feb3} :  NC, 0', file=f)
