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
for inverted in [False, True]:
  inv_name = '_inv' if inverted else ''
  inv_text = ', inverted' if inverted else ''
  for connector_type in ["l", "c", "r"]:
    desc = connector_desc[connector_type]
    with open(f'mm_tiger_board_alexeev_template_{connector_type}{inv_name}.txt', 'w') as f:
      print(f'# tigers connected to {desc} connector of MM{mm} with gem-tiger board{inv_text}', file=f)
      print(f'# connector (FEB), pin (on FEB3_pinout scheme) : detector, channel', file=f)
      params = mm_connector_strips[connector_type]
      for j1_pin in range(1, 131):
        mm_strip = -1
        if j1_pin == 2 or j1_pin == 139:
          continue
        elif j1_pin % 2 and j1_pin >= 121: # 121, 123, 125, 127
          continue
        elif not j1_pin % 2:
          mm_strip = params[0] + 1 + params[1] + (j1_pin - 4) // 2
        elif j1_pin < params[1] * 2:
          mm_strip = params[0] + params[1] - (j1_pin - 1) // 2
        elif j1_pin < (params[1] + params[2]) * 2:
          mm_strip = -1
        else:
          mm_strip = params[0] + params[1] + 64 + 1 + (127 - j1_pin) // 2

        ch_number = j1_pin + 1 if j1_pin%2 else j1_pin - 3
        
        j2_pin = -1
        if ch_number%2:
          j2_pin = 136 - (ch_number-1) // 2
        elif ch_number <= 8:
          j2_pin = 141 + (ch_number - 2) // 2
        elif ch_number <= 16:
          j2_pin = 137 + (16 - ch_number) // 2
        elif ch_number <= 28:
          j2_pin = 67 + (28 - ch_number) // 2
        elif ch_number <= 120:
          j2_pin = 1 + (120 - ch_number) // 2
        else:
          j2_pin = -1

        if j2_pin < 0:
          print(f'Feb pin for panasonic pin is below zero! [{j1_pin} -> {mm_strip} -> {ch_number}]')
          continue
        
        feb_pin = 1 + (j2_pin - 1) * 2 if j2_pin <= 72 else 2 + (j2_pin - 73) * 2
        if mm_strip < 0:
          print(f'{feb}, {feb_pin} : NC, 0', file=f)
        else:
          if(inverted):
            mm_strip = 361 - mm_strip
          print(f'{feb}, {feb_pin} : MM{mm}, {mm_strip}', file=f)
