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
for connector_type in ['l', 'c', 'r']: # 'l', 'c', 'r'
  desc = connector_desc[connector_type]
  with open(f'mm_tiger_board_bautin-2_template_{connector_type}.txt', 'w') as f:
    print(f'# tigers connected to {desc} connector of MM{mm} with B2-board', file=f)
    print(f'# connector (FEB), pin (on FEB3_pinout scheme) : detector, channel', file=f)
    params = mm_connector_strips[connector_type]
    for panasonic_pin in range(1, 131):
      mm_strip = -1
      mm_in_scheme = -1
      if panasonic_pin == 1 or panasonic_pin == 66:
        continue
      elif panasonic_pin >= 93 and panasonic_pin <= 98:
        continue
      elif panasonic_pin < 66:
        mm_in_scheme = 30 + panasonic_pin
        mm_strip = params[0] + params[1] - 1 + panasonic_pin
      elif panasonic_pin <= 66 + params[3]:
        mm_strip = params[0] + params[1] + 64 + 1 + (panasonic_pin - 67)
        mm_in_scheme = 30 + panasonic_pin - 1 if panasonic_pin < 99 else panasonic_pin - 99
      elif panasonic_pin <= 66 + params[3] + params[2]:
        mm_in_scheme = 30 + panasonic_pin - 1 if panasonic_pin < 99 else panasonic_pin - 99
        mm_strip = -1 # connected to GND
      else:
        mm_strip = params[0] + 1 + (panasonic_pin - 67 - params[2] - params[3])
        mm_in_scheme = 30 + panasonic_pin - 1 if panasonic_pin < 99 else panasonic_pin - 99
      feb_pin = -1
      if mm_in_scheme <= 5:
        feb_pin = 81 + mm_in_scheme * 2
      elif mm_in_scheme <= 11:
        feb_pin = 69 + (11 - mm_in_scheme) * 2
      elif mm_in_scheme <= 25:
        if mm_in_scheme % 2: # odd
          feb_pin = 2 + (mm_in_scheme - 13)
        else:
          feb_pin = 1 + (24 - mm_in_scheme)
      elif mm_in_scheme <= 31:
        feb_pin = 133 + (mm_in_scheme - 26) * 2
      elif mm_in_scheme <= 95:
        feb_pin = 142 - (mm_in_scheme - 32) * 2 # == 16 + (95 - mm_in_scheme) * 2
      else:
        feb_pin = 17 + (mm_in_scheme - 96) * 2
      if feb_pin < 0:
        print(f'Feb pin for panasonic pin is below zero! [{panasonic_pin} -> {mm_strip} -> {mm_in_scheme}]')
      elif mm_strip < 0:
        print(f'{feb}, {feb_pin} :  NC, 0', file=f)
      else:
        print(f'{feb}, {feb_pin} : MM{mm}, {mm_strip}', file=f)
      # print(f'{panasonic_pin} -> {mm_strip} -> {mm_in_scheme} -> {feb_pin}')
    lemo = [144, 15]
    for ch in lemo:
      print(f'{feb}, {ch} : lemo, {lemo.index(ch)} # lemo {lemo.index(ch)}', file=f)
