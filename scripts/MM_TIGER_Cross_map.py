# Since MM started from 1,
# First APV: MM strips 1 -- 122
# Second:
# third:
# Also, on channel MM_0 in MM_Mu2E_Cross.pdf set the ground, channel removed
with open('MM_Mu2E_Cross_map.txt', 'w') as f:
  print(f'# connector, pin : detector, channel', file=f)
  for i in(range(12,140,2)):
    channel = 1 + (128-6) + 27 + int((i-12)/2)
    print(f'1, {i} : MM0, {channel}', file=f)
  lemo = [11,17,23,29,35,41,47]
  for ch in lemo:
    print(f'1, {ch} : lemo, {lemo.index(ch)}', file=f)
