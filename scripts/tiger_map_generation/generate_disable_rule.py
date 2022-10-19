import os, sys

filename = sys.argv[1] if len(sys.argv) > 1 else ""
print(f'Input file: {filename}')

if not filename:
  exit(0)

active_channels = dict()

with open(filename, 'r') as f:
  lines = f.readlines()
  for l in lines:
    l = l.strip()
    if not l or l[0] == '#':
      continue
    (gr,t,c,d,_) = map(lambda x: int(x.strip()), l.split())
    if gr not in active_channels.keys():
      active_channels[gr] = dict()
    if t not in active_channels[gr].keys():
      active_channels[gr][t] = set()
    if d >= 0:
      active_channels[gr][t].add(c)
with open('ch_to_disable', 'w') as fout:
  print('#Comments line begin with "#"', file=fout)
  print('#GEMROC TIGER CH (also wrote first:last)', file=fout)
  for gr in active_channels.keys():
    for t in active_channels[gr].keys():
      for ch in range(64):
        if ch not in active_channels[gr][t]:
          print(f'{gr} {t} {ch}', file=fout)

