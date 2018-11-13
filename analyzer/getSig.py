import math
merge = []
merge2 = []
sig = []
tmp = []

with open("Nevt.txt") as fp:
  for line in fp.readlines():
    line.replace("\n", "")

    if '[' in line:
      sig.append(line.split("\', \'")[2])
      tmp = line.split("\', \'")[3]
      sig.append(tmp.replace("'] \n",""))
    if 'muta' in line:
      if 'hist' in line: continue
      tmp = line.split(" ")
      sig.append(float(tmp[5]))
      sig.append(float(tmp[6]))
      sig.append(float(tmp[7]))
    if 'ntotal' in line:
      merge.append(sig)
      sig = []

for items in merge:
  if items[0] in ['Ch0', 'Ch1', 'Ch3', 'Ch4', 'Ch5']: continue
  if math.isnan(items[2]) or math.isnan(items[3]) or math.isnan(items[4]): continue
  merge2.append(items)

print merge2

merge2.sort(key = lambda x: x[2], reverse = True)
print "high 10 muta: "
for i in range(0,10):
  print str(merge2[i])

merge2.sort(key = lambda x: x[3], reverse = True)
print "high 10 tata: "
for i in range(0,10):
  print str(merge2[i])

merge2.sort(key = lambda x: x[4], reverse = True)
print "high 10 nunu: "
for i in range(0,10):
  print str(merge2[i])
