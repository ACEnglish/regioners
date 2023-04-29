import math
n = 10
t = 4

csize = math.ceil(n / t)  #floor
print(csize)
for i in range(t):
    print('starting thread', i)
    start = i*csize
    for j in range(start, start + csize):
        print(j)
