from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np
import re


num_buses = 240

problem_type = { 1: 'plain', 2 : 'planar', 3 : 'topological'}

k = [2, 3, 4, 5, 6]

D = list(range(100, 1001, 100))

load_shed_values = {}
for i in k:
    load_shed_values[i] = {'x': [], 'y': []}

for i in k:
    for d in D:
        file = '../output/240-{}-planar-{}-heuristic.txt'.format(i, d)
        lines = open(file, 'r')
        for line in lines:
            if re.match("iterations(.*)", line):
                load_shed_values[i]['x'].append(d) 
                load_shed_values[i]['y'].append(float(re.split('\s+', line.rstrip())[-1]))



'''
for (key, value) in problem_type.items():
    print(value)
    for i in k:
        file = '../output/240-{}-{}-heuristic.txt'.format(i, value)
        if key == 2: 
            file = '../output/240-{}-{}-500-heuristic.txt'.format(i, value)
        lines = open(file, 'r')
        load_shed = 0.0 
        time = 0.0 
        iterations = 0
        opt_gap = 0.0
        for line in lines:
            if re.match("opt(.*)", line):
                opt_gap = float(re.split('\s+', line.rstrip())[-1])
            if re.match("load(.*)", line):
                load_shed = float(re.split('\s+', line.rstrip())[-1])
            if re.match("(.*)time(.*)", line):
                time = float(re.split('\s+', line.rstrip())[-1])
            if re.match("(.*)iterations(.*)", line):
                iterations = int(re.split('\s+', line.rstrip())[-1])
        print('{} & {} & {} & {} & {} \\\\'.format(i, load_shed, time, iterations, opt_gap))

print()
for i in k:
    for d in D:
        file = '../output/240-{}-planar-{}-heuristic.txt'.format(i, d)
        lines = open(file, 'r')
        load_shed = 0.0 
        time = 0.0 
        iterations = 0
        for line in lines:
            if re.match("load(.*)", line):
                load_shed = float(re.split('\s+', line.rstrip())[-1])
            if re.match("(.*)time(.*)", line):
                time = float(re.split('\s+', line.rstrip())[-1])
            if re.match("(.*)iterations(.*)", line):
                iterations = int(re.split('\s+', line.rstrip())[-1])
        print('{} & {} & {} & {} \\\\'.format(i, load_shed, time, iterations))
'''

plt.rcParams["font.family"] = "Times New Roman"
plt.rc('text', usetex=True)

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.view_init(elev=31, azim=-143)
ax.grid(False)
for c, z in zip(['r', 'g', 'b', 'y', 'c'], [2, 3, 4, 5, 6]):
    xs = [1,2,3,4,5,6,7,8,9,10]
    ys = load_shed_values[z]['y']

    # You can provide either a single color or an array. To demonstrate this,
    # the first bar of each set will be colored cyan.
    cs = [c] * len(xs)
    ax.bar(xs, ys, zs=z, zdir='y', color=cs, alpha=0.5)

xlabels = [r'$D=200$ km', r'$D=600$ km', r'$D=1000$ km']
ax.set_xticklabels(xlabels)
ax.set_xticks([2, 6, 10])


ylabels = [r'$k=2$', r'$k=3$', r'$k=4$', r'$k=5$', r'$k=6$']
ax.set_yticklabels(ylabels)
ax.set_yticks([2, 3, 4, 5, 6])

zlabels = [0, 10, 20, 30, 40]
ax.set_zticklabels(zlabels)
ax.set_zticks([0, 10, 20, 30, 40])


ax.set_zlabel('iterations')
plt.savefig('iterations.pdf', format='pdf')