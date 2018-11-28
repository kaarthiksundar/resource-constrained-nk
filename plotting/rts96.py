from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np
import re


num_buses = 24

problem_type = { 1: 'plain', 2 : 'planar', 3 : 'topological'}

k = [2, 3, 4, 5, 6]

D = list(range(10, 31, 2))

table_files = []

for (key, value) in problem_type.items():
    print(value)
    for i in k:
        file = '../output/24-{}-{}-heuristic.txt'.format(i, value)
        if key == 2: 
            file = '../output/24-{}-{}-5-heuristic.txt'.format(i, value)
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

for (key, value) in problem_type.items():
    for d in D:
        for j in k:
            file = '{}-{}-{}-{}-heuristic.txt'.format(num_buses, j, value, d)
            table_files.append(file)
