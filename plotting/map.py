import folium
import csv
import json
import branca
import re


attr = ('&copy; <a href="http://www.openstreetmap.org/copyright">OpenStreetMap</a> '
        'contributors, &copy; <a href="http://cartodb.com/attributions">CartoDB</a>')
tiles = 'http://{s}.basemaps.cartocdn.com/light_nolabels/{z}/{x}/{y}.png'
white_tile = branca.utilities.image_to_url([[1, 1], [1, 1]])
m = folium.Map(location=[42.397411, -114.068631],  zoom_start=5, tiles=tiles, attr=attr, zoom_control=False, width='38%', height='100%', prefer_canvas=True)

geo_json_data = json.load(open('canada-states.json'))

folium.GeoJson(geo_json_data, style_function=lambda feature: {
        'fillColor': 'yellow',
        'color': 'black',
        'opacity': 0.6,
        'weight': 1,
    }).add_to(m)

geo_json_data_us = json.load(open('us-states.json'))

folium.GeoJson(geo_json_data_us, style_function=lambda feature: {
        'fillColor': 'yellow',
        'color': 'black',
        'opacity': 0.6,
        'weight': 1,
    }).add_to(m)

f = open('../data/240_geo.csv', 'r')
reader = csv.reader(f)
rows = []  
for row in reader:
    rows.append(row)
position = {}
for i in range(1, len(rows)):
    row = rows[i]
    bus = int(row[0])
    lat = float(row[1])
    lon = float(row[2])
    position[bus] = (lat, lon)
    
f = open('../data/240_branches.csv', 'r')
reader = csv.reader(f)
rows = []
for row in reader:
    rows.append(row)
branches = {}
for i in range(0, len(rows)):
    row = rows[i]
    f_bus = int(row[0])
    t_bus = int(row[1])
    branches[i] = [list(position[f_bus]), list(position[t_bus])]

for key, value in branches.items():
    folium.PolyLine(value, color='green', weight=1.5).add_to(m)

for key, value in position.items():
    folium.CircleMarker(location=list(value), radius=2, weight=1.5, color='brown', fill_color='brown', fill_opacity=1).add_to(m)


m.save('wecc240.html')

# traditional 
k = 6 
file = '../output/240-{}-planar-500-heuristic.txt'.format(k)
f = open(file, 'r')
branches = []
bus_shed = {}
lines = list(f)
i = 0
while (i < len(lines)):
    line = lines[i]
    if re.match("interdicted(.*)", line):
        i += 1 
        while (bool(re.match("bus(.*)", lines[i])) == False):
            new_line = lines[i]
            parsed = re.split('\s+', new_line.rstrip())
            f_bus = int(parsed[0])
            t_bus = int(parsed[1])
            branches.append([f_bus, t_bus])
            i += 1
        i += 1
        while (bool(re.match("time(.*)", lines[i])) == False):
            new_line = lines[i]
            parsed = re.split('\s+', new_line.rstrip())
            bus_id = int(parsed[0])
            load_shed = float(parsed[1])
            bus_shed[bus_id] = load_shed
            i += 1
    i += 1

for branch in branches:
    f_bus = branch[0]
    t_bus = branch[1]
    value = [list(position[f_bus]), list(position[t_bus])]
    # folium.PolyLine(value, color='red', weight=2).add_to(m)

total_shed = sum(bus_shed.values())
for (bus_id, value) in bus_shed.items():
    folium.CircleMarker(
    location=list(position[bus_id]),
    radius=value/total_shed*100,
    color='black',
    weight=1,
    fill_opacity=0.6,
    opacity=1,
    fill_color='black').add_to(m)

m.save('wecc240-planar.html')
# planar 


# topological

'''
m = folium.Map(zoom_start=5, tiles=tiles, attr=attr, zoom_control=False, width='38%', height='95%', prefer_canvas=True)
f = open('../data/24_geo.csv', 'r')
reader = csv.reader(f)
rows = []  
for row in reader:
    rows.append(row)
position = {}
for i in range(1, len(rows)):
    row = rows[i]
    bus = int(row[0])
    lat = float(row[1])
    lon = float(row[2])
    position[bus] = (lat, lon)
    
f = open('../data/24_branches.csv', 'r')
reader = csv.reader(f)
rows = []
for row in reader:
    rows.append(row)
branches = {}
for i in range(0, len(rows)):
    row = rows[i]
    f_bus = int(row[0])
    t_bus = int(row[1])
    branches[i] = [list(position[f_bus]), list(position[t_bus])]

for key, value in branches.items():
    folium.PolyLine(value, color='green', weight=1.5).add_to(m)

for key, value in position.items():
    folium.CircleMarker(location=list(value), radius=2, weight=1.5, color='black', fill_color='brown', fill_opacity=1).add_to(m)

m.save('rts96.html')
'''