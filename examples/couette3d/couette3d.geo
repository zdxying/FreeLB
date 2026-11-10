[geometry]
grid = 16, 16, 10
domain = 0.0, 0.0, 0.0, 16.0, 16.0, 10.0
cell_len = 1.0
band_width = 1.5

[source.domain]
shape = box
entity = 1
material = 1
role = domain
box = 0.0, 0.0, 0.0, 16.0, 16.0, 10.0

[source.bottom]
shape = box
entity = 2
material = 2
role = solid
box = 0.0, 0.0, -1.0, 16.0, 16.0, 1.0

[source.top]
shape = box
entity = 3
material = 3
role = solid
box = 0.0, 0.0, 9.0, 16.0, 16.0, 11.0

[patch.bottom]
entity = 2
kind = wall
box = -1.0, -1.0, -1.0, 17.0, 17.0, 2.5

[patch.top]
entity = 3
kind = wall
box = -1.0, -1.0, 7.5, 17.0, 17.0, 12.0
