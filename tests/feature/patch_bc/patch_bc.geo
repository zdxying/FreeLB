[geometry]
grid = 12, 12, 8
domain = 0.0, 0.0, 0.0, 12.0, 12.0, 8.0
cell_len = 1.0
band_width = 1.5

[source.domain]
shape = box
entity = 1
material = 1
role = domain
box = 0.0, 0.0, 0.0, 12.0, 12.0, 8.0

[source.obstacle]
shape = box
entity = 2
material = 2
role = solid
box = 4.0, 4.0, 2.0, 8.0, 8.0, 6.0

[source.plate]
shape = box
entity = 3
material = 3
role = solid
box = 1.0, 1.0, 1.0, 3.0, 3.0, 7.0

[patch.obstacle]
entity = 2
kind = wall
box = 2.5, 2.5, 0.5, 9.5, 9.5, 7.5

[patch.plate]
entity = 3
kind = wall
box = -0.5, -0.5, -0.5, 4.5, 4.5, 8.5
