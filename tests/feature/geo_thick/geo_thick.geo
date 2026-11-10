# geo_thick.geo
# Thick-wall / multi-entity fixture for B3/B4.  Two 4-cell-thick solid plates
# with distinct entities and materials sit in a 32^3 fluid domain; the gap
# between them is fluid, so a wall patch or a link boundary can target either
# side without a one-cell staircase.

[geometry]
grid = 32, 32, 32
domain = 0, 0, 0, 32, 32, 32
cell_len = 1
band_width = 1.5

[source.domain]
shape = box
box = 0, 0, 0, 32, 32, 32
entity = 1
material = 1
role = domain

[source.wall_a]
shape = box
box = 8, 8, 8, 12, 24, 24
entity = 2
material = 2
role = solid

[source.wall_b]
shape = box
box = 20, 8, 8, 24, 24, 24
entity = 3
material = 3
role = solid
