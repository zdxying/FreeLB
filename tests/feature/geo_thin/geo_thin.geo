# geo_thin.geo
# Thin-wall fixture for B3/B4: two 1-cell-thick plates in a 32^3 fluid domain.
# The plates are deliberately one cell thick so a later link-BC branch can
# exercise a wall thinner than the D3Q19 stencil reach.

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

[source.plate_a]
shape = box
box = 8, 4, 4, 9, 28, 28
entity = 2
material = 2
role = solid

[source.plate_b]
shape = box
box = 20, 4, 4, 21, 28, 28
entity = 3
material = 3
role = solid
