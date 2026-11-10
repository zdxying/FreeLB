# geo_bulk.geo
# Bulk fixture for B3/B5: a 32^3 domain with a 24^3 thick solid block.
# The fluid fraction is 18944 / 32768 = 0.578, so block fill is not dominated
# by the thin-duct halo overhead.

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

[source.block]
shape = box
box = 4, 4, 4, 28, 28, 28
entity = 2
material = 2
role = solid
