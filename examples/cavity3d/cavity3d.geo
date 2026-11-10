# cavity3d.geo
# Lid-driven cavity with an internal solid obstacle.  The domain and the
# obstacle are the geometry file; the obstacle's Wall patch drives its
# bounce-back cells, while the domain-edge walls and the moving lid stay
# AABB-driven (they have no surface band to attach a patch to).
#
# The grid and cell_len must match [Mesh] in cavity3d.ini.

[geometry]
grid = 100, 100, 100
domain = 0.0, 0.0, 0.0, 100.0, 100.0, 100.0
cell_len = 1.0
band_width = 1.5

[source.domain]
shape = box
entity = 1
material = 1
role = domain
box = 0.0, 0.0, 0.0, 100.0, 100.0, 100.0

[source.obstacle]
shape = box
entity = 2
material = 2
role = solid
box = 40.0, 40.0, 0.0, 60.0, 60.0, 70.0

[patch.obstacle]
entity = 2
kind = wall
box = 38.5, 38.5, 0.0, 61.5, 61.5, 71.5
