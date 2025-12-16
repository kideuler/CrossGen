A purely 2D implementation of 

CrossGen: Learning and Generating Cross Fields for Quad Meshing

TODO:

STEP 1:
- [ ] Create some 2D geometries of interest in gmsh.
- [ ] scale and Mesh them and save meshes with 1/300 edge length in bounding box.
- [ ] Convert meshes to .obj file.
- [ ] Run .obj though Quadriflow and save .obj quad meshes.
- [ ] inspect meshes for validity in paraview.
- [ ] create python scripts to grab mesh data, rasterize, get crossfield data and save it. 

STEP 2:
- [ ] Get training script with pytorch up and running.
- [ ] Train and validate.
- [ ] Hook in quadwild to test on new meshes.
- [ ] Test on some interesting examples that network has not seen.

STEP 3: (not determined)