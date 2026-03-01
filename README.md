A purely 2D implementation of Cross-fields along with block decomposition

TODO:
- [X] Implement polyvectors to obtain crossfields.
- [X] Robustly detect and classify singularities.
- [X] Cut Mesh for MIQ.
- [X] Link Comisol for Greedy mixed integer solver.
- [X] construct UV integer grid map.
- [X] Some basic Visualization.
- [X] Reformat Code to better separate different pieces. (Mesh, Polyvectors, IGM)
- [X] In polyvector store the parallel transport for each triangle (0,1,2,3) rotation matrix
- [X] Add options to solve MBO method using formula $(M+\tau L)u^{k+1}=Mu^k$ where M is diagonal triangle areas. L is standard laplacian. Therefore $u^{k+1} = (M+\tau L)^{-1} M u^k$
- [ ] Implement separatrix tracing (no topology simplification, no edge-maps) 
- [ ] Create separatrix graph including t-junctions.
- [ ] Do topological simplification as listed in the viertel imr 2019 paper
- [ ] Show all steps in gui and get some example quad meshes.
- [ ] Spin block topology in 3D for 90 and 180 degrees.
- [ ] Smooth hex blocks using 3rd order TMOP on Q3 elements. 