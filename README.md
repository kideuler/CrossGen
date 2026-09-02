A purely 2D implementation of Cross-fields along with block decomposition

# TODO:
- [X] Implement polyvectors to obtain crossfields.
- [X] Robustly detect and classify singularities.
- [X] Cut Mesh for MIQ.
- [X] Link Comisol for Greedy mixed integer solver.
- [X] construct UV integer grid map.
- [X] Some basic Visualization.
- [X] Reformat Code to better separate different pieces. (Mesh, Polyvectors, IGM)
- [X] In polyvector store the parallel transport for each triangle (0,1,2,3) rotation matrix
- [X] Add options to solve MBO method using formula $(M+\tau L)u^{k+1}=Mu^k$ where M is diagonal triangle areas. L is standard laplacian. Therefore $u^{k+1} = (M+\tau L)^{-1} M u^k$
- [X] Make step operations on separatrix instead of on tracePoints
- [X] Test singular triangle, for 3,5 singularities, one with directions going towards singularity, 4 tests total.
- [X] Implement DG-SIPG to in P0 to get crosses on faces and singularities on nodes.
- [X] Implement basic polysquare stuff from 2023 paper
- [X] Go back and try to fix tracing and simplification
- [X] See if we can do a chord collapse operation for complex polysquares
- [x] Try medial axis again.
- [X] Change mesh cutting a little so singularities do not lie on the same cut. Should conditions better.
- [X] Tracing.
- [X] Multi-materials in MERIDIAN through alignment to feature lines in energy.
- [X] SIPG integration through TORSION
- [X] Disk material pinning rotational symmetry
- [X] Templated circular inclusions: excise them (Stage 0c) and fill each with an O-grid (Stage 11).
- [X] Get Laghos to output to paraview.
- [X] Get bubble working correctly (`--disk-templates --near-miss 0.04`, both pipelines)
- [ ] Create a QuadMesh class and do TMOP.
- [ ] Get tests ready.
- [ ] Try to parallelize (for a little)
