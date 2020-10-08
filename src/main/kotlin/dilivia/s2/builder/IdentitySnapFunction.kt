package dilivia.s2.builder

import dilivia.s2.Assertions
import dilivia.s2.S1Angle
import dilivia.s2.S2Point

// A SnapFunction that snaps every vertex to itself.  It should be used when
// vertices do not need to be snapped to a discrete set of locations (such as
// E7 lat/lngs), or when maximum accuracy is desired.
//
// If the given "snap_radius" is zero, then all input vertices are preserved
// exactly.  Otherwise, S2Builder merges nearby vertices to ensure that no
// vertex pair is closer than "snap_radius".  Furthermore, vertices are
// separated from non-incident edges by at least "min_edge_vertex_separation",
// equal to (0.5 * snap_radius).  For example, if the snap_radius is 1km, then
// vertices will be separated from non-incident edges by at least 500m.
// REQUIRES: snap_radius <= SnapFunction::kMaxSnapRadius()
class IdentitySnapFunction(snapRadius: S1Angle = S1Angle.zero) : SnapFunction() {

    override var snapRadius: S1Angle = snapRadius
        set(value) {
            Assertions.assertLE(value, kMaxSnapRadius)
            field = value
        }

    // For the identity snap function, all vertex pairs are separated by at
    // least snap_radius().
    // Since SnapFunction does not move the input point, output vertices are
    // separated by the full snap_radius().
    override fun minVertexSeparation(): S1Angle = snapRadius

    // For the identity snap function, edges are separated from all non-incident
    // vertices by at least 0.5 * snap_radius().
    // In the worst case configuration, the edge separation is half of the
    // vertex separation.
    override fun minEdgeVertexSeparation(): S1Angle  = snapRadius * 0.5

    override fun snapPoint(point: S2Point): S2Point = point

    override fun clone(): SnapFunction = IdentitySnapFunction(snapRadius)

}
