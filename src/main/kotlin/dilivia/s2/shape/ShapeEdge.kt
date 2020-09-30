package dilivia.s2.shape

// A class representing a ShapeEdgeId together with the two endpoints of that
// edge.  It should be passed by reference.
data class ShapeEdge(val id: ShapeEdgeId, val edge: S2Shape.Edge) {
    constructor(shapeId: Int, edgeId: Int, edge: S2Shape.Edge) : this(ShapeEdgeId(shapeId, edgeId), edge)
    constructor(shape: S2Shape, edge_id: Int) : this(shape.id, edge_id, shape.edge(edge_id))

    val v0 = edge.v0
    val v1 = edge.v1

}
