package dilivia.s2

import com.google.common.geometry.S2
import mu.KotlinLogging


// Defines functions related to determining whether two geodesic edges cross
// and for computing intersection points.
//
// The predicates CrossingSign(), VertexCrossing(), and EdgeOrVertexCrossing()
// are robust, meaning that they produce correct, consistent results even in
// pathological cases.  See s2predicates.h for additional robust predicates.
//
// See also S2EdgeCrosser (which efficiently tests an edge against a sequence
// of other edges) and S2CrossingEdgeQuery (which uses an index to speed up
// the process).
object S2EdgeCrossings {

    private val logger = KotlinLogging.logger {  }

    // kIntersectionError is an upper bound on the distance from the intersection
    // point returned by GetIntersection() to the true intersection point.
    val kIntersectionError: S1Angle = S1Angle.radians(8.0 * S2Predicates.kDoubleRoundingEpsilon)

    // This value can be used as the S2Builder snap_radius() to ensure that edges
    // that have been displaced by up to kIntersectionError are merged back
    // together again.  For example this can happen when geometry is intersected
    // with a set of tiles and then unioned.  It is equal to twice the
    // intersection error because input edges might have been displaced in
    // opposite directions.
    val kIntersectionMergeRadius: S1Angle = 2.0 * kIntersectionError

    // This function determines whether the edge AB intersects the edge CD.
    // Returns +1 if AB crosses CD at a point that is interior to both edges.
    // Returns  0 if any two vertices from different edges are the same.
    // Returns -1 otherwise.
    //
    // Note that if an edge is degenerate (A == B or C == D), the return value
    // is 0 if two vertices from different edges are the same and -1 otherwise.
    //
    // Properties of CrossingSign:
    //
    //  (1) CrossingSign(b,a,c,d) == CrossingSign(a,b,c,d)
    //  (2) CrossingSign(c,d,a,b) == CrossingSign(a,b,c,d)
    //  (3) CrossingSign(a,b,c,d) == 0 if a==c, a==d, b==c, b==d
    //  (3) CrossingSign(a,b,c,d) <= 0 if a==b or c==d (see above)
    //
    // This function implements an exact, consistent perturbation model such
    // that no three points are ever considered to be collinear.  This means
    // that even if you have 4 points A, B, C, D that lie exactly in a line
    // (say, around the equator), C and D will be treated as being slightly to
    // one side or the other of AB.  This is done in a way such that the
    // results are always consistent (see s2pred::Sign).
    //
    // Note that if you want to check an edge against a collection of other edges,
    // it is much more efficient to use an S2EdgeCrosser (see s2edge_crosser.h).
    fun crossingSign(a: S2Point, b: S2Point, c: S2Point, d: S2Point): Int {
        val crosser = S2EdgeCrosser(a, b, c)
        return crosser.crossingSign(d)
    }

    // Given two edges AB and CD where at least two vertices are identical
    // (i.e. CrossingSign(a,b,c,d) == 0), this function defines whether the
    // two edges "cross" in a such a way that point-in-polygon containment tests
    // can be implemented by counting the number of edge crossings.  The basic
    // rule is that a "crossing" occurs if AB is encountered after CD during a
    // CCW sweep around the shared vertex starting from a fixed reference point.
    //
    // Note that according to this rule, if AB crosses CD then in general CD
    // does not cross AB.  However, this leads to the correct result when
    // counting polygon edge crossings.  For example, suppose that A,B,C are
    // three consecutive vertices of a CCW polygon.  If we now consider the edge
    // crossings of a segment BP as P sweeps around B, the crossing number
    // changes parity exactly when BP crosses BA or BC.
    //
    // Useful properties of VertexCrossing (VC):
    //
    //  (1) VC(a,a,c,d) == VC(a,b,c,c) == false
    //  (2) VC(a,b,a,b) == VC(a,b,b,a) == true
    //  (3) VC(a,b,c,d) == VC(a,b,d,c) == VC(b,a,c,d) == VC(b,a,d,c)
    //  (3) If exactly one of a,b equals one of c,d, then exactly one of
    //      VC(a,b,c,d) and VC(c,d,a,b) is true
    //
    // It is an error to call this method with 4 distinct vertices.
    fun vertexCrossing(a: S2Point, b: S2Point, c: S2Point, d: S2Point): Boolean {
        // If A == B or C == D there is no intersection.  We need to check this
        // case first in case 3 or more input points are identical.
        if (a == b || c == d) return false

        // If any other pair of vertices is equal, there is a crossing if and only
        // if OrderedCCW() indicates that the edge AB is further CCW around the
        // shared vertex O (either A or B) than the edge CD, starting from an
        // arbitrary fixed reference point.
        //
        // Optimization: if AB=CD or AB=DC, we can avoid most of the calculations.
        if (a == c) return (b == d) || S2Predicates.orderedCCW(a.ortho(), d, b, a)
        if (b == d) return S2Predicates.orderedCCW(b.ortho(), c, a, b)

        if (a == d) return (b == c) || S2Predicates.orderedCCW(a.ortho(), c, b, a);
        if (b == c) return S2Predicates.orderedCCW(b.ortho(), d, a, b);

        logger.error("VertexCrossing called with 4 distinct vertices")
        return false;
    }

    // A convenience function that calls CrossingSign() to handle cases
    // where all four vertices are distinct, and VertexCrossing() to handle
    // cases where two or more vertices are the same.  This defines a crossing
    // function such that point-in-polygon containment tests can be implemented
    // by simply counting edge crossings.
    fun edgeOrVertexCrossing(a: S2Point, b: S2Point, c: S2Point, d: S2Point): Boolean {
        TODO()
    }

    // Given two edges AB and CD such that CrossingSign(A, B, C, D) > 0, returns
    // their intersection point.  Useful properties of GetIntersection (GI):
    //
    //  (1) GI(b,a,c,d) == GI(a,b,d,c) == GI(a,b,c,d)
    //  (2) GI(c,d,a,b) == GI(a,b,c,d)
    //
    // The returned intersection point X is guaranteed to be very close to the
    // true intersection point of AB and CD, even if the edges intersect at a
    // very small angle.  See "kIntersectionError" below for details.
    fun getIntersection(a: S2Point, b: S2Point, c: S2Point, d: S2Point): S2Point {
        TODO()
    }

}