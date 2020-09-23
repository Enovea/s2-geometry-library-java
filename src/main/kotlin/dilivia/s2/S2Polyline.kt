package dilivia.s2

import com.google.common.geometry.S2.M_PI
import dilivia.s2.Assertions.assert
import dilivia.s2.Assertions.assertGE
import dilivia.s2.Assertions.assertLE
import dilivia.s2.Assertions.assertLT
import mu.KotlinLogging
import kotlin.math.*

// An S2Polyline represents a sequence of zero or more vertices connected by
// straight edges (geodesics).  Edges of length 0 and 180 degrees are not
// allowed, i.e. adjacent vertices should not be identical or antipodal.
class S2Polyline(val vertices: List<S2Point>) : S2Region {

    // Creates an empty S2Polyline that should be initialized by calling Init()
    // or Decode().
    constructor() : this(emptyList())

    init {
        assert { isValid() }
    }

    // Return true if the given vertices form a valid polyline.
    fun isValid(): Boolean {
        val error = findValidationError()
        if (!error.isOk()) {
            logger.error { error.message }
            return false;
        }
        return true;
    }

    // Returns true if this is *not* a valid polyline and sets "error"
    // appropriately.  Otherwise returns false and leaves "error" unchanged.
    //
    // REQUIRES: error != nullptr
    fun findValidationError(): S2Error {
        // All vertices must be unit length.
        for (i in 0 until numVertices()) {
            if (!S2Point.isUnitLength(vertex(i))) {
                return S2Error(code = S2Error.NOT_UNIT_LENGTH, "Vertex %d is not unit length".format(i))
            }
        }
        // Adjacent vertices must not be identical or antipodal.
        for (i in 1 until numVertices()) {
            if (vertex(i - 1) == vertex(i)) {
                return S2Error(code = S2Error.DUPLICATE_VERTICES, "Vertices %d and %d are identical".format(i - 1, i))
            }
            if (vertex(i - 1) == -vertex(i)) {
                return S2Error(code = S2Error.ANTIPODAL_VERTICES, "Vertices %d and %d are antipodal".format(i - 1, i))
            }
        }
        return S2Error()
    }

    fun numVertices(): Int = vertices.size
    fun vertex(k: Int): S2Point {
        assertGE(k, 0)
        Assertions.assertLT(k, numVertices())
        return vertices[k]
    }

    // Return the length of the polyline.
    fun getLength(): S1Angle = getLength(vertices)

    // Return the true centroid of the polyline multiplied by the length of the
    // polyline (see s2centroids.h for details on centroids).  The result is not
    // unit length, so you may want to normalize it.
    //
    // Prescaling by the polyline length makes it easy to compute the centroid
    // of several polylines (by simply adding up their centroids).
    fun getCentroid(): S2Point = getCentroid(vertices)

    // Return the point whose distance from vertex 0 along the polyline is the
    // given fraction of the polyline's total length.  Fractions less than zero
    // or greater than one are clamped.  The return value is unit length.  This
    // cost of this function is currently linear in the number of vertices.
    // The polyline must not be empty.
    fun interpolate(fraction: Double): S2Point = getSuffix(fraction).first

    // Like Interpolate(), but also return the index of the next polyline
    // vertex after the interpolated point P.  This allows the caller to easily
    // construct a given suffix of the polyline by concatenating P with the
    // polyline vertices starting at "next_vertex".  Note that P is guaranteed
    // to be different than vertex(*next_vertex), so this will never result in
    // a duplicate vertex.
    //
    // The polyline must not be empty.  Note that if "fraction" >= 1.0, then
    // "next_vertex" will be set to num_vertices() (indicating that no vertices
    // from the polyline need to be appended).  The value of "next_vertex" is
    // always between 1 and num_vertices().
    //
    // This method can also be used to construct a prefix of the polyline, by
    // taking the polyline vertices up to "next_vertex - 1" and appending the
    // returned point P if it is different from the last vertex (since in this
    // case there is no guarantee of distinctness).
    fun getSuffix(fraction: Double): Pair<S2Point, Int /* next vertex */> {
        Assertions.assertGT(numVertices(), 0)
        var next_vertex = 0
        // We intentionally let the (fraction >= 1) case fall through, since
        // we need to handle it in the loop below in any case because of
        // possible roundoff errors.
        if (fraction <= 0) {
            return vertex(0) to 1
        }
        var length_sum = S1Angle()
        for (i in 1 until numVertices()) {
            length_sum += S1Angle(vertex(i - 1), vertex(i))
        }
        var target = fraction * length_sum
        for (i in 1 until numVertices()) {
            val length = S1Angle(vertex(i - 1), vertex(i))
            if (target < length) {
                // This interpolates with respect to arc length rather than
                // straight-line distance, and produces a unit-length result.
                val result = S2EdgeDistances.interpolateAtDistance(target, vertex(i - 1), vertex(i));
                // It is possible that (result == vertex(i)) due to rounding errors.
                next_vertex = if (result == vertex(i)) (i + 1) else i;
                return result to next_vertex
            }
            target -= length;
        }
        next_vertex = numVertices()
        return vertex(numVertices() - 1) to next_vertex
    }

    // The inverse operation of GetSuffix/Interpolate.  Given a point on the
    // polyline, returns the ratio of the distance to the point from the
    // beginning of the polyline over the length of the polyline.  The return
    // value is always betwen 0 and 1 inclusive.  See GetSuffix() for the
    // meaning of "next_vertex".
    //
    // The polyline should not be empty.  If it has fewer than 2 vertices, the
    // return value is zero.
    fun unInterpolate(point: S2Point, nextVertex: Int): Double {
        Assertions.assertGT(numVertices(), 0);
        if (numVertices() < 2) {
            return 0.0
        }
        var length_sum = S1Angle()
        for (i in 1 until nextVertex) {
            length_sum += S1Angle(vertex(i - 1), vertex(i))
        }
        val length_to_point = length_sum + S1Angle(vertex(nextVertex - 1), point)
        for (i in nextVertex until numVertices()) {
            length_sum += S1Angle(vertex(i - 1), vertex(i));
        }
        // The ratio can be greater than 1.0 due to rounding errors or because the
        // point is not exactly on the polyline.
        return min(1.0, (length_to_point / length_sum).radians);
    }

    // Given a point, returns a point on the polyline that is closest to the given
    // point.  See GetSuffix() for the meaning of "next_vertex", which is chosen
    // here w.r.t. the projected point as opposed to the interpolated point in
    // GetSuffix().
    //
    // The polyline must be non-empty.
    fun project(point: S2Point): Pair<S2Point, Int /*next vertex*/> {
        Assertions.assertGT(numVertices(), 0);

        if (numVertices() == 1) {
            // If there is only one vertex, it is always closest to any given point.
            return vertex(0) to 1
        }

        // Initial value larger than any possible distance on the unit sphere.
        var min_distance = S1Angle.radians(10);
        var min_index = -1;

        // Find the line segment in the polyline that is closest to the point given.
        for (i in 1 until numVertices()) {
            val distance_to_segment = S2EdgeDistances.getDistance(point, vertex(i - 1), vertex(i));
            if (distance_to_segment < min_distance) {
                min_distance = distance_to_segment;
                min_index = i;
            }
        }
        Assertions.assertNE(min_index, -1);

        // Compute the point on the segment found that is closest to the point given.
        val closest_point = S2EdgeDistances.project(point, vertex(min_index - 1), vertex(min_index))

        val next_vertex = min_index + (if (closest_point == vertex(min_index)) 1 else 0)
        return closest_point to next_vertex
    }

    // Returns true if the point given is on the right hand side of the polyline,
    // using a naive definition of "right-hand-sideness" where the point is on
    // the RHS of the polyline iff the point is on the RHS of the line segment in
    // the polyline which it is closest to.
    //
    // The polyline must have at least 2 vertices.
    fun isOnRight(point: S2Point): Boolean {
        assertGE(numVertices(), 2);

        var (closest_point, next_vertex) = project(point)

        assertGE(next_vertex, 1)
        assertLE(next_vertex, numVertices())

        // If the closest point C is an interior vertex of the polyline, let B and D
        // be the previous and next vertices.  The given point P is on the right of
        // the polyline (locally) if B, P, D are ordered CCW around vertex C.
        if (closest_point == vertex(next_vertex - 1) && next_vertex > 1 && next_vertex < numVertices()) {
            if (point == vertex(next_vertex - 1))
                return false;  // Polyline vertices are not on the RHS.
            return S2Predicates.orderedCCW(vertex(next_vertex - 2), point, vertex(next_vertex), vertex(next_vertex - 1));
        }

        // Otherwise, the closest point C is incident to exactly one polyline edge.
        // We test the point P against that edge.
        if (next_vertex == numVertices())
            --next_vertex;

        return S2Predicates.sign(point, vertex(next_vertex), vertex(next_vertex - 1)) > 0
    }

    // Return true if this polyline intersects the given polyline. If the
    // polylines share a vertex they are considered to be intersecting. When a
    // polyline endpoint is the only intersection with the other polyline, the
    // function may return true or false arbitrarily.
    //
    // The running time is quadratic in the number of vertices.  (To intersect
    // polylines more efficiently, or compute the actual intersection geometry,
    // use S2BooleanOperation.)
    fun intersects(line: S2Polyline): Boolean {
        if (numVertices() <= 0 || line.numVertices() <= 0) {
            return false
        }

        if (!rectBound.intersects(line.rectBound)) {
            return false;
        }

        // TODO(ericv): Use S2ShapeIndex here.
        for (i in 1 until numVertices()) {
            val crosser = S2EdgeCrosser(vertex(i - 1), vertex(i), line.vertex(0))
            for (j in 1 until line.numVertices()) {
                if (crosser.crossingSign(line.vertex(j)) >= 0) {
                    return true;
                }
            }
        }
        return false;
    }

    // Reverse the order of the polyline vertices.
    fun reversed(): S2Polyline = S2Polyline(vertices.reversed())

    // Return a subsequence of vertex indices such that the polyline connecting
    // these vertices is never further than "tolerance" from the original
    // polyline.  Provided the first and last vertices are distinct, they are
    // always preserved; if they are not, the subsequence may contain only a
    // single index.
    //
    // Some useful properties of the algorithm:
    //
    //  - It runs in linear time.
    //
    //  - The output is always a valid polyline.  In particular, adjacent
    //    output vertices are never identical or antipodal.
    //
    //  - The method is not optimal, but it tends to produce 2-3% fewer
    //    vertices than the Douglas-Peucker algorithm with the same tolerance.
    //
    //  - The output is *parametrically* equivalent to the original polyline to
    //    within the given tolerance.  For example, if a polyline backtracks on
    //    itself and then proceeds onwards, the backtracking will be preserved
    //    (to within the given tolerance).  This is different than the
    //    Douglas-Peucker algorithm, which only guarantees geometric equivalence.
    //
    // See also S2PolylineSimplifier, which uses the same algorithm but is more
    // efficient and supports more features, and also S2Builder, which can
    // simplify polylines and polygons, supports snapping (e.g. to E7 lat/lng
    // coordinates or S2CellId centers), and can split polylines at intersection
    // points.
    fun subsampleVertices(tolerance: S1Angle): List<Int> {
        logger.trace { """
            | subsampleVertices(tolerance = $tolerance)
            | --------------------------------------------------------------
            | vertices: $vertices
        """.trimMargin() }
        val indices = mutableListOf<Int>()
        if (numVertices() == 0) return indices

        indices.add(0)
        val clampedTolerance = maxOf(tolerance, S1Angle.radians(0))
        var index = 0
        logger.trace { "Clamped tolerance = $clampedTolerance" }
        while (index + 1 < numVertices()) {
            logger.trace { "Process index $index" }
            val nextIndex = findEndVertex(this, clampedTolerance, index)
            logger.trace { "next index = $nextIndex" }
            // Don't create duplicate adjacent vertices.
            if (vertex(nextIndex) != vertex(index)) {
                indices.add(nextIndex)
            }
            index = nextIndex;
        }
        return indices
    }

    // Return true if two polylines are exactly the same.
    override fun equals(other: Any?): Boolean {
        if (this === other) return true
        if (other !is S2Polyline) return false

        if (vertices != other.vertices) return false

        return true
    }

    override fun hashCode(): Int {
        return vertices.hashCode()
    }


    // Return true if two polylines have the same number of vertices, and
    // corresponding vertex pairs are separated by no more than "max_error".
    // (For testing purposes.)
    fun approxEquals(b: S2Polyline, max_error: S1Angle = S1Angle.radians(1e-15)): Boolean {
        if (numVertices() != b.numVertices()) return false;
        for (offset in 0 until numVertices()) {
            if (!S2Point.approxEquals(vertex(offset), b.vertex(offset), max_error)) {
                return false;
            }
        }
        return true;
    }

    // Return true if "covered" is within "max_error" of a contiguous subpath of
    // this polyline over its entire length.  Specifically, this method returns
    // true if this polyline has parameterization a:[0,1] -> S^2, "covered" has
    // parameterization b:[0,1] -> S^2, and there is a non-decreasing function
    // f:[0,1] -> [0,1] such that distance(a(f(t)), b(t)) <= max_error for all t.
    //
    // You can think of this as testing whether it is possible to drive a car
    // along "covered" and a car along some subpath of this polyline such that no
    // car ever goes backward, and the cars are always within "max_error" of each
    // other.
    //
    // This function is well-defined for empty polylines:
    //    anything.covers(empty) = true
    //    empty.covers(nonempty) = false
    fun nearlyCovers(covered: S2Polyline, max_error: S1Angle): Boolean {
        // NOTE: This algorithm is described assuming that adjacent vertices in a
        // polyline are never at the same point.  That is, the ith and i+1th vertices
        // of a polyline are never at the same point in space.  The implementation
        // does not make this assumption.

        // DEFINITIONS:
        //   - edge "i" of a polyline is the edge from the ith to i+1th vertex.
        //   - covered_j is a polyline consisting of edges 0 through j of "covered."
        //   - this_i is a polyline consisting of edges 0 through i of this polyline.
        //
        // A search state is represented as an (int, int, bool) tuple, (i, j,
        // i_in_progress).  Using the "drive a car" analogy from the header comment, a
        // search state signifies that you can drive one car along "covered" from its
        // first vertex through a point on its jth edge, and another car along this
        // polyline from some point on or before its ith edge to a to a point on its
        // ith edge, such that no car ever goes backward, and the cars are always
        // within "max_error" of each other.  If i_in_progress is true, it means that
        // you can definitely drive along "covered" through the jth vertex (beginning
        // of the jth edge). Otherwise, you can definitely drive along "covered"
        // through the point on the jth edge of "covered" closest to the ith vertex of
        // this polyline.
        //
        // The algorithm begins by finding all edges of this polyline that are within
        // "max_error" of the first vertex of "covered," and adding search states
        // representing all of these possible starting states to the stack of
        // "pending" states.
        //
        // The algorithm proceeds by popping the next pending state,
        // (i,j,i_in_progress), off of the stack.  First it checks to see if that
        // state represents finding a valid covering of "covered" and returns true if
        // so.  Next, if the state represents reaching the end of this polyline
        // without finding a successful covering, the algorithm moves on to the next
        // state in the stack.  Otherwise, if state (i+1,j,false) is valid, it is
        // added to the stack of pending states.  Same for state (i,j+1,true).
        //
        // We need the stack because when "i" and "j" can both be incremented,
        // sometimes only one choice leads to a solution.  We use a set to keep track
        // of visited states to avoid duplicating work.  With the set, the worst-case
        // number of states examined is O(n+m) where n = this->num_vertices() and m =
        // covered.num_vertices().  Without it, the amount of work could be as high as
        // O((n*m)^2).  Using set, the running time is O((n*m) log (n*m)).
        //
        // TODO(user): Benchmark this, and see if the set is worth it.

        if (covered.numVertices() == 0) return true;
        if (numVertices() == 0) return false;

        val pending = mutableListOf<SearchState>()
        val done = mutableSetOf<SearchState>()

        // Find all possible starting states.
        var i = 0
        var next_i = nextDistinctVertex(this, 0)
        var next_next_i = 0
        while (next_i < numVertices()) {
            next_next_i = nextDistinctVertex(this, next_i)
            val closest_point = S2EdgeDistances.project(covered.vertex(0), this.vertex(i), vertex(next_i))

            // In order to avoid duplicate starting states, we exclude the end vertex
            // of each edge *except* for the last non-degenerate edge.
            if ((next_next_i == this.numVertices() || closest_point != this.vertex(next_i)) &&
                    S1Angle(closest_point, covered.vertex(0)) <= max_error) {
            pending.add(SearchState(i, 0, true))
        }
            i = next_i
            next_i = next_next_i
        }

        while (pending.isNotEmpty()) {
            val state = pending.removeLast()
            if (!done.add(state)) continue

            val next_i = nextDistinctVertex(this, state.i)
            val next_j = nextDistinctVertex(covered, state.j)
            if (next_j == covered.numVertices()) {
                return true;
            } else if (next_i == this.numVertices()) {
                continue;
            }

            val i_begin: S2Point
            val j_begin: S2Point
            if (state.i_in_progress) {
                j_begin = covered.vertex(state.j);
                i_begin = S2EdgeDistances.project(j_begin, this.vertex(state.i), this.vertex(next_i))
            } else {
                i_begin = this.vertex(state.i)
                j_begin = S2EdgeDistances.project(i_begin, covered.vertex(state.j), covered.vertex(next_j))
            }

            if (S2EdgeDistances.isEdgeBNearEdgeA(j_begin, covered.vertex(next_j), i_begin, this.vertex(next_i), max_error)) {
                pending.add(SearchState(next_i, state.j, false))
            }
            if (S2EdgeDistances.isEdgeBNearEdgeA(i_begin, this.vertex(next_i), j_begin, covered.vertex(next_j), max_error)) {
                pending.add(SearchState(state.i, next_j, true))
            }
        }
        return false
    }

    ////////////////////////////////////////////////////////////////////////
    // S2Region interface (see s2region.h for details):

    override fun clone(): S2Polyline {
        return S2Polyline(vertices.toList())
    }

    override val capBound: S2Cap
        get() = rectBound.capBound

    override val rectBound: S2LatLngRect
        get() {
            val bounder = S2LatLngRectBounder()
            for (i in 0 until numVertices()) {
                bounder.addPoint(vertex(i))
            }
            return bounder.getBound()
        }

    override fun contains(cell: S2Cell): Boolean = false

    override fun mayIntersect(cell: S2Cell): Boolean {
        if (numVertices() == 0) return false;

        // We only need to check whether the cell contains vertex 0 for correctness,
        // but these tests are cheap compared to edge crossings so we might as well
        // check all the vertices.
        for (i in 0 until numVertices()) {
            if (cell.contains(vertex(i))) return true
        }
        val cell_vertices = Array<S2Point>(4) { S2Point() }
        for (i in 0..3) {
            cell_vertices[i] = cell.getVertex(i)
        }
        for (j in 0..3) {
            val crosser = S2EdgeCrosser(cell_vertices[j], cell_vertices[(j + 1) and 3], vertex(0))
            for (i in 1 until numVertices()) {
                if (crosser.crossingSign(vertex(i)) >= 0) {
                    // There is a proper crossing, or two vertices were the same.
                    return true;
                }
            }
        }
        return false;
    }

    // Always return false, because "containment" is not numerically
    // well-defined except at the polyline vertices.
    override fun contains(p: S2Point): Boolean = false

    // Wrapper class for indexing a polyline (see S2ShapeIndex).  Once this
    // object is inserted into an S2ShapeIndex it is owned by that index, and
    // will be automatically deleted when no longer needed by the index.  Note
    // that this class does not take ownership of the polyline itself (see
    // OwningShape below).  You can also subtype this class to store additional
    // data (see S2Shape for details).
    class Shape(id: Int = 0, val polyline: S2Polyline) : S2Shape(id) {

        // S2Shape interface:

        override val numEdges: Int = max(0, polyline.numVertices() - 1)

        override fun edge(edgeId: Int): Edge = Edge(polyline.vertex(edgeId), polyline.vertex(edgeId + 1))

        override val dimension: Int = 1

        override fun getReferencePoint(): ReferencePoint = ReferencePoint(contained = false)

        override val numChains: Int = min(1, numEdges)

        override fun chain(chain_id: Int): Chain {
            Assertions.assertEQ(chain_id, 0)
            return Chain(0, numEdges)
        }

        override fun chainEdge(chainId: Int, offset: Int): Edge {
            Assertions.assertEQ(chainId, 0);
            return Edge(polyline.vertex(offset), polyline.vertex(offset + 1));
        }

        override fun chainPosition(edgeId: Int): ChainPosition = ChainPosition(0, edgeId)

        override val typeTag: TypeTag = kTypeTag

        companion object {
            val kTypeTag: TypeTag = 2U
        }
    };


    companion object {

        private val logger = KotlinLogging.logger(S2Polyline::class.java.name)

        fun fromLatLng(vertices: List<S2LatLng>): S2Polyline = S2Polyline(vertices.map { it.toPoint() })


        // Returns the length of the polyline.  Returns zero for polylines with fewer
        // than two vertices.
        fun getLength(polyline: List<S2Point>): S1Angle {
            val length = MutableS1Angle()
            for (i in 1 until polyline.size) {
                length.radians += S1Angle(polyline[i - 1], polyline[i]).radians
            }
            return length;
        }

        // Returns the true centroid of the polyline multiplied by the length of the
        // polyline (see s2centroids.h for details on centroids).  The result is not
        // unit length, so you may want to normalize it.
        //
        // Scaling by the polyline length makes it easy to compute the centroid of
        // several polylines (by simply adding up their centroids).
        //
        // CAVEAT: Returns S2Point() for degenerate polylines (e.g., AA).  [Note that
        // this answer is correct; the result of this function is a line integral over
        // the polyline, whose value is always zero if the polyline is degenerate.]
        fun getCentroid(polyline: List<S2Point>): S2Point {
            var centroid = S2Point()
            for (i in 1 until polyline.size) {
                centroid += S2Centroids.trueCentroid(polyline[i - 1], polyline[i])
            }
            return centroid;
        }

        // Given a polyline, a tolerance distance, and a start index, this function
        // returns the maximal end index such that the line segment between these two
        // vertices passes within "tolerance" of all interior vertices, in order.
        private fun findEndVertex(polyline: S2Polyline, tolerance: S1Angle, index: Int): Int {
            logger.trace { "findEndVertex(tolerance = $tolerance, index = $index)" }
            assertGE(tolerance.radians, 0.0)
            assertLT((index + 1), polyline.numVertices())

            // The basic idea is to keep track of the "pie wedge" of angles from the
            // starting vertex such that a ray from the starting vertex at that angle
            // will pass through the discs of radius "tolerance" centered around all
            // vertices processed so far.

            // First we define a "coordinate frame" for the tangent and normal spaces
            // at the starting vertex.  Essentially this means picking three
            // orthonormal vectors X,Y,Z such that X and Y span the tangent plane at
            // the starting vertex, and Z is "up".  We use the coordinate frame to
            // define a mapping from 3D direction vectors to a one-dimensional "ray
            // angle" in the range (-Pi, Pi].  The angle of a direction vector is
            // computed by transforming it into the X,Y,Z basis, and then calculating
            // atan2(y,x).  This mapping allows us to represent a wedge of angles as a
            // 1D interval.  Since the interval wraps around, we represent it as an
            // S1Interval, i.e. an interval on the unit circle.
            val origin = polyline.vertex(index)
            val frame = S2Point.getFrame(origin)

            // As we go along, we keep track of the current wedge of angles and the
            // distance to the last vertex (which must be non-decreasing).
            var current_wedge = S1Interval.full
            var last_distance = 0.0

            var targetIdx = index + 1
            while (targetIdx < polyline.numVertices()) {
                val candidate = polyline.vertex(index)
                val distance = origin.angle(candidate)

                // We don't allow simplification to create edges longer than 90 degrees,
                // to avoid numeric instability as lengths approach 180 degrees.  (We do
                // need to allow for original edges longer than 90 degrees, though.)
                if (distance > M_PI/2 && last_distance > 0) break

                // Vertices must be in increasing order along the ray, except for the
                // initial disc around the origin.
                if (distance < last_distance && last_distance > tolerance.radians) break
                last_distance = distance

                // Points that are within the tolerance distance of the origin do not
                // constrain the ray direction, so we can ignore them.
                if (distance <= tolerance.radians) {
                    targetIdx++
                    continue
                }

                // If the current wedge of angles does not contain the angle to this
                // vertex, then stop right now.  Note that the wedge of possible ray
                // angles is not necessarily empty yet, but we can't continue unless we
                // are willing to backtrack to the last vertex that was contained within
                // the wedge (since we don't create new vertices).  This would be more
                // complicated and also make the worst-case running time more than linear.
                val direction = S2Point.toFrame(frame, candidate)
                val center = atan2(direction.y(), direction.x())
                if (!current_wedge.contains(center)) break

                // To determine how this vertex constrains the possible ray angles,
                // consider the triangle ABC where A is the origin, B is the candidate
                // vertex, and C is one of the two tangent points between A and the
                // spherical cap of radius "tolerance" centered at B.  Then from the
                // spherical law of sines, sin(a)/sin(A) = sin(c)/sin(C), where "a" and
                // "c" are the lengths of the edges opposite A and C.  In our case C is a
                // 90 degree angle, therefore A = asin(sin(a) / sin(c)).  Angle A is the
                // half-angle of the allowable wedge.
                val half_angle = asin(sin(tolerance.radians) / sin(distance))
                val target = S1Interval.fromPoint(center).expanded(half_angle)
                current_wedge = current_wedge.intersection(target)
                assert { !current_wedge.isEmpty }
                targetIdx++
            }
            // We break out of the loop when we reach a vertex index that can't be
            // included in the line segment, so back up by one vertex.
            return targetIdx - 1;
        }

        // Return the first i > "index" such that the ith vertex of "pline" is not at
        // the same point as the "index"th vertex.  Returns pline.num_vertices() if
        // there is no such value of i.
        private fun nextDistinctVertex(pline: S2Polyline, index: Int): Int {
            var result = index
            val initial = pline.vertex(result)
            do {
                ++result
            } while (index < pline.numVertices() && pline.vertex(index) == initial)
            return index;
        }

        // This struct represents a search state in the NearlyCovers algorithm
        // below.  See the description of the algorithm for details.
        data class SearchState(val i: Int, val j: Int, val i_in_progress: Boolean) : Comparable<SearchState> {
            // This operator is needed for storing SearchStates in a set.  The ordering
            // chosen has no special meaning.
            override fun compareTo(other: SearchState): Int {
                if (i != other.i) return i.compareTo(other.i)
                if (j != other.j) return j.compareTo(other.j)
                return i_in_progress.compareTo(other.i_in_progress)
            }
        }

    }
}

