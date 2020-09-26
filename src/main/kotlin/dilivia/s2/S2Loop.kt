package dilivia.s2

import Matrix3x3
import com.google.common.geometry.S2.M_PI
import com.google.common.geometry.S2.M_PI_2
import dilivia.s2.shape.MutableS2ShapeIndex
import dilivia.s2.shape.S2Shape
import dilivia.s2.shape.S2ShapeUtil
import dilivia.s2.shape.TypeTag
import mu.KotlinLogging
import java.util.concurrent.atomic.AtomicInteger

// An S2Loop represents a simple spherical polygon.  It consists of a single
// chain of vertices where the first vertex is implicitly connected to the
// last. All loops are defined to have a CCW orientation, i.e. the interior of
// the loop is on the left side of the edges.  This implies that a clockwise
// loop enclosing a small area is interpreted to be a CCW loop enclosing a
// very large area.
//
// Loops are not allowed to have any duplicate vertices (whether adjacent or
// not).  Non-adjacent edges are not allowed to intersect, and furthermore edges
// of length 180 degrees are not allowed (i.e., adjacent vertices cannot be
// antipodal).  Loops must have at least 3 vertices (except for the empty and
// full loops discussed below).  Although these restrictions are not enforced
// in optimized code, you may get unexpected results if they are violated.
//
// There are two special loops: the "empty loop" contains no points, while the
// "full loop" contains all points.  These loops do not have any edges, but to
// preserve the invariant that every loop can be represented as a vertex
// chain, they are defined as having exactly one vertex each (see kEmpty and
// kFull).
//
// Point containment of loops is defined such that if the sphere is subdivided
// into faces (loops), every point is contained by exactly one face.  This
// implies that loops do not necessarily contain their vertices.
//
// Note: The reason that duplicate vertices and intersecting edges are not
// allowed is that they make it harder to define and implement loop
// relationships, e.g. whether one loop contains another.  If your data does
// not satisfy these restrictions, you can use S2Builder to normalize it.

// @constructor
// Convenience constructor to disable the automatic validity checking
// controlled by the --s2debug flag.  Example:
//
//   S2Loop* loop = new S2Loop(vertices, S2Debug::DISABLE)
//
// This is equivalent to:
//
//   S2Loop* loop = new S2Loop
//   loop->set_s2debug_override(S2Debug::DISABLE)
//   loop->Init(vertices)
//
// The main reason to use this constructor is if you intend to call
// IsValid() explicitly.  See set_s2debug_override() for details.

// @property depth
// The nesting depth, if this field belongs to an S2Polygon.  We define it
// here to optimize field packing.
// The depth of a loop is defined as its nesting level within its containing
// polygon.  "Outer shell" loops have depth 0, holes within those loops have
// depth 1, shells within those holes have depth 2, etc.  This field is only
// used by the S2Polygon implementation.
class S2Loop internal constructor(val vertices: List<S2Point>, val depth: Int = 0, check: Boolean, initOriginAndBound: Boolean = true) : S2Region {

    // "bound_" is a conservative bound on all points contained by this loop:
    // if A.Contains(P), then A.bound_.Contains(S2LatLng(P)).
    private var bound: S2LatLngRect = S2LatLngRect.empty

    // Since "bound_" is not exact, it is possible that a loop A contains
    // another loop B whose bounds are slightly larger.  "subregion_bound_"
    // has been expanded sufficiently to account for this error, i.e.
    // if A.Contains(B), then A.subregion_bound_.Contains(B.bound_).
    private var subregion_bound: S2LatLngRect = S2LatLngRect.empty

    // Spatial index for this loop.
    private val index: MutableS2ShapeIndex = MutableS2ShapeIndex()

    // In general we build the index the first time it is needed, but we make an
    // exception for Contains(S2Point) because this method has a simple brute
    // force implementation that is also relatively cheap.  For this one method
    // we keep track of the number of calls made and only build the index once
    // enough calls have been made that we think an index would be worthwhile.
    val unindexed_contains_calls = AtomicInteger(0)

    private var origin_inside = false  // Does the loop contain S2::Origin()?

    // Initialize a loop with given vertices.  The last vertex is implicitly
    // connected to the first.  All points should be unit length.  Loops must
    // have at least 3 vertices (except for the empty and full loops, see
    // kEmpty and kFull).  This method may be called multiple times.
    init {
        if (initOriginAndBound) {
            initOriginAndBound(check)
        }
    }

    private fun initOriginAndBound(check: Boolean) {
        if (numVertices() < 3) {
            // Check for the special empty and full loops (which have one vertex).
            if (!is_empty_or_full()) {
                origin_inside = false
                return  // Bail out without trying to access non-existent vertices.
            }
            // If the vertex is in the southern hemisphere then the loop is full,
            // otherwise it is empty.
            origin_inside = (vertex(0).z() < 0)
        } else {
            // Point containment testing is done by counting edge crossings starting
            // at a fixed point on the sphere (S2::Origin()).  Historically this was
            // important, but it is now no longer necessary, and it may be worthwhile
            // experimenting with using a loop vertex as the reference point.  In any
            // case, we need to know whether the reference point (S2::Origin) is
            // inside or outside the loop before we can construct the S2ShapeIndex.
            // We do this by first guessing that it is outside, and then seeing
            // whether we get the correct containment result for vertex 1.  If the
            // result is incorrect, the origin must be inside the loop.
            //
            // A loop with consecutive vertices A,B,C contains vertex B if and only if
            // the fixed vector R = S2::Ortho(B) is contained by the wedge ABC.  The
            // wedge is closed at A and open at C, i.e. the point B is inside the loop
            // if A=R but not if C=R.  This convention is required for compatibility
            // with S2::VertexCrossing.  (Note that we can't use S2::Origin()
            // as the fixed vector because of the possibility that B == S2::Origin().)
            //
            // TODO(ericv): Investigate using vertex(0) as the reference point.

            origin_inside = false  // Initialize before calling Contains().
            val v1_inside = S2Predicates.orderedCCW(S2Point.ortho(vertex(1)), vertex(0), vertex(2), vertex(1))
            // Note that Contains(S2Point) only does a bounds check once InitIndex()
            // has been called, so it doesn't matter that bound_ is undefined here.
            if (v1_inside != contains(vertex(1))) {
                origin_inside = true
            }
        }
        // We *must* call InitBound() before InitIndex(), because InitBound() calls
        // Contains(S2Point), and Contains(S2Point) does a bounds check whenever the
        // index is not fresh (i.e., the loop has been added to the index but the
        // index has not been updated yet).
        //
        // TODO(ericv): When fewer S2Loop methods depend on internal bounds checks,
        // consider computing the bound on demand as well.
        initBound()
        initIndex(check)
    }

    private fun initBound() {
        // Check for the special empty and full loops.
        if (is_empty_or_full()) {
            if (is_empty()) {
                subregion_bound = S2LatLngRect.empty
                bound = S2LatLngRect.empty
            } else {
                subregion_bound = S2LatLngRect.full
                bound = S2LatLngRect.full
            }
            return
        }

        // The bounding rectangle of a loop is not necessarily the same as the
        // bounding rectangle of its vertices.  First, the maximal latitude may be
        // attained along the interior of an edge.  Second, the loop may wrap
        // entirely around the sphere (e.g. a loop that defines two revolutions of a
        // candy-cane stripe).  Third, the loop may include one or both poles.
        // Note that a small clockwise loop near the equator contains both poles.

        val bounder = S2LatLngRectBounder()
        for (i in 0 .. numVertices()) {
            bounder.addPoint(vertex(i))
        }
        var b = bounder.getBound()
        if (contains(S2Point(0, 0, 1))) {
            b = S2LatLngRect(R1Interval(b.lat.lo, M_PI_2), S1Interval.full)
        }
        // If a loop contains the south pole, then either it wraps entirely
        // around the sphere (full longitude range), or it also contains the
        // north pole in which case b.lng().is_full() due to the test above.
        // Either way, we only need to do the south pole containment test if
        // b.lng().is_full().
        if (b.lng.isFull && contains(S2Point(0, 0, -1))) {
            b = S2LatLngRect(
                    R1Interval(-M_PI_2, b.lat.hi),
                    b.lng
            )
        }
        bound = b
        subregion_bound = S2LatLngRectBounder.expandForSubregions(bound)
    }

    private fun initIndex(check: Boolean) {
        index.add(Shape(loop = this))
     /*  TODO if (!FLAGS_s2loop_lazy_indexing) {
            index.ForceBuild()
        }*/
        if (check) {
            // Note that FLAGS_s2debug is false in optimized builds (by default).
            Assertions.assert { isValid() }
        }
    }

    // Default constructor.  The loop must be initialized by calling Init() or
    // Decode() before it is used.
    constructor() : this(emptyList(), check = true)

    // Convenience constructor that calls Init() with the given vertices.
    constructor(vertices: List<S2Point>) : this(vertices, check = true)

    // Construct a loop corresponding to the given cell.
    //
    // Note that the loop and cell *do not* contain exactly the same set of
    // points, because S2Loop and S2Cell have slightly different definitions of
    // point containment.  For example, an S2Cell vertex is contained by all
    // four neighboring S2Cells, but it is contained by exactly one of four
    // S2Loops constructed from those cells.  As another example, the S2Cell
    // coverings of "cell" and "S2Loop(cell)" will be different, because the
    // loop contains points on its boundary that actually belong to other cells
    // (i.e., the covering will include a layer of neighboring cells).
    constructor(cell: S2Cell): this((0..3).map { cell.getVertex(it) }, 0, true)

    // Returns true if this is a valid loop.  Note that validity is checked
    // automatically during initialization when --s2debug is enabled (true by
    // default in debug binaries).
    fun isValid(): Boolean {
        val error = findValidationError()
        if (!error.isOk()) {
            logger.error { error.message }
            return false
        }
        return true
    }

    // Returns true if this is *not* a valid loop and sets "error"
    // appropriately.  Otherwise returns false and leaves "error" unchanged.
    //
    // REQUIRES: error != nullptr
    fun findValidationError(): S2Error {
        var error = findValidationErrorNoIndex()
        if (error.isOk()) {
            error = S2ShapeUtil.findSelfIntersection(index)
        }
        return error
    }

    fun numVertices(): Int = vertices.size

    // For convenience, we make two entire copies of the vertex list available:
    // vertex(n..2*n-1) is mapped to vertex(0..n-1), where n == num_vertices().
    //
    // REQUIRES: 0 <= i < 2 * num_vertices()
    fun vertex(i: Int): S2Point {
        Assertions.assertGE(i, 0)
        Assertions.assertLT(i, 2 * numVertices())
        val j = i - numVertices()
        return vertices[if (j < 0) i else j]
    }

    // Like vertex(), but this method returns vertices in reverse order if the
    // loop represents a polygon hole.  For example, arguments 0, 1, 2 are
    // mapped to vertices n-1, n-2, n-3, where n == num_vertices().  This
    // ensures that the interior of the polygon is always to the left of the
    // vertex chain.
    //
    // REQUIRES: 0 <= i < 2 * num_vertices()
    fun oriented_vertex(i: Int): S2Point {
        Assertions.assertGE(i, 0)
        Assertions.assertLT(i, 2 * numVertices())
        var j = i - numVertices()
        if (j < 0) j = i
        if (is_hole()) j = numVertices() - 1 - j
        return vertices[j]
    }

    // Returns true if this is the special empty loop that contains no points.
    fun is_empty(): Boolean = is_empty_or_full() && !contains_origin()

    // Returns true if this is the special full loop that contains all points.
    fun is_full(): Boolean = is_empty_or_full() && contains_origin()

    // Returns true if this loop is either empty or full.
    fun is_empty_or_full(): Boolean = numVertices() == 1

    // Returns true if this loop represents a hole in its containing polygon.
    fun is_hole(): Boolean {
        return (depth and 1) != 0
    }

    // The sign of a loop is -1 if the loop represents a hole in its containing
    // polygon, and +1 otherwise.
    fun sign(): Int {
        return if (is_hole()) -1 else 1
    }

    // Returns true if the loop area is at most 2*Pi.  Degenerate loops are
    // handled consistently with s2pred::Sign(), i.e., if a loop can be
    // expressed as the union of degenerate or nearly-degenerate CCW triangles,
    // then it will always be considered normalized.
    fun isNormalized(): Boolean {
        // Optimization: if the longitude span is less than 180 degrees, then the
        // loop covers less than half the sphere and is therefore normalized.
        if (bound.lng.length < M_PI) return true

        return S2LoopMeasures.isNormalized(verticesSpan())
    }

    // Invert the loop if necessary so that the area enclosed by the loop is at
    // most 2*Pi.
    fun normalize(): Unit = TODO()

    // Reverse the order of the loop vertices, effectively complementing the
    // region represented by the loop.  For example, the loop ABCD (with edges
    // AB, BC, CD, DA) becomes the loop DCBA (with edges DC, CB, BA, AD).
    // Notice that the last edge is the same in both cases except that its
    // direction has been reversed.
    fun invert(): Unit = TODO()

    // Returns the area of the loop interior, i.e. the region on the left side of
    // the loop.  The return value is between 0 and 4*Pi.  (Note that the return
    // value is not affected by whether this loop is a "hole" or a "shell".)
    fun getArea(): Double = TODO()

    // Returns the true centroid of the loop multiplied by the area of the loop
    // (see s2centroids.h for details on centroids).  The result is not unit
    // length, so you may want to normalize it.  Also note that in general, the
    // centroid may not be contained by the loop.
    //
    // We prescale by the loop area for two reasons: (1) it is cheaper to
    // compute this way, and (2) it makes it easier to compute the centroid of
    // more complicated shapes (by splitting them into disjoint regions and
    // adding their centroids).
    //
    // Note that the return value is not affected by whether this loop is a
    // "hole" or a "shell".
    fun getCentroid(): S2Point = TODO()

    // Returns the geodesic curvature of the loop, defined as the sum of the turn
    // angles at each vertex (see S2::TurnAngle).  The result is positive if the
    // loop is counter-clockwise, negative if the loop is clockwise, and zero if
    // the loop is a great circle.  The geodesic curvature is equal to 2*Pi minus
    // the area of the loop.
    //
    // Degenerate and nearly-degenerate loops are handled consistently with
    // s2pred::Sign().  So for example, if a loop has zero area (i.e., it is a
    // very small CCW loop) then its geodesic curvature will always be positive.
    fun getCurvature(): Double = TODO()

    // Returns the maximum error in GetCurvature().  The return value is not
    // constant; it depends on the loop.
    fun getCurvatureMaxError(): Double = TODO()

    // Returns the distance from the given point to the loop interior.  If the
    // loop is empty, return S1Angle::Infinity().  "x" should be unit length.
    fun getDistance(x: S2Point): S1Angle = TODO()

    // Returns the distance from the given point to the loop boundary.  If the
    // loop is empty or full, return S1Angle::Infinity() (since the loop has no
    // boundary).  "x" should be unit length.
    fun getDistanceToBoundary(x: S2Point): S1Angle = TODO()

    // If the given point is contained by the loop, return it.  Otherwise return
    // the closest point on the loop boundary.  If the loop is empty, return the
    // input argument.  Note that the result may or may not be contained by the
    // loop.  "x" should be unit length.
    fun project(x: S2Point): S2Point = TODO()

    // Returns the closest point on the loop boundary to the given point.  If the
    // loop is empty or full, return the input argument (since the loop has no
    // boundary).  "x" should be unit length.
    fun projectToBoundary(x: S2Point): S2Point = TODO()

    // Returns true if the region contained by this loop is a superset of the
    // region contained by the given other loop.
    fun contains(b: S2Loop): Boolean = TODO()

    // Returns true if the region contained by this loop intersects the region
    // contained by the given other loop.
    fun intersects(b: S2Loop): Boolean = TODO()

    // Returns true if two loops have the same vertices in the same linear order
    // (i.e., cyclic rotations are not allowed).
    fun equals(b: S2Loop): Boolean = TODO()

    // Returns true if two loops have the same boundary.  This is true if and
    // only if the loops have the same vertices in the same cyclic order (i.e.,
    // the vertices may be cyclically rotated).  The empty and full loops are
    // considered to have different boundaries.
    fun boundaryEquals(b: S2Loop): Boolean = TODO()

    // Returns true if two loops have the same boundary except for vertex
    // perturbations.  More precisely, the vertices in the two loops must be in
    // the same cyclic order, and corresponding vertex pairs must be separated
    // by no more than "max_error".
    fun boundaryApproxEquals(b: S2Loop, max_error: S1Angle = S1Angle.radians(1e-15)): Boolean = TODO()

    // Returns true if the two loop boundaries are within "max_error" of each
    // other along their entire lengths.  The two loops may have different
    // numbers of vertices.  More precisely, this method returns true if the two
    // loops have parameterizations a:[0,1] -> S^2, b:[0,1] -> S^2 such that
    // distance(a(t), b(t)) <= max_error for all t.  You can   of this as
    // testing whether it is possible to drive two cars all the way around the
    // two loops such that no car ever goes backward and the cars are always
    // within "max_error" of each other.
    fun boundaryNear(b: S2Loop, max_error: S1Angle = S1Angle.radians(1e-15)): Boolean = TODO()

    // This method computes the oriented surface integral of some quantity f(x)
    // over the loop interior, given a function f_tri(A,B,C) that returns the
    // corresponding integral over the spherical triangle ABC.  Here "oriented
    // surface integral" means:
    //
    // (1) f_tri(A,B,C) must be the integral of f if ABC is counterclockwise,
    //     and the integral of -f if ABC is clockwise.
    //
    // (2) The result of this function is *either* the integral of f over the
    //     loop interior, or the integral of (-f) over the loop exterior.
    //
    // Note that there are at least two common situations where it easy to work
    // around property (2) above:
    //
    //  - If the integral of f over the entire sphere is zero, then it doesn't
    //    matter which case is returned because they are always equal.
    //
    //  - If f is non-negative, then it is easy to detect when the integral over
    //    the loop exterior has been returned, and the integral over the loop
    //    interior can be obtained by adding the integral of f over the entire
    //    unit sphere (a constant) to the result.
    //
    // Also requires that the default constructor for T must initialize the
    // value to zero.  (This is true for built-in types such as "double".)

    /* TODO() template <class T>
    T GetSurfaceIntegral(T f_tri(const S2Point&, const S2Point&, const S2Point&))
       = S2::GetSurfaceIntegral(vertices_span(), f_tri) */

    ////////////////////////////////////////////////////////////////////////
    // S2Region interface (see s2region.h for details):

    override fun clone(): S2Loop {
        val loop = S2Loop(this.vertices, depth, false, false)
        loop.origin_inside = origin_inside
        loop.bound = bound
        loop.subregion_bound = subregion_bound
        loop.unindexed_contains_calls.set(0)
        loop.initIndex(false)
        return  loop
    }

    // GetRectBound() returns essentially tight results, while GetCapBound()
    // might have a lot of extra padding.  Both bounds are conservative in that
    // if the loop contains a point P, then the bound contains P also.

    override val capBound: S2Cap
        get() = TODO("Not yet implemented")
    override val rectBound: S2LatLngRect
        get() = bound

    override fun contains(cell: S2Cell): Boolean {
        TODO("Not yet implemented")
    }

    override fun mayIntersect(cell: S2Cell): Boolean {
        TODO("Not yet implemented")
    }

    // The point 'p' does not need to be normalized.

    override fun contains(p: S2Point): Boolean {
        // NOTE(ericv): A bounds check slows down this function by about 50%.  It is
        // worthwhile only when it might allow us to delay building the index.
        if (!index.isFresh() && !bound.contains(p)) return false;

        // For small loops it is faster to just check all the crossings.  We also
        // use this method during loop initialization because InitOriginAndBound()
        // calls Contains() before InitIndex().  Otherwise, we keep track of the
        // number of calls to Contains() and only build the index when enough calls
        // have been made so that we think it is worth the effort.  Note that the
        // code below is structured so that if many calls are made in parallel only
        // one thread builds the index, while the rest continue using brute force
        // until the index is actually available.
        //
        // The constants below were tuned using the benchmarks.  It turns out that
        // building the index costs roughly 50x as much as Contains().  (The ratio
        // increases slowly from 46x with 64 points to 61x with 256k points.)  The
        // textbook approach to this problem would be to wait until the cumulative
        // time we would have saved with an index approximately equals the cost of
        // building the index, and then build it.  (This gives the optimal
        // competitive ratio of 2; look up "competitive algorithms" for details.)
        // We set the limit somewhat lower than this (20 rather than 50) because
        // building the index may be forced anyway by other API calls, and so we
        // want to err on the side of building it too early.

        val kMaxBruteForceVertices = 32;
        val kMaxUnindexedContainsCalls = 20;  // See notes above.
        if (index.numShapeIds() == 0 ||  // InitIndex() not called yet
                numVertices() <= kMaxBruteForceVertices ||
                (!index.isFresh() &&
                        unindexed_contains_calls.incrementAndGet() != kMaxUnindexedContainsCalls)) {
            return bruteForceContains(p);
        }
        // Otherwise we look up the S2ShapeIndex cell containing this point.  Note
        // the index is built automatically the first time an iterator is created.
        val it = index.iterator()
        if (!it.locate(p)) return false;
        return contains(it, p);
    }

    ////////////////////////////////////////////////////////////////////////
    // Methods intended primarily for use by the S2Polygon implementation:

    // Given two loops of a polygon, return true if A contains B.  This version
    // of Contains() is cheap because it does not test for edge intersections.
    // The loops must meet all the S2Polygon requirements; for example this
    // implies that their boundaries may not cross or have any shared edges
    // (although they may have shared vertices).
    fun containsNested(b: S2Loop): Boolean = TODO()

    // Returns +1 if A contains the boundary of B, -1 if A excludes the boundary
    // of B, and 0 if the boundaries of A and B cross.  Shared edges are handled
    // as follows: If XY is a shared edge, define Reversed(XY) to be true if XY
    // appears in opposite directions in A and B.  Then A contains XY if and
    // only if Reversed(XY) == B->is_hole().  (Intuitively, this checks whether
    // A contains a vanishingly small region extending from the boundary of B
    // toward the interior of the polygon to which loop B belongs.)
    //
    // This method is used for testing containment and intersection of
    // multi-loop polygons.  Note that this method is not symmetric, since the
    // result depends on the direction of loop A but not on the direction of
    // loop B (in the absence of shared edges).
    //
    // REQUIRES: neither loop is empty.
    // REQUIRES: if b->is_full(), then !b->is_hole().
    fun compareBoundary(b: S2Loop): Int = TODO()

    // Given two loops whose boundaries do not cross (see CompareBoundary),
    // return true if A contains the boundary of B.  If "reverse_b" is true, the
    // boundary of B is reversed first (which only affects the result when there
    // are shared edges).  This method is cheaper than CompareBoundary() because
    // it does not test for edge intersections.
    //
    // REQUIRES: neither loop is empty.
    // REQUIRES: if b->is_full(), then reverse_b == false.
    fun containsNonCrossingBoundary(b: S2Loop, reverse_b: S2Loop): Boolean = TODO()

    // Wrapper class for indexing a loop (see S2ShapeIndex).  Once this object
    // is inserted into an S2ShapeIndex it is owned by that index, and will be
    // automatically deleted when no longer needed by the index.  Note that this
    // class does not take ownership of the loop itself (see OwningShape below).
    // You can also subtype this class to store additional data (see S2Shape for
    // details).
    class Shape(id: Int = 0, val loop: S2Loop) : S2Shape(id) {

        // S2Shape interface:

        override val numEdges: Int = if (loop.is_empty_or_full()) 0 else loop.numVertices()

        override fun edge(edgeId: Int): Edge = Edge(loop.vertex(edgeId), loop.vertex(edgeId + 1))

        override val dimension: Int = 2

        override fun getReferencePoint(): ReferencePoint = ReferencePoint(S2Point.origin(), loop.contains_origin())

        override val numChains: Int
            get() = TODO("Not yet implemented")

        override fun chain(chain_id: Int): Chain {
            TODO("Not yet implemented")
        }

        override fun chainEdge(chainId: Int, offset: Int): Edge {
            Assertions.assertEQ(chainId, 0)
            return Edge(loop.vertex(offset), loop.vertex(offset + 1))
        }

        override fun chainPosition(edgeId: Int): ChainPosition = ChainPosition(0, edgeId)

        override val typeTag: TypeTag = kTypeTag

        companion object {
            val kTypeTag: TypeTag = 0U
        }

    }

    // List of vertices with operator[] maps index values in the range
    // [n, 2*n-1] to the range [0, n-1] by subtracting n (where n == size()).
    // In other words, two full copies of the vertex array are available.  (This
    // is a compromise between convenience and efficiency, since computing the
    // index modulo "n" is surprisingly expensive.)
    //
    // This property is useful for implementing algorithms where the elements of
    // the span represent the vertices of a loop.
    class S2PointLoopSpan(val points: List<S2Point> = emptyList()): List<S2Point> by points {

        // Like operator[], but allows index values in the range [0, 2*size()-1]
        // where each index i >= size() is mapped to i - size().
        override operator fun get(index: Int): S2Point {
            Assertions.assertGE(index, 0)
            Assertions.assertLT(index, 2 * points.size)
            val j = index - points.size
            return points[if (j < 0) index else j]
        }

        fun rotated(startIndex: Int): S2PointLoopSpan {
            Assertions.assertGE(startIndex, 0)
            Assertions.assertLT(startIndex, points.size)
            val vertices = mutableListOf<S2Point>()
            for (i in startIndex until (startIndex+points.size)) {
                vertices.add(this[i])
            }
            return S2PointLoopSpan(vertices)
        }

        override fun toString(): String = points.toString()
    }

    // Returns an S2PointLoopSpan containing the loop vertices, for use with the
    // functions defined in s2loop_measures.h.
    fun verticesSpan(): S2PointLoopSpan = S2PointLoopSpan(vertices)

    // Returns true if this loop contains S2::Origin().
    private fun contains_origin(): Boolean = origin_inside

    // A version of Contains(S2Point) that does not use the S2ShapeIndex.
    // Used by the S2Polygon implementation.
    private fun bruteForceContains(p: S2Point): Boolean = TODO()

    // Like FindValidationError(), but skips any checks that would require
    // building the S2ShapeIndex (i.e., self-intersection tests).  This is used
    // by the S2Polygon implementation, which uses its own index to check for
    // loop self-intersections.
    private fun findValidationErrorNoIndex(): S2Error {
        // subregion_bound_ must be at least as large as bound_.  (This is an
        // internal consistency check rather than a test of client data.)
        Assertions.assert { subregion_bound.contains(bound) }

        // All vertices must be unit length.  (Unfortunately this check happens too
        // late in debug mode, because S2Loop construction calls s2pred::Sign which
        // expects vertices to be unit length.  But it is still a useful check in
        // optimized builds.)
        for (i in 0 until numVertices()) {
            if (!S2Point.isUnitLength(vertex(i))) {
                return S2Error(S2Error.NOT_UNIT_LENGTH,  "Vertex %d is not unit length".format(i))
            }
        }
        // Loops must have at least 3 vertices (except for the empty and full loops).
        if (numVertices() < 3) {
            if (is_empty_or_full()) {
                return S2Error(S2Error.OK);  // Skip remaining tests.
            }
            return S2Error(S2Error.LOOP_NOT_ENOUGH_VERTICES,
            "Non-empty, non-full loops must have at least 3 vertices")
        }
        // Loops are not allowed to have any duplicate vertices or edge crossings.
        // We split this check into two parts.  First we check that no edge is
        // degenerate (identical endpoints).  Then we check that there are no
        // intersections between non-adjacent edges (including at vertices).  The
        // second part needs the S2ShapeIndex, so it does not fall within the scope
        // of this method.
        for (i in 0 until numVertices()) {
            if (vertex(i) == vertex(i+1)) {
                return S2Error(S2Error.DUPLICATE_VERTICES,
                    "Edge %d is degenerate (duplicate vertex)".format(i))
            }
            if (vertex(i) == -vertex(i + 1)) {
                return S2Error(S2Error.ANTIPODAL_VERTICES,
                    "Vertices %d and %d are antipodal".format(i, (i + 1) % numVertices()))
            }
        }
        return S2Error(S2Error.OK)
    }

    // Converts the loop vertices to the S2XYZFaceSiTi format and store the result
    // in the given array, which must be large enough to store all the vertices.
    fun getXYZFaceSiTiVertices(): List<S2XYZFaceSiTi> = TODO()

  // Given an iterator that is already positioned at the S2ShapeIndexCell
  // containing "p", returns Contains(p).
  fun contains(it: MutableS2ShapeIndex.Iterator, p: S2Point): Boolean = TODO()
/*
  // Returns true if the loop boundary intersects "target".  It may also
  // return true when the loop boundary does not intersect "target" but
  // some edge comes within the worst-case error tolerance.
  //
  // REQUIRES: it.id().contains(target.id())
  // [This condition is true whenever it.Locate(target) returns INDEXED.]
  bool BoundaryApproxIntersects(const MutableS2ShapeIndex::Iterator& it,
                                const S2Cell& target) const

  // Returns an index "first" and a direction "dir" such that the vertex
  // sequence (first, first + dir, ..., first + (n - 1) * dir) does not change
  // when the loop vertex order is rotated or reversed.  This allows the loop
  // vertices to be traversed in a canonical order.
  S2::LoopOrder GetCanonicalLoopOrder() const
*/
  // Returns the index of a vertex at point "p", or -1 if not found.
  // The return value is in the range 1..num_vertices_ if found.
  private fun findVertex(p: S2Point): Int {
    if (numVertices() < 10) {
        // Exhaustive search.  Return value must be in the range [1..N].
        for (i in 1 .. numVertices()) {
            if (vertex(i) == p) return i
        }
        return -1
    }
    TODO("Implementation S2ShapeIndex")
    /*MutableS2ShapeIndex::Iterator it(&index_)
    if (!it.Locate(p)) return -1

    const S2ClippedShape& a_clipped = it.cell().clipped(0)
    for (int i = a_clipped.num_edges() - 1; i >= 0; --i) {
        int ai = a_clipped.edge(i)
        // Return value must be in the range [1..N].
        if (vertex(ai) == p) return (ai == 0) ? num_vertices() : ai
        if (vertex(ai+1) == p) return ai+1
    }
    return -1;*/
    }

    // When the loop is modified (Invert(), or Init() called again) then the
    // indexing structures need to be cleared since they become invalid.
    private fun clearIndex(): Unit = TODO()

    companion object {

        private val logger = KotlinLogging.logger(S2Loop::class.java.name)

        // The single vertex in the "empty loop" vertex chain.
        val kEmptyVertex = S2Point(0, 0, 1)

        // The single vertex in the "full loop" vertex chain.
        val kFullVertex = S2Point(0, 0, -1)

        // Any single-vertex loop is interpreted as being either the empty loop or the
        // full loop, depending on whether the vertex is in the northern or southern
        // hemisphere respectively.
        // A special vertex chain of length 1 that creates an empty loop (i.e., a
        // loop with no edges that contains no points).  Example usage:
        //
        //    S2Loop empty(S2Loop::kEmpty())
        //
        // The loop may be safely encoded lossily (e.g. by snapping it to an S2Cell
        // center) as long as its position does not move by 90 degrees or more.
        val kEmpty = listOf(kEmptyVertex)

        // A special vertex chain of length 1 that creates a full loop (i.e., a loop
        // with no edges that contains all points).  See kEmpty() for details.
        val kFull = listOf(kFullVertex)

        // Constructs a regular polygon with the given number of vertices, all
        // located on a circle of the specified radius around "center".  The radius
        // is the actual distance from "center" to each vertex.
        fun makeRegularLoop(center: S2Point, radius: S2Point, num_vertices: Int): S2Loop = TODO()

        // Like the function above, but this version constructs a loop centered
        // around the z-axis of the given coordinate frame, with the first vertex in
        // the direction of the positive x-axis.  (This allows the loop to be
        // rotated for testing purposes.)
        fun makeRegularLoop(frame: Matrix3x3, radius: S1Angle, num_vertices: Int): S2Loop = TODO()

        // This method checks all edges of loop A for intersection against all edges
        // of loop B.  If there is any shared vertex, the wedges centered at this
        // vertex are sent to "relation".
        //
        // If the two loop boundaries cross, this method is guaranteed to return
        // true.  It also returns true in certain cases if the loop relationship is
        // equivalent to crossing.  For example, if the relation is Contains() and a
        // point P is found such that B contains P but A does not contain P, this
        // method will return true to indicate that the result is the same as though
        // a pair of crossing edges were found (since Contains() returns false in
        // both cases).
        //
        // See Contains(), Intersects() and CompareBoundary() for the three uses of
        // this function.
        fun hasCrossingRelation(a: S2Loop, b: S2Loop, relation: LoopRelation): Boolean = TODO()
    }
}


// LoopRelation is an abstract class that defines a relationship between two
// loops (Contains, Intersects, or CompareBoundary).
abstract class LoopRelation {

    // Optionally, a_target() and b_target() can specify an early-exit condition
    // for the loop relation.  If any point P is found such that
    //
    //   A.Contains(P) == a_crossing_target() &&
    //   B.Contains(P) == b_crossing_target()
    //
    // then the loop relation is assumed to be the same as if a pair of crossing
    // edges were found.  For example, the Contains() relation has
    //
    //   a_crossing_target() == 0
    //   b_crossing_target() == 1
    //
    // because if A.Contains(P) == 0 (false) and B.Contains(P) == 1 (true) for
    // any point P, then it is equivalent to finding an edge crossing (i.e.,
    // since Contains() returns false in both cases).
    //
    // Loop relations that do not have an early-exit condition of this form
    // should return -1 for both crossing targets.
    abstract fun a_crossing_target(): Int
    abstract fun b_crossing_target(): Int

    // Given a vertex "ab1" that is shared between the two loops, return true if
    // the two associated wedges (a0, ab1, b2) and (b0, ab1, b2) are equivalent
    // to an edge crossing.  The loop relation is also allowed to maintain its
    // own internal state, and can return true if it observes any sequence of
    // wedges that are equivalent to an edge crossing.
    fun wedgesCross(a0: S2Point, ab1: S2Point, a2: S2Point, b0: S2Point, b2: S2Point): Boolean = TODO()

}
