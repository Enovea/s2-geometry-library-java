/**
 * This project is a kotlin port of the Google s2 geometry library (Copyright 2005 Google Inc. All Rights Reserved.):
 *                                 https://github.com/google/s2geometry.git
 *
 * Copyright © 2020 Dilivia (contact@dilivia.com)
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
package dilivia.s2.region

import Matrix3x3
import com.google.common.geometry.S2.M_PI
import com.google.common.geometry.S2.M_PI_2
import dilivia.s2.Assertions
import dilivia.s2.Assertions.assert
import dilivia.s2.Assertions.assertEQ
import dilivia.s2.R1Interval
import dilivia.s2.S1Angle
import dilivia.s2.S1Interval
import dilivia.s2.S2CellId
import dilivia.s2.S2EdgeClipping
import dilivia.s2.S2EdgeCrosser
import dilivia.s2.S2EdgeDistances
import dilivia.s2.S2Error
import dilivia.s2.S2LatLngRect
import dilivia.s2.S2LatLngRectBounder
import dilivia.s2.S2PaddedCell
import dilivia.s2.S2Point
import dilivia.s2.S2Predicates
import dilivia.s2.S2WedgeRelations
import dilivia.s2.coords.S2Coords
import dilivia.s2.coords.S2XYZFaceSiTi
import dilivia.s2.index.S2CrossingEdgePairsScanner
import dilivia.s2.index.S2CrossingEdgeQuery
import dilivia.s2.index.MutableS2ShapeIndex
import dilivia.s2.index.RangeIterator
import dilivia.s2.shape.S2ClippedShape
import dilivia.s2.shape.S2Shape
import dilivia.s2.index.S2ShapeIndex
import dilivia.s2.index.S2ShapeIndexCell
import dilivia.s2.shape.Edge
import dilivia.s2.shape.TypeTag
import mu.KotlinLogging
import java.util.concurrent.atomic.AtomicInteger
import kotlin.math.cos
import kotlin.math.sin

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
class S2Loop internal constructor(vertices: List<S2Point>, val depth: Int = 0, check: Boolean, initOriginAndBound: Boolean = true) : S2Region {

    private val vertices: MutableList<S2Point> = vertices.toMutableList()

    // "bound_" is a conservative bound on all points contained by this loop:
    // if A.Contains(P), then A.bound_.Contains(S2LatLng(P)).
    private var bound: S2LatLngRect = S2LatLngRect.empty

    // Since "bound_" is not exact, it is possible that a loop A contains
    // another loop B whose bounds are slightly larger.  "subregion_bound_"
    // has been expanded sufficiently to account for this error, i.e.
    // if A.Contains(B), then A.subregion_bound_.Contains(B.bound_).
    private var subregionBound: S2LatLngRect = S2LatLngRect.empty

    // Spatial index for this loop.
    val index: MutableS2ShapeIndex = MutableS2ShapeIndex()

    // In general we build the index the first time it is needed, but we make an
    // exception for Contains(S2Point) because this method has a simple brute
    // force implementation that is also relatively cheap.  For this one method
    // we keep track of the number of calls made and only build the index once
    // enough calls have been made that we think an index would be worthwhile.
    val unindexedContainsCalls = AtomicInteger(0)

    private var originInside = false  // Does the loop contain S2::Origin()?

    // Initialize a loop with given vertices.  The last vertex is implicitly
    // connected to the first.  All points should be unit length.  Loops must
    // have at least 3 vertices (except for the empty and full loops, see
    // kEmpty and kFull).  This method may be called multiple times.
    init {
        logger.trace { "Create loop: $vertices" }
        if (initOriginAndBound) {
            initOriginAndBound(check)
        }
    }

    private fun initOriginAndBound(check: Boolean) {
        if (numVertices() < 3) {
            // Check for the special empty and full loops (which have one vertex).
            if (!isEmptyOrFull()) {
                originInside = false
                return  // Bail out without trying to access non-existent vertices.
            }
            // If the vertex is in the southern hemisphere then the loop is full,
            // otherwise it is empty.
            originInside = (vertex(0).z() < 0)
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

            originInside = false  // Initialize before calling Contains().
            val v1Inside = S2Predicates.orderedCCW(S2Point.ortho(vertex(1)), vertex(0), vertex(2), vertex(1))
            // Note that Contains(S2Point) only does a bounds check once InitIndex()
            // has been called, so it doesn't matter that bound_ is undefined here.
            if (v1Inside != contains(vertex(1))) {
                originInside = true
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
        if (isEmptyOrFull()) {
            if (isEmpty()) {
                subregionBound = S2LatLngRect.empty
                bound = S2LatLngRect.empty
                logger.trace { "initBound(): Loop is empty => bound = $bound, subregion = $subregionBound" }
            } else {
                subregionBound = S2LatLngRect.full
                bound = S2LatLngRect.full
                logger.trace { "initBound(): Loop is full => bound = $bound, subregion = $subregionBound" }
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
        for (i in 0..numVertices()) {
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
        subregionBound = S2LatLngRectBounder.expandForSubregions(bound)
    }

    private fun initIndex(check: Boolean) {
        index.add(Shape(loop = this))
        /*  TODO if (!FLAGS_s2loop_lazy_indexing) {
               index.ForceBuild()
           }*/
        if (check) {
            // Note that FLAGS_s2debug is false in optimized builds (by default).
            assert { isValid() }
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
    constructor(cell: S2Cell) : this((0..3).map { cell.getVertex(it) }, 0, true)

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
            error = S2CrossingEdgePairsScanner.findSelfIntersection(index)
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
    fun orientedVertex(i: Int): S2Point {
        Assertions.assertGE(i, 0)
        Assertions.assertLT(i, 2 * numVertices())
        var j = i - numVertices()
        if (j < 0) j = i
        if (isHole()) j = numVertices() - 1 - j
        return vertices[j]
    }

    // Returns true if this is the special empty loop that contains no points.
    fun isEmpty(): Boolean = isEmptyOrFull() && !containsOrigin()

    // Returns true if this is the special full loop that contains all points.
    fun isFull(): Boolean = isEmptyOrFull() && containsOrigin()

    // Returns true if this loop is either empty or full.
    fun isEmptyOrFull(): Boolean = numVertices() == 1

    // Returns true if this loop represents a hole in its containing polygon.
    fun isHole(): Boolean = (depth and 1) != 0

    // The sign of a loop is -1 if the loop represents a hole in its containing
    // polygon, and +1 otherwise.
    fun sign(): Int = if (isHole()) -1 else 1

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
    fun normalize(): Unit {
        if (!isNormalized()) invert()
        assert { isNormalized() }
    }

    // Reverse the order of the loop vertices, effectively complementing the
    // region represented by the loop.  For example, the loop ABCD (with edges
    // AB, BC, CD, DA) becomes the loop DCBA (with edges DC, CB, BA, AD).
    // Notice that the last edge is the same in both cases except that its
    // direction has been reversed.
    fun invert() {
        clearIndex()
        if (isEmptyOrFull()) {
            vertices[0] = if (isFull()) kEmptyVertex else kFullVertex
        } else {
            vertices.reverse()
        }
        // origin_inside_ must be set correctly before building the S2ShapeIndex.
        originInside = originInside xor true
        if (bound.lat.lo > -M_PI_2 && bound.lat.hi < M_PI_2) {
            // The complement of this loop contains both poles.
            subregionBound = S2LatLngRect.full
            bound = S2LatLngRect.full
        } else {
            initBound()
        }
        initIndex(check = true)
    }

    // Returns the area of the loop interior, i.e. the region on the left side of
    // the loop.  The return value is between 0 and 4*Pi.  (Note that the return
    // value is not affected by whether this loop is a "hole" or a "shell".)
    fun getArea(): Double {
        // S2Loop has its own convention for empty and full loops.
        if (isEmptyOrFull()) {
            return if (containsOrigin()) (4 * M_PI) else 0.0
        }
        return S2LoopMeasures.getArea(verticesSpan())
    }

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
    fun getCentroid(): S2Point = S2LoopMeasures.getCentroid(verticesSpan())

    // Returns the geodesic curvature of the loop, defined as the sum of the turn
    // angles at each vertex (see S2::TurnAngle).  The result is positive if the
    // loop is counter-clockwise, negative if the loop is clockwise, and zero if
    // the loop is a great circle.  The geodesic curvature is equal to 2*Pi minus
    // the area of the loop.
    //
    // Degenerate and nearly-degenerate loops are handled consistently with
    // s2pred::Sign().  So for example, if a loop has zero area (i.e., it is a
    // very small CCW loop) then its geodesic curvature will always be positive.
    fun getCurvature(): Double {
        // S2Loop has its own convention for empty and full loops.  For such loops,
        // we return the limit value as the area approaches 0 or 4*Pi respectively.
        if (isEmptyOrFull()) {
            return if (containsOrigin()) (-2 * M_PI) else (2 * M_PI)
        }
        return S2LoopMeasures.getCurvature(verticesSpan())
    }

    // Returns the maximum error in GetCurvature().  The return value is not
    // constant; it depends on the loop.
    fun getCurvatureMaxError(): Double = S2LoopMeasures.getCurvatureMaxError(verticesSpan())

    // Returns the distance from the given point to the loop interior.  If the
    // loop is empty, return S1Angle::Infinity().  "x" should be unit length.
    fun getDistance(x: S2Point): S1Angle {
        // Note that S2Loop::Contains(S2Point) is slightly more efficient than the
        // generic version used by S2ClosestEdgeQuery.
        if (contains(x)) return S1Angle.zero
        return getDistanceToBoundary(x)
    }

    // Returns the distance from the given point to the loop boundary.  If the
    // loop is empty or full, return S1Angle::Infinity() (since the loop has no
    // boundary).  "x" should be unit length.
    fun getDistanceToBoundary(x: S2Point): S1Angle {
        TODO()
        /*S2ClosestEdgeQuery::Options options;
        options.set_include_interiors(false);
        S2ClosestEdgeQuery::PointTarget t(x);
        return S2ClosestEdgeQuery(&index_, options).GetDistance(&t).ToAngle();*/
    }

    // If the given point is contained by the loop, return it.  Otherwise return
    // the closest point on the loop boundary.  If the loop is empty, return the
    // input argument.  Note that the result may or may not be contained by the
    // loop.  "x" should be unit length.
    fun project(x: S2Point): S2Point {
        if (contains(x)) return x;
        return projectToBoundary(x);
    }

    // Returns the closest point on the loop boundary to the given point.  If the
    // loop is empty or full, return the input argument (since the loop has no
    // boundary).  "x" should be unit length.
    fun projectToBoundary(x: S2Point): S2Point {
        TODO()
        /*S2ClosestEdgeQuery::Options options;
        options.set_include_interiors(false);
        S2ClosestEdgeQuery q(&index_, options);
        S2ClosestEdgeQuery::PointTarget target(x);
        S2ClosestEdgeQuery::Result edge = q.FindClosestEdge(&target);
        return q.Project(x, edge);*/
    }

    // Returns true if the region contained by this loop is a superset of the
    // region contained by the given other loop.
    fun contains(b: S2Loop): Boolean {
        // For this loop A to contains the given loop B, all of the following must
        // be true:
        //
        //  (1) There are no edge crossings between A and B except at vertices.
        //
        //  (2) At every vertex that is shared between A and B, the local edge
        //      ordering implies that A contains B.
        //
        //  (3) If there are no shared vertices, then A must contain a vertex of B
        //      and B must not contain a vertex of A.  (An arbitrary vertex may be
        //      chosen in each case.)
        //
        // The second part of (3) is necessary to detect the case of two loops whose
        // union is the entire sphere, i.e. two loops that contains each other's
        // boundaries but not each other's interiors.
        if (!subregionBound.contains(b.bound)) return false;

        // Special cases to handle either loop being empty or full.
        if (isEmptyOrFull() || b.isEmptyOrFull()) {
            return isFull() || b.isEmpty()
        }

        // Check whether there are any edge crossings, and also check the loop
        // relationship at any shared vertices.
        val relation = ContainsRelation()
        if (hasCrossingRelation(this, b, relation)) return false

        // There are no crossings, and if there are any shared vertices then A
        // contains B locally at each shared vertex.
        if (relation.found_shared_vertex) return true

        // Since there are no edge intersections or shared vertices, we just need to
        // test condition (3) above.  We can skip this test if we discovered that A
        // contains at least one point of B while checking for edge crossings.
        if (!contains(b.vertex(0))) return false

        // We still need to check whether (A union B) is the entire sphere.
        // Normally this check is very cheap due to the bounding box precondition.
        if ((b.subregionBound.contains(bound) ||
                        b.bound.union(bound).isFull) && b.contains(vertex(0))) {
            return false
        }
        return true
    }

    // Returns true if the region contained by this loop intersects the region
    // contained by the given other loop.
    fun intersects(b: S2Loop): Boolean {
        // a->Intersects(b) if and only if !a->Complement()->Contains(b).
        // This code is similar to Contains(), but is optimized for the case
        // where both loops enclose less than half of the sphere.
        if (!bound.intersects(b.bound)) return false

        // Check whether there are any edge crossings, and also check the loop
        // relationship at any shared vertices.
        val relation = IntersectsRelation()
        if (hasCrossingRelation(this, b, relation)) return true
        if (relation.found_shared_vertex) return false

        // Since there are no edge intersections or shared vertices, the loops
        // intersect only if A contains B, B contains A, or the two loops contain
        // each other's boundaries.  These checks are usually cheap because of the
        // bounding box preconditions.  Note that neither loop is empty (because of
        // the bounding box check above), so it is safe to access vertex(0).

        // Check whether A contains B, or A and B contain each other's boundaries.
        // (Note that A contains all the vertices of B in either case.)
        if (subregionBound.contains(b.bound) || bound.union(b.bound).isFull) {
            if (contains(b.vertex(0))) return true
        }
        // Check whether B contains A.
        if (b.subregionBound.contains(bound)) {
            if (b.contains(vertex(0))) return true;
        }
        return false
    }

    // Returns true if two loops have the same vertices in the same linear order
    // (i.e., cyclic rotations are not allowed).
    override fun equals(other: Any?): Boolean {
        if (this === other) return true
        if (other !is S2Loop) return false

        if (vertices != other.vertices) return false

        return true
    }

    override fun hashCode(): Int {
        return vertices.hashCode()
    }

    // Returns true if two loops have the same boundary.  This is true if and
    // only if the loops have the same vertices in the same cyclic order (i.e.,
    // the vertices may be cyclically rotated).  The empty and full loops are
    // considered to have different boundaries.
    fun boundaryEquals(b: S2Loop): Boolean {
        if (numVertices() != b.numVertices()) return false

        // Special case to handle empty or full loops.  Since they have the same
        // number of vertices, if one loop is empty/full then so is the other.
        if (isEmptyOrFull()) return isEmpty() == b.isEmpty()

        for (offset in 0 until numVertices()) {
            if (vertex(offset) == b.vertex(0)) {
                // There is at most one starting offset since loop vertices are unique.
                for (i in 0 until numVertices()) {
                    if (vertex(i + offset) != b.vertex(i)) return false
                }
                return true;
            }
        }
        return false;
    }

    // Returns true if two loops have the same boundary except for vertex
    // perturbations.  More precisely, the vertices in the two loops must be in
    // the same cyclic order, and corresponding vertex pairs must be separated
    // by no more than "max_error".
    fun boundaryApproxEquals(b: S2Loop, max_error: S1Angle = S1Angle.radians(1e-15)): Boolean {
        if (numVertices() != b.numVertices()) return false

        // Special case to handle empty or full loops.  Since they have the same
        // number of vertices, if one loop is empty/full then so is the other.
        if (isEmptyOrFull()) return isEmpty() == b.isEmpty()

        for (offset in 0 until numVertices()) {
            if (S2Point.approxEquals(vertex(offset), b.vertex(0), max_error)) {
                var success = true
                for (i in 0 until numVertices()) {
                    if (!S2Point.approxEquals(vertex(i + offset), b.vertex(i), max_error)) {
                        success = false
                        break
                    }
                }
                if (success) return true
                // Otherwise continue looping.  There may be more than one candidate
                // starting offset since vertices are only matched approximately.
            }
        }
        return false
    }

    // Returns true if the two loop boundaries are within "max_error" of each
    // other along their entire lengths.  The two loops may have different
    // numbers of vertices.  More precisely, this method returns true if the two
    // loops have parameterizations a:[0,1] -> S^2, b:[0,1] -> S^2 such that
    // distance(a(t), b(t)) <= max_error for all t.  You can   of this as
    // testing whether it is possible to drive two cars all the way around the
    // two loops such that no car ever goes backward and the cars are always
    // within "max_error" of each other.
    fun boundaryNear(b: S2Loop, max_error: S1Angle = S1Angle.radians(1e-15)): Boolean {
        // Special case to handle empty or full loops.
        if (isEmptyOrFull() || b.isEmptyOrFull()) {
            return (isEmpty() && b.isEmpty()) || (isFull() && b.isFull())
        }

        for (a_offset in 0 until numVertices()) {
            if (matchBoundaries(this, b, a_offset, max_error)) return true;
        }
        return false;
    }

    ////////////////////////////////////////////////////////////////////////
    // S2Region interface (see s2region.h for details):

    public override fun clone(): S2Loop {
        val loop = S2Loop(this.vertices.toMutableList(), depth, check = false, initOriginAndBound = false)
        loop.originInside = originInside
        loop.bound = bound
        loop.subregionBound = subregionBound
        loop.unindexedContainsCalls.set(0)
        loop.initIndex(false)
        return loop
    }

    // GetRectBound() returns essentially tight results, while GetCapBound()
    // might have a lot of extra padding.  Both bounds are conservative in that
    // if the loop contains a point P, then the bound contains P also.

    override val capBound: S2Cap
        get() = bound.capBound

    override val rectBound: S2LatLngRect
        get() = bound

    override fun contains(cell: S2Cell): Boolean {
        val iterator = index.iterator()
        val relation = iterator.locate(cell.id())

        // If "target" is disjoint from all index cells, it is not contained.
        // Similarly, if "target" is subdivided into one or more index cells then it
        // is not contained, since index cells are subdivided only if they (nearly)
        // intersect a sufficient number of edges.  (But note that if "target" itself
        // is an index cell then it may be contained, since it could be a cell with
        // no edges in the loop interior.)
        if (relation != S2ShapeIndex.CellRelation.INDEXED) return false

        // Otherwise check if any edges intersect "target".
        if (boundaryApproxIntersects(iterator, cell)) return false

        // Otherwise check if the loop contains the center of "target".
        return contains(iterator, cell.getCenter())
    }

    override fun mayIntersect(cell: S2Cell): Boolean {
        val iterator = index.iterator()
        val relation = iterator.locate(cell.id())

        // If "target" does not overlap any index cell, there is no intersection.
        if (relation == S2ShapeIndex.CellRelation.DISJOINT) return false

        // If "target" is subdivided into one or more index cells, there is an
        // intersection to within the S2ShapeIndex error bound (see Contains).
        if (relation == S2ShapeIndex.CellRelation.SUBDIVIDED) return true

        // If "target" is an index cell, there is an intersection because index cells
        // are created only if they have at least one edge or they are entirely
        // contained by the loop.
        if (iterator.id() == cell.id()) return true

        // Otherwise check if any edges intersect "target".
        if (boundaryApproxIntersects(iterator, cell)) return true

        // Otherwise check if the loop contains the center of "target".
        return contains(iterator, cell.getCenter())
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
                        unindexedContainsCalls.incrementAndGet() != kMaxUnindexedContainsCalls)) {
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
    fun containsNested(b: S2Loop): Boolean {
        if (!subregionBound.contains(b.bound)) return false

        // Special cases to handle either loop being empty or full.  Also bail out
        // when B has no vertices to avoid heap overflow on the vertex(1) call
        // below.  (This method is called during polygon initialization before the
        // client has an opportunity to call IsValid().)
        if (isEmptyOrFull() || b.numVertices() < 2) {
            return isFull() || b.isEmpty()
        }

        // We are given that A and B do not share any edges, and that either one
        // loop contains the other or they do not intersect.
        val m = findVertex(b.vertex(1))
        if (m < 0) {
            // Since b->vertex(1) is not shared, we can check whether A contains it.
            return contains(b.vertex(1))
        }
        // Check whether the edge order around b->vertex(1) is compatible with
        // A containing B.
        return S2WedgeRelations.wedgeContains(vertex(m - 1), vertex(m), vertex(m + 1), b.vertex(0), b.vertex(2))
    }

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
    fun compareBoundary(b: S2Loop): Int {
        assert { !isEmpty() && !b.isEmpty() }
        assert { !b.isFull() || !b.isHole() }

        // The bounds must intersect for containment or crossing.
        if (!bound.intersects(b.bound)) return -1

        // Full loops are handled as though the loop surrounded the entire sphere.
        if (isFull()) return 1
        if (b.isFull()) return -1

        // Check whether there are any edge crossings, and also check the loop
        // relationship at any shared vertices.
        val relation = CompareBoundaryRelation(b.isHole())
        if (hasCrossingRelation(this, b, relation)) return 0
        if (relation.found_shared_vertex) {
            return if (relation.contains_edge) 1 else -1
        }

        // There are no edge intersections or shared vertices, so we can check
        // whether A contains an arbitrary vertex of B.
        return if (contains(b.vertex(0))) 1 else -1
    }

    // Given two loops whose boundaries do not cross (see CompareBoundary),
    // return true if A contains the boundary of B.  If "reverse_b" is true, the
    // boundary of B is reversed first (which only affects the result when there
    // are shared edges).  This method is cheaper than CompareBoundary() because
    // it does not test for edge intersections.
    //
    // REQUIRES: neither loop is empty.
    // REQUIRES: if b->is_full(), then reverse_b == false.
    fun containsNonCrossingBoundary(b: S2Loop, reverse_b: Boolean): Boolean {
        assert { !isEmpty() && !b.isEmpty() }
        assert { !b.isFull() || !reverse_b }

        // The bounds must intersect for containment.
        if (!bound.intersects(b.bound)) return false

        // Full loops are handled as though the loop surrounded the entire sphere.
        if (isFull()) return true
        if (b.isFull()) return false

        val m = findVertex(b.vertex(0))
        if (m < 0) {
            // Since vertex b0 is not shared, we can check whether A contains it.
            return contains(b.vertex(0))
        }
        // Otherwise check whether the edge (b0, b1) is contained by A.
        return wedgeContainsSemiwedge(vertex(m - 1), vertex(m), vertex(m + 1), b.vertex(1), reverse_b)
    }

    // Wrapper class for indexing a loop (see S2ShapeIndex).  Once this object
    // is inserted into an S2ShapeIndex it is owned by that index, and will be
    // automatically deleted when no longer needed by the index.  Note that this
    // class does not take ownership of the loop itself (see OwningShape below).
    // You can also subtype this class to store additional data (see S2Shape for
    // details).
    class Shape(id: Int = 0, val loop: S2Loop) : S2Shape(id) {

        // S2Shape interface:

        override val numEdges: Int
            get() = if (loop.isEmptyOrFull()) 0 else loop.numVertices()

        override fun edge(edgeId: Int): Edge = Edge(loop.vertex(edgeId), loop.vertex(edgeId + 1))

        override val dimension: Int = 2

        override fun getReferencePoint(): ReferencePoint = ReferencePoint(S2Point.origin(), loop.containsOrigin())

        override val numChains: Int
            get() = if (loop.isEmpty()) 0 else 1

        override fun chain(chain_id: Int): Chain {
            assertEQ(chain_id, 0)
            return Chain(0, numEdges)
        }

        override fun chainEdge(chainId: Int, offset: Int): Edge {
            assertEQ(chainId, 0)
            return Edge(loop.vertex(offset), loop.vertex(offset + 1))
        }

        override fun chainPosition(edgeId: Int): ChainPosition = ChainPosition(0, edgeId)

        override val typeTag: TypeTag = kNoTypeTag


    }

    // List of vertices with operator[] maps index values in the range
    // [n, 2*n-1] to the range [0, n-1] by subtracting n (where n == size()).
    // In other words, two full copies of the vertex array are available.  (This
    // is a compromise between convenience and efficiency, since computing the
    // index modulo "n" is surprisingly expensive.)
    //
    // This property is useful for implementing algorithms where the elements of
    // the span represent the vertices of a loop.
    class S2PointLoopSpan(val points: List<S2Point> = emptyList()) : List<S2Point> by points {

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
            for (i in startIndex until (startIndex + points.size)) {
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
    private fun containsOrigin(): Boolean = originInside

    // A version of Contains(S2Point) that does not use the S2ShapeIndex.
    // Used by the S2Polygon implementation.
    private fun bruteForceContains(p: S2Point): Boolean {
        // Empty and full loops don't need a special case, but invalid loops with
        // zero vertices do, so we might as well handle them all at once.
        if (numVertices() < 3) return originInside

        val origin = S2Point.origin()
        val crosser = S2EdgeCrosser(origin, p, vertex(0))
        var inside = originInside
        for (i in 1..numVertices()) {
            inside = inside xor crosser.edgeOrVertexCrossing(vertex(i))
        }
        return inside;
    }

    // Like FindValidationError(), but skips any checks that would require
    // building the S2ShapeIndex (i.e., self-intersection tests).  This is used
    // by the S2Polygon implementation, which uses its own index to check for
    // loop self-intersections.
    private fun findValidationErrorNoIndex(): S2Error {
        // subregion_bound_ must be at least as large as bound_.  (This is an
        // internal consistency check rather than a test of client data.)
        assert { subregionBound.contains(bound) }

        // All vertices must be unit length.  (Unfortunately this check happens too
        // late in debug mode, because S2Loop construction calls s2pred::Sign which
        // expects vertices to be unit length.  But it is still a useful check in
        // optimized builds.)
        for (i in 0 until numVertices()) {
            if (!S2Point.isUnitLength(vertex(i))) {
                return S2Error(S2Error.NOT_UNIT_LENGTH, "Vertex %d is not unit length".format(i))
            }
        }
        // Loops must have at least 3 vertices (except for the empty and full loops).
        if (numVertices() < 3) {
            if (isEmptyOrFull()) {
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
            if (vertex(i) == vertex(i + 1)) {
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
    fun getXYZFaceSiTiVertices(): List<S2XYZFaceSiTi> {
        return vertices.asSequence().map { p ->
            val (cellLevel, faceSiTi) = S2Coords.xyzToFaceSiTi(p)
            S2XYZFaceSiTi(xyz = p, cellLevel = cellLevel, face = faceSiTi.face, si = faceSiTi.si, ti = faceSiTi.ti)
        }.toList()
    }

    // Given an iterator that is already positioned at the S2ShapeIndexCell
    // containing "p", returns Contains(p).
    fun contains(iterator: MutableS2ShapeIndex.Iterator, p: S2Point): Boolean {
        // Test containment by drawing a line segment from the cell center to the
        // given point and counting edge crossings.
        val a_clipped = iterator.cell()!!.clipped(0)
        var inside = a_clipped.containsCenter
        val a_num_edges = a_clipped.numEdges()
        if (a_num_edges > 0) {
            val center = iterator.center()
            val crosser = S2EdgeCrosser(center, p)
            var ai_prev = -2
            for (i in 0 until a_num_edges) {
                val ai = a_clipped.edge(i)
                if (ai != ai_prev + 1) crosser.restartAt(vertex(ai))
                ai_prev = ai
                inside = inside xor crosser.edgeOrVertexCrossing(vertex(ai + 1))
            }
        }
        return inside
    }

    // Returns true if the loop boundary intersects "target".  It may also
    // return true when the loop boundary does not intersect "target" but
    // some edge comes within the worst-case error tolerance.
    //
    // REQUIRES: it.id().contains(target.id())
    // [This condition is true whenever it.Locate(target) returns INDEXED.]
    fun boundaryApproxIntersects(iterator: MutableS2ShapeIndex.Iterator, target: S2Cell): Boolean {
        assert { iterator.id().contains(target.id()) }
        val a_clipped = iterator.cell()!!.clipped(0)
        val a_num_edges = a_clipped.numEdges()

        // If there are no edges, there is no intersection.
        if (a_num_edges == 0) return false

        // We can save some work if "target" is the index cell itself.
        if (iterator.id() == target.id()) return true

        // Otherwise check whether any of the edges intersect "target".
        val kMaxError = (S2EdgeClipping.kFaceClipErrorUVCoord + S2EdgeClipping.kIntersectsRectErrorUVDist)
        val bound = target.boundUV().expanded(kMaxError)
        for (i in 0 until a_num_edges) {
            val ai = a_clipped.edge(i)
            val clippedEdge = S2EdgeClipping.clipToPaddedFace(vertex(ai), vertex(ai + 1), target.face(), kMaxError)
            if (clippedEdge != null && S2EdgeClipping.intersectsRect(clippedEdge.first, clippedEdge.second, bound)) {
                return true
            }
        }
        return false
    }

    // Returns an index "first" and a direction "dir" such that the vertex
    // sequence (first, first + dir, ..., first + (n - 1) * dir) does not change
    // when the loop vertex order is rotated or reversed.  This allows the loop
    // vertices to be traversed in a canonical order.
    fun getCanonicalLoopOrder(): S2LoopMeasures.LoopOrder = S2LoopMeasures.getCanonicalLoopOrder(verticesSpan())

    // Returns the index of a vertex at point "p", or -1 if not found.
    // The return value is in the range 1..num_vertices_ if found.
    private fun findVertex(p: S2Point): Int {
        if (numVertices() < 10) {
            // Exhaustive search.  Return value must be in the range [1..N].
            for (i in 1..numVertices()) {
                if (vertex(i) == p) return i
            }
            return -1
        }
        val iter = index.iterator()
        if (!iter.locate(p)) return -1

        val aClipped = iter.cell()!!.clipped(0)
        var i = aClipped.numEdges() - 1
        while (i >= 0) {
            val ai = aClipped.edge(i)
            // Return value must be in the range [1..N].
            if (vertex(ai) == p) return if (ai == 0) numVertices() else ai
            if (vertex(ai + 1) == p) return ai + 1
            --i
        }
        return -1;
    }

    // When the loop is modified (Invert(), or Init() called again) then the
    // indexing structures need to be cleared since they become invalid.
    private fun clearIndex(): Unit {
        unindexedContainsCalls.set(0)
        index.clear()
    }

    override fun toString(): String {
        return "S2Loop(vertices=$vertices)"
    }


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
        fun makeRegularLoop(center: S2Point, radius: S1Angle, num_vertices: Int): S2Loop {
            return makeRegularLoop(S2Point.getFrame(center), radius, num_vertices)
        }

        // Like the function above, but this version constructs a loop centered
        // around the z-axis of the given coordinate frame, with the first vertex in
        // the direction of the positive x-axis.  (This allows the loop to be
        // rotated for testing purposes.)
        fun makeRegularLoop(frame: Matrix3x3, radius: S1Angle, num_vertices: Int): S2Loop {
            // We construct the loop in the given frame coordinates, with the center at
            // (0, 0, 1).  For a loop of radius "r", the loop vertices have the form
            // (x, y, z) where x^2 + y^2 = sin(r) and z = cos(r).  The distance on the
            // sphere (arc length) from each vertex to the center is acos(cos(r)) = r.
            val z = cos(radius.radians)
            val r = sin(radius.radians)
            val radian_step = 2 * M_PI / num_vertices
            val vertices = mutableListOf<S2Point>()
            for (i in 0 until num_vertices) {
                val angle = i * radian_step
                val p = S2Point(r * cos(angle), r * sin(angle), z)
                vertices.add(S2Point.fromFrame(frame, p).normalize())
            }
            return S2Loop(vertices)
        }

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
        fun hasCrossingRelation(a: S2Loop, b: S2Loop, relation: LoopRelation): Boolean {
            // We look for S2CellId ranges where the indexes of A and B overlap, and
            // then test those edges for crossings.
            val ai = RangeIterator(a.index)
            val bi = RangeIterator(b.index)
            val ab = LoopCrosser(a, b, relation, false)  // Tests edges of A against B
            val ba = LoopCrosser(b, a, relation, true)   // Tests edges of B against A
            while (!ai.done() || !bi.done()) {
                if (ai.rangeMax() < bi.rangeMin()) {
                    // The A and B cells don't overlap, and A precedes B.
                    ai.seekTo(bi);
                } else if (bi.rangeMax() < ai.rangeMin()) {
                    // The A and B cells don't overlap, and B precedes A.
                    bi.seekTo(ai);
                } else {
                    // One cell contains the other.  Determine which cell is larger.
                    val ab_relation = ai.id().lsb() - bi.id().lsb()
                    if (ab_relation > 0UL) {
                        // A's index cell is larger.
                        if (ab.hasCrossingRelation(ai, bi)) return true
                    } else if (ab_relation < 0UL) {
                        // B's index cell is larger.
                        if (ba.hasCrossingRelation(bi, ai)) return true
                    } else {
                        // The A and B cells are the same.  Since the two cells have the same
                        // center point P, check whether P satisfies the crossing targets.
                        if (ai.containsCenter(0) && ab.a_crossing_target == 1 &&
                                bi.containsCenter(0) && ab.b_crossing_target == 1) {
                            return true
                        }
                        // Otherwise test all the edge crossings directly.
                        if (ai.numEdges(0) > 0 && bi.numEdges(0) > 0 &&
                                ab.cellCrossesCell(ai.clipped(0), bi.clipped(0))) {
                            return true
                        }
                        ai.next()
                        bi.next()
                    }
                }
            }
            return false;
        }

        fun matchBoundaries(a: S2Loop, b: S2Loop, a_offset: Int, max_error: S1Angle): Boolean {
            // The state consists of a pair (i,j).  A state transition consists of
            // incrementing either "i" or "j".  "i" can be incremented only if
            // a(i+1+a_offset) is near the edge from b(j) to b(j+1), and a similar rule
            // applies to "j".  The function returns true iff we can proceed all the way
            // around both loops in this way.
            //
            // Note that when "i" and "j" can both be incremented, sometimes only one
            // choice leads to a solution.  We handle this using a stack and
            // backtracking.  We also keep track of which states have already been
            // explored to avoid duplicating work.

            val pending = mutableListOf<Pair<Int, Int>>()
            val done = mutableSetOf<Pair<Int, Int>>()
            pending.add(Pair(0, 0))
            while (pending.isNotEmpty()) {
                val i = pending.last().first
                val j = pending.last().second
                pending.removeLast()
                if (i == a.numVertices() && j == b.numVertices()) {
                    return true
                }
                done.add(Pair(i, j))

                // If (i == na && offset == na-1) where na == a->num_vertices(), then
                // then (i+1+offset) overflows the [0, 2*na-1] range allowed by vertex().
                // So we reduce the range if necessary.
                var io = i + a_offset
                if (io >= a.numVertices()) io -= a.numVertices()

                if (i < a.numVertices() &&
                        !done.contains(Pair(i + 1, j)) &&
                        S2EdgeDistances.getDistance(a.vertex(io + 1), b.vertex(j), b.vertex(j + 1)) <= max_error
                ) {
                    pending.add(Pair(i + 1, j))
                }
                if (j < b.numVertices() &&
                        !done.contains(Pair(i, j + 1)) &&
                        S2EdgeDistances.getDistance(b.vertex(j + 1), a.vertex(io), a.vertex(io + 1)) <= max_error) {
                    pending.add(Pair(i, j + 1))
                }
            }
            return false
        }

        // Returns true if the wedge (a0, ab1, a2) contains the "semiwedge" defined as
        // any non-empty open set of rays immediately CCW from the edge (ab1, b2).  If
        // "reverse_b" is true, then substitute "clockwise" for "CCW"; this simulates
        // what would happen if the direction of loop B was reversed.
        fun wedgeContainsSemiwedge(a0: S2Point, ab1: S2Point, a2: S2Point, b2: S2Point, reverse_b: Boolean): Boolean {
            if (b2 == a0 || b2 == a2) {
                // We have a shared or reversed edge.
                return (b2 == a0) == reverse_b;
            } else {
                return S2Predicates.orderedCCW(a0, a2, b2, ab1)
            }
        }
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
    abstract fun wedgesCross(a0: S2Point, ab1: S2Point, a2: S2Point, b0: S2Point, b2: S2Point): Boolean

}

// Loop relation for Contains().
class ContainsRelation : LoopRelation() {
    var found_shared_vertex: Boolean = false

    // If A.Contains(P) == false && B.Contains(P) == true, it is equivalent to
    // having an edge crossing (i.e., Contains returns false).
    override fun a_crossing_target(): Int = 0
    override fun b_crossing_target(): Int = 1

    override fun wedgesCross(a0: S2Point, ab1: S2Point, a2: S2Point, b0: S2Point, b2: S2Point): Boolean {
        found_shared_vertex = true
        return !S2WedgeRelations.wedgeContains(a0, ab1, a2, b0, b2)
    }

}

// Loop relation for Intersects().
class IntersectsRelation : LoopRelation() {
    var found_shared_vertex: Boolean = false

    // If A.Contains(P) == true && B.Contains(P) == true, it is equivalent to
    // having an edge crossing (i.e., Intersects returns true).
    override fun a_crossing_target(): Int = 1
    override fun b_crossing_target(): Int = 1

    override fun wedgesCross(a0: S2Point, ab1: S2Point, a2: S2Point, b0: S2Point, b2: S2Point): Boolean {
        found_shared_vertex = true;
        return S2WedgeRelations.wedgeIntersects(a0, ab1, a2, b0, b2)
    }

};


// Loop relation for CompareBoundary().
class CompareBoundaryRelation(reversed_b: Boolean) : LoopRelation() {

    var reverse_b: Boolean = reversed_b       // True if loop B should be reversed.
    var found_shared_vertex: Boolean = false  // True if any wedge was processed.
    var contains_edge: Boolean = false        // True if any edge of B is contained by A.
    var excludes_edge: Boolean = false        // True if any edge of B is excluded by A.

    // The CompareBoundary relation does not have a useful early-exit condition,
    // so we return -1 for both crossing targets.
    //
    // Aside: A possible early exit condition could be based on the following.
    //   If A contains a point of both B and ~B, then A intersects Boundary(B).
    //   If ~A contains a point of both B and ~B, then ~A intersects Boundary(B).
    //   So if the intersections of {A, ~A} with {B, ~B} are all non-empty,
    //   the return value is 0, i.e., Boundary(A) intersects Boundary(B).
    // Unfortunately it isn't worth detecting this situation because by the
    // time we have seen a point in all four intersection regions, we are also
    // guaranteed to have seen at least one pair of crossing edges.
    override fun a_crossing_target(): Int = -1
    override fun b_crossing_target(): Int = -1

    override fun wedgesCross(a0: S2Point, ab1: S2Point, a2: S2Point, b0: S2Point, b2: S2Point): Boolean {
        // Because we don't care about the interior of B, only its boundary, it is
        // sufficient to check whether A contains the semiwedge (ab1, b2).
        found_shared_vertex = true
        if (S2Loop.wedgeContainsSemiwedge(a0, ab1, a2, b2, reverse_b)) {
            contains_edge = true
        } else {
            excludes_edge = true
        }
        return contains_edge and excludes_edge
    }

}

// LoopCrosser is a helper class for determining whether two loops cross.
// It is instantiated twice for each pair of loops to be tested, once for the
// pair (A,B) and once for the pair (B,A), in order to be able to process
// edges in either loop nesting order.
class LoopCrosser(val a: S2Loop, val b: S2Loop, val relation: LoopRelation, val swapped: Boolean) {

    val a_crossing_target: Int = if (swapped) relation.b_crossing_target() else relation.a_crossing_target()
    val b_crossing_target: Int = if (swapped) relation.a_crossing_target() else relation.b_crossing_target()

    // State maintained by StartEdge() and EdgeCrossesCell().
    private val crosser: S2EdgeCrosser = S2EdgeCrosser()
    private var aj: Int = -1
    private var bjPrev: Int = -1

    // Temporary data declared here to avoid repeated memory allocations.
    private val bQuery: S2CrossingEdgeQuery = S2CrossingEdgeQuery(b.index)
    private val bCells = mutableListOf<S2ShapeIndexCell>()

    // Given two iterators positioned such that ai->id().Contains(bi->id()),
    // return true if there is a crossing relationship anywhere within ai->id().
    // Specifically, this method returns true if there is an edge crossing, a
    // wedge crossing, or a point P that matches both "crossing targets".
    // Advances both iterators past ai->id().
    fun hasCrossingRelation(ai: RangeIterator, bi: RangeIterator): Boolean {
        assert { ai.id().contains(bi.id()) }
        if (ai.numEdges(0) == 0) {
            if (ai.containsCenter(0) && a_crossing_target == 1) {
            // All points within ai->id() satisfy the crossing target for A, so it's
            // worth iterating through the cells of B to see whether any cell
            // centers also satisfy the crossing target for B.
            do {
                if (bi.containsCenter(0) && b_crossing_target == 1) return true
                bi.next()
            } while (bi.id() <= ai.rangeMax())
        } else {
            // The crossing target for A is not satisfied, so we skip over the cells
            // of B using binary search.
            bi.seekBeyond(ai)
        }
        } else {
            // The current cell of A has at least one edge, so check for crossings.
            if (hasCrossing(ai, bi)) return true
        }
        ai.next()
        return false
    }

    // Given two index cells, return true if there are any edge crossings or
    // wedge crossings within those cells.
    fun cellCrossesCell(a_clipped: S2ClippedShape, b_clipped: S2ClippedShape): Boolean {
        // Test all edges of "a_clipped" against all edges of "b_clipped".
        val aNumEdges = a_clipped.numEdges()
        for (i in 0 until aNumEdges) {
            startEdge(a_clipped.edge(i))
            if (edgeCrossesCell(b_clipped)) return true
        }
        return false
    }

    // Given two iterators positioned such that ai->id().Contains(bi->id()),
    // return true if there is an edge crossing or wedge crosssing anywhere
    // within ai->id().  Advances "bi" (only) past ai->id().
    private fun hasCrossing(ai: RangeIterator, bi: RangeIterator): Boolean {
        assert { ai.id().contains(bi.id()) }
        // If ai->id() intersects many edges of B, then it is faster to use
        // S2CrossingEdgeQuery to narrow down the candidates.  But if it intersects
        // only a few edges, it is faster to check all the crossings directly.
        // We handle this by advancing "bi" and keeping track of how many edges we
        // would need to test.

        val kEdgeQueryMinEdges = 20;  // Tuned using benchmarks.
        var totalEdges = 0
        bCells.clear()
        do {
            if (bi.numEdges(0) > 0) {
                totalEdges += bi.cell()!!.num_edges()
                if (totalEdges >= kEdgeQueryMinEdges) {
                    // There are too many edges to test them directly, so use
                    // S2CrossingEdgeQuery.
                    if (cellCrossesAnySubcell(ai.clipped(0), ai.id())) return true
                    bi.seekBeyond(ai)
                    return false;
                }
                bCells.add(bi.cell()!!)
            }
            bi.next()
        } while (bi.id() <= ai.rangeMax())

        // Test all the edge crossings directly.
        for (b_cell in bCells) {
            if (cellCrossesCell(ai.clipped(0), b_cell.clipped(0))) {
            return true
        }
        }
        return false
    }

    // Given an index cell of A, return true if there are any edge or wedge
    // crossings with any index cell of B contained within "b_id".
    private fun cellCrossesAnySubcell(a_clipped: S2ClippedShape, b_id: S2CellId): Boolean {
        // Test all edges of "a_clipped" against all edges of B.  The relevant B
        // edges are guaranteed to be children of "b_id", which lets us find the
        // correct index cells more efficiently.
        val bRoot = S2PaddedCell(b_id, 0.0)
        val aNumEdges = a_clipped.numEdges()
        for (i in 0 until aNumEdges) {
            val aj = a_clipped.edge(i)
            // Use an S2CrossingEdgeQuery starting at "b_root" to find the index cells
            // of B that might contain crossing edges.
            bQuery.getCells(a.vertex(aj), a.vertex(aj + 1), bRoot, bCells)
            if (bCells.isEmpty()) continue
            startEdge(aj)
            for (b_cell in bCells) {
                if (edgeCrossesCell(b_cell.clipped(0))) return true
            }
        }
        return false;
    }

    // Prepare to check the given edge of loop A for crossings.
    private fun startEdge(aj: Int) {
        // Start testing the given edge of A for crossings.
        crosser.init(a.vertex(aj), a.vertex(aj + 1));
        this.aj = aj;
        bjPrev = -2;
    }

    // Check the current edge of loop A for crossings with all edges of the
    // given index cell of loop B.
    private fun edgeCrossesCell(b_clipped: S2ClippedShape): Boolean {
        // Test the current edge of A against all edges of "b_clipped".
        val bNumEdges = b_clipped.numEdges()
        for (j in 0 until bNumEdges) {
            val bj = b_clipped.edge(j)
            if (bj != bjPrev + 1) crosser.restartAt(b.vertex(bj))
            bjPrev = bj
            val crossing = crosser.crossingSign(b.vertex(bj + 1))
            if (crossing < 0) continue
            if (crossing > 0) return true
            // We only need to check each shared vertex once, so we only
            // consider the case where a_vertex(aj_+1) == b_.vertex(bj+1).
            if (a.vertex(aj + 1) == b.vertex(bj + 1)) {
                if (swapped) {
                    if (relation.wedgesCross(b.vertex(bj), b.vertex(bj + 1), b.vertex(bj + 2), a.vertex(aj), a.vertex(aj + 2))) {
                        return true
                    }
                } else if (relation.wedgesCross(a.vertex(aj), a.vertex(aj + 1), a.vertex(aj + 2), b.vertex(bj), b.vertex(bj + 2))) {
                    return true
                }
            }
        }
        return false;
    }

};
