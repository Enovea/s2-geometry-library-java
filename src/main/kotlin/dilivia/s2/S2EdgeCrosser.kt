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
package dilivia.s2

import com.google.common.geometry.S2
import com.google.common.geometry.S2.DBL_EPSILON

/**
 * This class allows edges to be efficiently tested for intersection with a given fixed edge AB.  It is especially
 * efficient when testing for intersection with an edge chain connecting vertices v0, v1, v2, ...
 *
 * Example usage:
 *
 *   fun countIntersections(a: S2Point, b: S2Point, edges: List<Pair<S2Point, S2Point>>) {
 *     val count: Int = 0
 *     val crosser = S2EdgeCrosser(a, b)
 *     for (edge : edges) {
 *       if (crosser.crossingSign(edge.first, edge.second) >= 0) {
 *         ++count;
 *       }
 *     }
 *     return count;
 *   }
 *
 * This class expects that the client already has all the necessary vertices stored in memory, so that this class can
 * refer to them with pointers and does not need to make its own copies.  If this is not the case (e.g., you
 * want to pass temporary objects as vertices), see S2CopyingEdgeCrosser.
 *
 * @constructor Default constructor; must be followed by a call to init().
 */
class S2EdgeCrosser() {

    // The fields below are constant after the call to Init().
    lateinit var a: S2Point
    lateinit var b: S2Point
    private lateinit var a_cross_b: S2Point

    // To reduce the number of calls to s2pred::ExpensiveSign(), we compute an
    // outward-facing tangent at A and B if necessary.  If the plane
    // perpendicular to one of these tangents separates AB from CD (i.e., one
    // edge on each side) then there is no intersection.
    private var have_tangents: Boolean = false        // True if the tangents have been computed.
    private lateinit var a_tangent: S2Point           // Outward-facing tangent at A.
    private lateinit var b_tangent: S2Point           // Outward-facing tangent at B.

    // The fields below are updated for each vertex in the chain.
    var c: S2Point? = null                            // Previous vertex in the vertex chain.
    private var acb: Int = 0                          // The orientation of triangle ACB.

    // The field below is a temporary used by CrossingSignInternal().
    private var bda: Int = 0                          // The orientation of triangle BDA.

    // Convenience constructor that calls Init() with the given fixed edge AB.
    // The arguments "a" and "b" must point to values that persist for the
    // lifetime of the S2EdgeCrosser object (or until the next Init() call).
    constructor(a: S2Point, b: S2Point): this() {
        init(a, b)
    }

    // Convenience constructor that uses AB as the fixed edge, and C as the
    // first vertex of the vertex chain (equivalent to calling RestartAt(c)).
    //
    // The arguments must point to values that persist until the next call.
    constructor(a: S2Point, b: S2Point, c: S2Point): this() {
        init(a, b)
        restartAt(c)
    }

    constructor(crosser: S2EdgeCrosser): this() {

    }

    // Initialize the S2EdgeCrosser with the given fixed edge AB.  The arguments
    // "a" and "b" must point to values that persist for the lifetime of the
    // S2EdgeCrosser object (or until the next init() call).
    fun init(a: S2Point, b: S2Point) {
        Assertions.assertPointIsUnitLength(a)
        Assertions.assertPointIsUnitLength(b)
        this.a = a
        this.b = b
        this.a_cross_b = a.crossProd(b)
        have_tangents = false
        c = null
    }

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
    // Note that if you want to check an edge against a chain of other edges,
    // it is slightly more efficient to use the single-argument version of
    // CrossingSign below.
    //
    // The arguments must point to values that persist until the next call.
    fun crossingSign(c: S2Point, d: S2Point): Int {
        if (this.c != c) restartAt(c)
        return crossingSign(d)
    }

    // This method extends the concept of a "crossing" to the case where AB
    // and CD have a vertex in common.  The two edges may or may not cross,
    // according to the rules defined in VertexCrossing() below.  The rules
    // are designed so that point containment tests can be implemented simply
    // by counting edge crossings.  Similarly, determining whether one edge
    // chain crosses another edge chain can be implemented by counting.
    //
    // Returns true if CrossingSign(c, d) > 0, or AB and CD share a vertex
    // and VertexCrossing(a, b, c, d) returns true.
    //
    // The arguments must point to values that persist until the next call.
    fun edgeOrVertexCrossing(c: S2Point, d: S2Point): Boolean {
        if (this.c != c) restartAt(c)
        return edgeOrVertexCrossing(d)
    }

    ///////////////////////// Edge Chain Methods ///////////////////////////
    //
    // You don't need to use these unless you're trying to squeeze out every
    // last drop of performance.  Essentially all you are saving is a test
    // whether the first vertex of the current edge is the same as the second
    // vertex of the previous edge.  Example usage:
    //
    //   vector<S2Point> chain;
    //   crosser.RestartAt(&chain[0]);
    //   for (int i = 1; i < chain.size(); ++i) {
    //     if (crosser.EdgeOrVertexCrossing(&chain[i])) { ++count; }
    //   }

    // Call this method when your chain 'jumps' to a new place.
    // The argument must point to a value that persists until the next call.
    fun restartAt(c: S2Point) {
        Assertions.assertPointIsUnitLength(c)
        this.c = c
        this.acb = - S2Predicates.triageSign(a, b, c, a_cross_b)
    }

    // Like CrossingSign above, but uses the last vertex passed to one of
    // the crossing methods (or RestartAt) as the first vertex of the current
    // edge.
    //
    // The argument must point to a value that persists until the next call.
    fun crossingSign(d: S2Point): Int {
        Assertions.assertPointIsUnitLength(d)
        // For there to be an edge crossing, the triangles ACB, CBD, BDA, DAC must
        // all be oriented the same way (CW or CCW).  We keep the orientation of ACB
        // as part of our state.  When each new point D arrives, we compute the
        // orientation of BDA and check whether it matches ACB.  This checks whether
        // the points C and D are on opposite sides of the great circle through AB.

        // Recall that TriageSign is invariant with respect to rotating its
        // arguments, i.e. ABD has the same orientation as BDA.
        val bda = S2Predicates.triageSign(a, b, d, a_cross_b)
        if (acb == -bda && bda != 0) {
            // The most common case -- triangles have opposite orientations.  Save the
            // current vertex D as the next vertex C, and also save the orientation of
            // the new triangle ACB (which is opposite to the current triangle BDA).
            this.c = d
            this.acb = -bda
            return -1;
        }
        this.bda = bda
        return crossingSignInternal(d)
    }

    // Like EdgeOrVertexCrossing above, but uses the last vertex passed to one
    // of the crossing methods (or RestartAt) as the first vertex of the
    // current edge.
    //
    // The argument must point to a value that persists until the next call.
    fun edgeOrVertexCrossing(d: S2Point): Boolean {
        // We need to copy c_ since it is clobbered by CrossingSign().
        val c = this.c!!
        val crossing = crossingSign(d)
        if (crossing < 0) return false;
        if (crossing > 0) return true;
        return S2EdgeCrossings.vertexCrossing(a, b, c, d)
    }

    // These functions handle the "slow path" of CrossingSign().
    private fun crossingSignInternal(d: S2Point): Int {
        // Compute the actual result, and then save the current vertex D as the next
        // vertex C, and save the orientation of the next triangle ACB (which is
        // opposite to the current triangle BDA).
        val result = crossingSignInternal2(d)
        this.c = d
        this.acb = -this.bda
        return result;
    }

    private fun crossingSignInternal2(d: S2Point): Int {
        val c = this.c!!

        // At this point, a very common situation is that A,B,C,D are four points on
        // a line such that AB does not overlap CD.  (For example, this happens when
        // a line or curve is sampled finely, or when geometry is constructed by
        // computing the union of S2CellIds.)  Most of the time, we can determine
        // that AB and CD do not intersect by computing the two outward-facing
        // tangents at A and B (parallel to AB) and testing whether AB and CD are on
        // opposite sides of the plane perpendicular to one of these tangents.  This
        // is moderately expensive but still much cheaper than s2pred::ExpensiveSign.
        if (!have_tangents) {
            val norm = S2Point.robustCrossProd(a, b).normalize()
            a_tangent = a.crossProd(norm)
            b_tangent = norm.crossProd(b)
            have_tangents = true;
        }
        // The error in RobustCrossProd() is insignificant.  The maximum error in
        // the call to CrossProd() (i.e., the maximum norm of the error vector) is
        // (0.5 + 1/sqrt(3)) * DBL_EPSILON.  The maximum error in each call to
        // DotProd() below is DBL_EPSILON.  (There is also a small relative error
        // term that is insignificant because we are comparing the result against a
        // constant that is very close to zero.)
        val kError = (1.5 + 1/ S2.M_SQRT3) * DBL_EPSILON
        if ((c.dotProd(a_tangent) > kError && d.dotProd(a_tangent) > kError)
                || (c.dotProd(b_tangent) > kError && d.dotProd(b_tangent) > kError)) {
            return -1;
        }

        // Otherwise, eliminate the cases where two vertices from different edges
        // are equal.  (These cases could be handled in the code below, but we would
        // rather avoid calling ExpensiveSign whenever possible.)
        if (a == c || a == d || b == c || b == d) return 0;

        // Eliminate cases where an input edge is degenerate.  (Note that in most
        // cases, if CD is degenerate then this method is not even called because
        // acb_ and bda have different signs.)
        if (a == b || c == d) return -1;

        // Otherwise it's time to break out the big guns.
        if (acb == 0) acb = -S2Predicates.expensiveSign(a, b, c)
        assert(acb != 0)
        if (bda == 0) bda = S2Predicates.expensiveSign(a, b, d)
        assert(bda != 0)
        if (bda != acb) return -1;

        val c_cross_d = c.crossProd(d)
        val cbd = -S2Predicates.sign(c, d, b, c_cross_d)
        assert(cbd != 0)
        if (cbd != acb) return -1;
        val dac = S2Predicates.sign(c, d, a, c_cross_d)
        assert(dac != 0);
        return if(dac != acb) -1 else 1
    }


}