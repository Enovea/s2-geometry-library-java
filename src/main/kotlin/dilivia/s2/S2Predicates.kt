package dilivia.s2

import com.google.common.geometry.S2.DBL_EPSILON
import com.google.common.geometry.S2.M_SQRT2
import dilivia.s2.math.R3VectorExactFloat
import java.math.MathContext
import kotlin.math.abs
import kotlin.math.pow
import kotlin.math.sqrt

// This class contains various predicates that are guaranteed to produce
// correct, consistent results.  They are also relatively efficient.  This is
// achieved by computing conservative error bounds and falling back to high
// precision or even exact arithmetic when the result is uncertain.  Such
// predicates are useful in implementing robust algorithms.
//
// See also S2EdgeCrosser, which implements various exact
// edge-crossing predicates more efficiently than can be done here.

// S2EdgeUtil contains the following exact predicates that test for edge
// crossings.  (Usually you will want to use S2EdgeCrosser, which
// implements them much more efficiently.)
//
// int CrossingSign(const S2Point& a0, const S2Point& a1,
//                  const S2Point& b0, const S2Point& b1);
//
// bool EdgeOrVertexCrossing(const S2Point& a0, const S2Point& a1,
//                           const S2Point& b0, const S2Point& b1);
@ExperimentalUnsignedTypes
object S2Predicates {

    // Returns 2 ** (-digits).  This could be implemented using "ldexp" except
    // that std::ldexp is not constexpr in C++11.
    private fun epsilonForDigits(digits: Int): Double {
        return if (digits < 64) 2.0.pow(-digits) else epsilonForDigits(digits - 63) / (1.toULong() shl 63).toDouble()
    }

    val kDoubleRoundingEpsilon = epsilonForDigits(53)

    val kDebug: Boolean = true

    // A predefined S1ChordAngle representing (approximately) 45 degrees.
    val k45Degrees = S1ChordAngle.fromLength2(2 - M_SQRT2);

    // Returns +1 if the points A, B, C are counterclockwise, -1 if the points
    // are clockwise, and 0 if any two points are the same.  This function is
    // essentially like taking the sign of the determinant of ABC, except that
    // it has additional logic to make sure that the above properties hold even
    // when the three points are coplanar, and to deal with the limitations of
    // floating-point arithmetic.
    //
    // Sign satisfies the following conditions:
    //
    //  (1) Sign(a,b,c) == 0 if and only if a == b, b == c, or c == a
    //  (2) Sign(b,c,a) == Sign(a,b,c) for all a,b,c
    //  (3) Sign(c,b,a) == -Sign(a,b,c) for all a,b,c
    //
    // In other words:
    //
    //  (1) The result is zero if and only if two points are the same.
    //  (2) Rotating the order of the arguments does not affect the result.
    //  (3) Exchanging any two arguments inverts the result.
    //
    // On the other hand, note that it is not true in general that
    // Sign(-a,b,c) == -Sign(a,b,c), or any similar identities
    // involving antipodal points.
    fun sign(a: S2Point, b: S2Point, c: S2Point): Int {
        // We don't need RobustCrossProd() here because Sign() does its own
        // error estimation and calls ExpensiveSign() if there is any uncertainty
        // about the result.
        return sign(a, b, c, a.crossProd(b))
    }

    // Given 4 points on the unit sphere, return true if the edges OA, OB, and
    // OC are encountered in that order while sweeping CCW around the point O.
    // You can think of this as testing whether A <= B <= C with respect to the
    // CCW ordering around O that starts at A, or equivalently, whether B is
    // contained in the range of angles (inclusive) that starts at A and extends
    // CCW to C.  Properties:
    //
    //  (1) If OrderedCCW(a,b,c,o) && OrderedCCW(b,a,c,o), then a == b
    //  (2) If OrderedCCW(a,b,c,o) && OrderedCCW(a,c,b,o), then b == c
    //  (3) If OrderedCCW(a,b,c,o) && OrderedCCW(c,b,a,o), then a == b == c
    //  (4) If a == b or b == c, then OrderedCCW(a,b,c,o) is true
    //  (5) Otherwise if a == c, then OrderedCCW(a,b,c,o) is false
    fun orderedCCW(a: S2Point, b: S2Point, c: S2Point, o: S2Point) {
        TODO("Not yet implemented")
    }

    // Returns -1, 0, or +1 according to whether AX < BX, A == B, or AX > BX
    // respectively.  Distances are measured with respect to the positions of X,
    // A, and B as though they were reprojected to lie exactly on the surface of
    // the unit sphere.  Furthermore, this method uses symbolic perturbations to
    // ensure that the result is non-zero whenever A != B, even when AX == BX
    // exactly, or even when A and B project to the same point on the sphere.
    // Such results are guaranteed to be self-consistent, i.e. if AB < BC and
    // BC < AC, then AB < AC.
    fun compareDistances(x: S2Point, a: S2Point, b: S2Point): Int {
        TODO("Not yet implemented")
    }

    // Returns -1, 0, or +1 according to whether the distance XY is less than,
    // equal to, or greater than "r" respectively.  Distances are measured with
    // respect the positions of all points as though they are projected to lie
    // exactly on the surface of the unit sphere.
    fun compareDistance(x: S2Point, y: S2Point, r: S1ChordAngle): Int {
        TODO("Not yet implemented")
    }

    // Returns -1, 0, or +1 according to whether the distance from the point X to
    // the edge A is less than, equal to, or greater than "r" respectively.
    // Distances are measured with respect the positions of all points as though
    // they were projected to lie exactly on the surface of the unit sphere.
    //
    // REQUIRES: A0 and A1 do not project to antipodal points (e.g., A0 == -A1).
    //           This requires that (A0 != C * A1) for any constant C < 0.
    //
    // NOTE(ericv): All of the predicates defined here could be extended to handle
    // edges consisting of antipodal points by implementing additional symbolic
    // perturbation logic (similar to Sign) in order to rigorously define the
    // direction of such edges.
    fun compareEdgeDistance(x: S2Point, a0: S2Point, a1: S2Point, r: S1ChordAngle): Int {
        TODO("Not yet implemented")
    }

    // Returns -1, 0, or +1 according to whether the normal of edge A has
    // negative, zero, or positive dot product with the normal of edge B.  This
    // essentially measures whether the edges A and B are closer to proceeding in
    // the same direction or in opposite directions around the sphere.
    //
    // This method returns an exact result, i.e. the result is zero if and only if
    // the two edges are exactly perpendicular or at least one edge is degenerate.
    // (i.e., both edge endpoints project to the same point on the sphere).
    //
    // CAVEAT: This method does not use symbolic perturbations.  Therefore it can
    // return zero even when A0 != A1 and B0 != B1, e.g. if (A0 == C * A1) exactly
    // for some constant C > 0 (which is possible even when both points are
    // considered "normalized").
    //
    // REQUIRES: Neither edge can consist of antipodal points (e.g., A0 == -A1)
    //           (see comments in CompareEdgeDistance).
    fun compareEdgeDirections(a0: S2Point, a1: S2Point, b0: S2Point, b1: S2Point): Int {
        TODO("Not yet implemented")
    }

    // Returns Sign(X0, X1, Z) where Z is the circumcenter of triangle ABC.
    // The return value is +1 if Z is to the left of edge X, and -1 if Z is to the
    // right of edge X.  The return value is zero if A == B, B == C, or C == A
    // (exactly), and also if X0 and X1 project to identical points on the sphere
    // (e.g., X0 == X1).
    //
    // The result is determined with respect to the positions of all points as
    // though they were projected to lie exactly on the surface of the unit
    // sphere.  Furthermore this method uses symbolic perturbations to compute a
    // consistent non-zero result even when Z lies exactly on edge X.
    //
    // REQUIRES: X0 and X1 do not project to antipodal points (e.g., X0 == -X1)
    //           (see comments in CompareEdgeDistance).
    fun edgeCircumcenterSign(x0: S2Point, x1: S2Point, a: S2Point, b: S2Point, c: S2Point): Int {
        TODO("Not yet implemented")
    }

    // This is a specialized method that is used to compute the intersection of an
    // edge X with the Voronoi diagram of a set of points, where each Voronoi
    // region is intersected with a disc of fixed radius "r".
    //
    // Given two sites A and B and an edge (X0, X1) such that d(A,X0) < d(B,X0)
    // and both sites are within the given distance "r" of edge X, this method
    // intersects the Voronoi region of each site with a disc of radius r and
    // determines whether either region has an empty intersection with edge X.  It
    // returns FIRST if site A has an empty intersection, SECOND if site B has an
    // empty intersection, NEITHER if neither site has an empty intersection, or
    // UNCERTAIN if A == B exactly.  Note that it is not possible for both
    // intersections to be empty because of the requirement that both sites are
    // within distance r of edge X.  (For example, the only reason that Voronoi
    // region A can have an empty intersection with X is that site B is closer to
    // all points on X that are within radius r of site A.)
    //
    // The result is determined with respect to the positions of all points as
    // though they were projected to lie exactly on the surface of the unit
    // sphere.  Furthermore this method uses symbolic perturbations to compute a
    // consistent non-zero result even when A and B lie on opposite sides of X
    // such that the Voronoi edge between them exactly coincides with edge X, or
    // when A and B are distinct but project to the same point on the sphere
    // (i.e., they are linearly dependent).
    //
    // REQUIRES: r < S1ChordAngle::Right() (90 degrees)
    // REQUIRES: s2pred::CompareDistances(x0, a, b) < 0
    // REQUIRES: s2pred::CompareEdgeDistance(a, x0, x1, r) <= 0
    // REQUIRES: s2pred::CompareEdgeDistance(b, x0, x1, r) <= 0
    // REQUIRES: X0 and X1 do not project to antipodal points (e.g., X0 == -X1)
    //           (see comments in CompareEdgeDistance).
    enum class Excluded { FIRST, SECOND, NEITHER, UNCERTAIN }

    fun getVoronoiSiteExclusion(a: S2Point, b: S2Point, x0: S2Point, x1: S2Point, r: S1ChordAngle): Excluded {
        TODO("Not yet implemented")
    }

    /////////////////////////// Low-Level Methods ////////////////////////////
    //
    // Most clients will not need the following methods.  They can be slightly
    // more efficient but are harder to use, since they require the client to do
    // all the actual crossing tests.

    // A more efficient version of Sign that allows the precomputed
    // cross-product of A and B to be specified.  (Unlike the 3 argument
    // version this method is also inlined.)
    fun sign(a: S2Point, b: S2Point, c: S2Point, aCrossB: S2Point): Int {
        var sign = triageSign(a, b, c, aCrossB)
        if (sign == 0) sign = expensiveSign(a, b, c)
        return sign
    }

    // This version of Sign returns +1 if the points are definitely CCW, -1 if
    // they are definitely CW, and 0 if two points are identical or the result
    // is uncertain.  Uncertain cases can be resolved, if desired, by calling
    // ExpensiveSign.
    //
    // The purpose of this method is to allow additional cheap tests to be done,
    // where possible, in order to avoid calling ExpensiveSign unnecessarily.
    fun triageSign(a: S2Point, b: S2Point, c: S2Point, aCrossB: S2Point): Int {
        // kMaxDetError is the maximum error in computing (AxB).C where all vectors
        // are unit length.  Using standard inequalities, it can be shown that
        //
        //  fl(AxB) = AxB + D where |D| <= (|AxB| + (2/sqrt(3))*|A|*|B|) * e
        //
        // where "fl()" denotes a calculation done in floating-point arithmetic,
        // |x| denotes either absolute value or the L2-norm as appropriate, and
        // e = 0.5*DBL_EPSILON.  Similarly,
        //
        //  fl(B.C) = B.C + d where |d| <= (1.5*|B.C| + 1.5*|B|*|C|) * e .
        //
        // Applying these bounds to the unit-length vectors A,B,C and neglecting
        // relative error (which does not affect the sign of the result), we get
        //
        //  fl((AxB).C) = (AxB).C + d where |d| <= (2.5 + 2/sqrt(3)) * e
        //
        // which is about 3.6548 * e, or 1.8274 * DBL_EPSILON.
        val kMaxDetError = 1.8274 * DBL_EPSILON;
        Assertions.assertPointIsUnitLength(a)
        Assertions.assertPointIsUnitLength(b)
        Assertions.assertPointIsUnitLength(c)
        val det = aCrossB.dotProd(c)

        // Double-check borderline cases in debug mode.
        Assertions.assert {
            !kDebug || (
                    abs(det) <= kMaxDetError
                            || abs(det) >= 100 * kMaxDetError
                            || det * expensiveSign(a, b, c) > 0.0
                    )
        }

        if (det > kMaxDetError) return 1;
        if (det < -kMaxDetError) return -1;
        return 0;
    }

    // This function is invoked by Sign() if the sign of the determinant is
    // uncertain.  It always returns a non-zero result unless two of the input
    // points are the same.  It uses a combination of multiple-precision
    // arithmetic and symbolic perturbations to ensure that its results are
    // always self-consistent (cf. Simulation of Simplicity, Edelsbrunner and
    // Muecke).  The basic idea is to assign an infinitesimal symbolic
    // perturbation to every possible S2Point such that no three S2Points are
    // collinear and no four S2Points are coplanar.  These perturbations are so
    // small that they do not affect the sign of any determinant that was
    // non-zero before the perturbations.  If "perturb" is false, then instead
    // the exact sign of the unperturbed input points is returned, which can be
    // zero even when all three points are distinct.
    //
    // Unlike Sign(), this method does not require the input points to be
    // normalized.
    fun expensiveSign(a: S2Point, b: S2Point, c: S2Point, perturb: Boolean = true): Int {
        // Return zero if and only if two points are the same.  This ensures (1).
        if (a == b || b == c || c == a) return 0;

        // Next we try recomputing the determinant still using floating-point
        // arithmetic but in a more precise way.  This is more expensive than the
        // simple calculation done by TriageSign(), but it is still *much* cheaper
        // than using arbitrary-precision arithmetic.  This optimization is able to
        // compute the correct determinant sign in virtually all cases except when
        // the three points are truly collinear (e.g., three points on the equator).
        val det_sign = stableSign(a, b, c)
        if (det_sign != 0) return det_sign;

        // TODO(ericv): Create a templated version of StableSign so that we can
        // retry in "long double" precision before falling back to ExactFloat.

        // TODO(ericv): Optimize ExactFloat so that it stores up to 32 bytes of
        // mantissa inline (without requiring memory allocation).

        // Otherwise fall back to exact arithmetic and symbolic permutations.
        return exactSign(a, b, c, perturb);
    }

    // Compute the determinant in a numerically stable way.  Unlike TriageSign(),
    // this method can usually compute the correct determinant sign even when all
    // three points are as collinear as possible.  For example if three points are
    // spaced 1km apart along a random line on the Earth's surface using the
    // nearest representable points, there is only a 0.4% chance that this method
    // will not be able to find the determinant sign.  The probability of failure
    // decreases as the points get closer together; if the collinear points are
    // 1 meter apart, the failure rate drops to 0.0004%.
    //
    // This method could be extended to also handle nearly-antipodal points (and
    // in fact an earlier version of this code did exactly that), but antipodal
    // points are rare in practice so it seems better to simply fall back to
    // exact arithmetic in that case.
    private fun stableSign(a: S2Point, b: S2Point, c: S2Point): Int {
        val ab = b - a
        val bc = c - b
        val ca = a - c
        val ab2 = ab.norm2()
        val bc2 = bc.norm2()
        val ca2 = ca.norm2()

        // Now compute the determinant ((A-C)x(B-C)).C, where the vertices have been
        // cyclically permuted if necessary so that AB is the longest edge.  (This
        // minimizes the magnitude of cross product.)  At the same time we also
        // compute the maximum error in the determinant.  Using a similar technique
        // to the one used for kMaxDetError, the error is at most
        //
        //   |d| <= (3 + 6/sqrt(3)) * |A-C| * |B-C| * e
        //
        // where e = 0.5 * DBL_EPSILON.  If the determinant magnitude is larger than
        // this value then we know its sign with certainty.
        val kDetErrorMultiplier = 3.2321 * DBL_EPSILON;  // see above
        val det: Double
        val max_error: Double
        if (ab2 >= bc2 && ab2 >= ca2) {
            // AB is the longest edge, so compute (A-C)x(B-C).C.
            det = -(ca.crossProd(bc).dotProd(c));
            max_error = kDetErrorMultiplier * sqrt(ca2 * bc2);
        } else if (bc2 >= ca2) {
            // BC is the longest edge, so compute (B-A)x(C-A).A.
            det = -(ab.crossProd(ca).dotProd(a));
            max_error = kDetErrorMultiplier * sqrt(ab2 * ca2);
        } else {
            // CA is the longest edge, so compute (C-B)x(A-B).B.
            det = -(bc.crossProd(ab).dotProd(b));
            max_error = kDetErrorMultiplier * sqrt(bc2 * ab2);
        }
        return if(abs(det) <= max_error) 0 else if(det > 0) 1 else -1
    }

    // Compute the determinant using exact arithmetic and/or symbolic
    // permutations.  Requires that the three points are distinct.
    private fun exactSign(a: S2Point, b: S2Point, c: S2Point, perturb: Boolean): Int {
        assert(a != b && b != c && c != a)

        // Sort the three points in lexicographic order, keeping track of the sign
        // of the permutation.  (Each exchange inverts the sign of the determinant.)
        var perm_sign = 1
        var pa: S2Point = a
        var pb: S2Point = b
        var pc: S2Point = c
        if (pa > pb) { val temp = pa; pa = pb; pb = temp; perm_sign = -perm_sign; }
        if (pb > pc) { val temp = pb; pb = pc; pc = temp; perm_sign = -perm_sign; }
        if (pa > pb) { val temp = pa; pa = pb; pb = temp; perm_sign = -perm_sign; }
        assert(pa < pb && pb < pc)

        // Construct multiple-precision versions of the sorted points and compute
        // their exact 3x3 determinant.
        val xa = R3VectorExactFloat(pa.x, pa.y, pa.z)
        val xb = R3VectorExactFloat(pb.x, pb.y, pb.z)
        val xc = R3VectorExactFloat(pc.x, pc.y, pc.z)
        val xb_cross_xc = xb.crossProd(xc)
        val det = xa.dotProd(xb_cross_xc)

        // If the exact determinant is non-zero, we're done.
        var det_sign = det.signum()
        if (det_sign == 0 && perturb) {
            // Otherwise, we need to resort to symbolic perturbations to resolve the
            // sign of the determinant.
            det_sign = symbolicallyPerturbedSign(xa, xb, xc, xb_cross_xc)
            assert(0 != det_sign)
        }
        return perm_sign * det_sign;
    }
    
    // The following function returns the sign of the determinant of three points
    // A, B, C under a model where every possible S2Point is slightly perturbed by
    // a unique infinitesmal amount such that no three perturbed points are
    // collinear and no four points are coplanar.  The perturbations are so small
    // that they do not change the sign of any determinant that was non-zero
    // before the perturbations, and therefore can be safely ignored unless the
    // determinant of three points is exactly zero (using multiple-precision
    // arithmetic).
    //
    // Since the symbolic perturbation of a given point is fixed (i.e., the
    // perturbation is the same for all calls to this method and does not depend
    // on the other two arguments), the results of this method are always
    // self-consistent.  It will never return results that would correspond to an
    // "impossible" configuration of non-degenerate points.
    //
    // Requirements:
    //   The 3x3 determinant of A, B, C must be exactly zero.
    //   The points must be distinct, with A < B < C in lexicographic order.
    //
    // Returns:
    //   +1 or -1 according to the sign of the determinant after the symbolic
    // perturbations are taken into account.
    //
    // Reference:
    //   "Simulation of Simplicity" (Edelsbrunner and Muecke, ACM Transactions on
    //   Graphics, 1990).
    //
    private fun symbolicallyPerturbedSign(a: R3VectorExactFloat, b: R3VectorExactFloat, c: R3VectorExactFloat, b_cross_c: R3VectorExactFloat): Int {
        // This method requires that the points are sorted in lexicographically
        // increasing order.  This is because every possible S2Point has its own
        // symbolic perturbation such that if A < B then the symbolic perturbation
        // for A is much larger than the perturbation for B.
        //
        // Alternatively, we could sort the points in this method and keep track of
        // the sign of the permutation, but it is more efficient to do this before
        // converting the inputs to the multi-precision representation, and this
        // also lets us re-use the result of the cross product B x C.
        assert(a < b && b < c)

        // Every input coordinate x[i] is assigned a symbolic perturbation dx[i].
        // We then compute the sign of the determinant of the perturbed points,
        // i.e.
        //               | a[0]+da[0]  a[1]+da[1]  a[2]+da[2] |
        //               | b[0]+db[0]  b[1]+db[1]  b[2]+db[2] |
        //               | c[0]+dc[0]  c[1]+dc[1]  c[2]+dc[2] |
        //
        // The perturbations are chosen such that
        //
        //   da[2] > da[1] > da[0] > db[2] > db[1] > db[0] > dc[2] > dc[1] > dc[0]
        //
        // where each perturbation is so much smaller than the previous one that we
        // don't even need to consider it unless the coefficients of all previous
        // perturbations are zero.  In fact, it is so small that we don't need to
        // consider it unless the coefficient of all products of the previous
        // perturbations are zero.  For example, we don't need to consider the
        // coefficient of db[1] unless the coefficient of db[2]*da[0] is zero.
        //
        // The follow code simply enumerates the coefficients of the perturbations
        // (and products of perturbations) that appear in the determinant above, in
        // order of decreasing perturbation magnitude.  The first non-zero
        // coefficient determines the sign of the result.  The easiest way to
        // enumerate the coefficients in the correct order is to pretend that each
        // perturbation is some tiny value "eps" raised to a power of two:
        //
        // eps**    1      2      4      8     16     32     64     128    256
        //        da[2]  da[1]  da[0]  db[2]  db[1]  db[0]  dc[2]  dc[1]  dc[0]
        //
        // Essentially we can then just count in binary and test the corresponding
        // subset of perturbations at each step.  So for example, we must test the
        // coefficient of db[2]*da[0] before db[1] because eps**12 > eps**16.
        //
        // Of course, not all products of these perturbations appear in the
        // determinant above, since the determinant only contains the products of
        // elements in distinct rows and columns.  Thus we don't need to consider
        // da[2]*da[1], db[1]*da[1], etc.  Furthermore, sometimes different pairs of
        // perturbations have the same coefficient in the determinant; for example,
        // da[1]*db[0] and db[1]*da[0] have the same coefficient (c[2]).  Therefore
        // we only need to test this coefficient the first time we encounter it in
        // the binary order above (which will be db[1]*da[0]).
        //
        // The sequence of tests below also appears in Table 4-ii of the paper
        // referenced above, if you just want to look it up, with the following
        // translations: [a,b,c] -> [i,j,k] and [0,1,2] -> [1,2,3].  Also note that
        // some of the signs are different because the opposite cross product is
        // used (e.g., B x C rather than C x B).

        var det_sign = b_cross_c[2].signum();            // da[2]
        if (det_sign != 0) return det_sign;
        det_sign = b_cross_c[1].signum();                // da[1]
        if (det_sign != 0) return det_sign;
        det_sign = b_cross_c[0].signum();                // da[0]
        if (det_sign != 0) return det_sign;

        det_sign = (c[0]*a[1] - c[1]*a[0]).signum();     // db[2]
        if (det_sign != 0) return det_sign;
        det_sign = c[0].signum();                        // db[2] * da[1]
        if (det_sign != 0) return det_sign;
        det_sign = -(c[1].signum());                     // db[2] * da[0]
        if (det_sign != 0) return det_sign;
        det_sign = (c[2]*a[0] - c[0]*a[2]).signum();     // db[1]
        if (det_sign != 0) return det_sign;
        det_sign = c[2].signum();                        // db[1] * da[0]
        if (det_sign != 0) return det_sign;
        // The following test is listed in the paper, but it is redundant because
        // the previous tests guarantee that C == (0, 0, 0).
        assert(0 == (c[1]*a[2] - c[2]*a[1]).signum());  // db[0]

        det_sign = (a[0]*b[1] - a[1]*b[0]).signum();     // dc[2]
        if (det_sign != 0) return det_sign;
        det_sign = -(b[0].signum());                     // dc[2] * da[1]
        if (det_sign != 0) return det_sign;
        det_sign = b[1].signum();                        // dc[2] * da[0]
        if (det_sign != 0) return det_sign;
        det_sign = a[0].signum();                        // dc[2] * db[1]
        if (det_sign != 0) return det_sign;
        return 1;                                     // dc[2] * db[1] * da[0]
    }

}
