/*
 * Copyright 2005 Google Inc.
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
package com.google.common.geometry;

import com.google.common.annotations.VisibleForTesting;
import dilivia.s2.S2Point;
import dilivia.s2.math.R2Vector;

public final strictfp class S2 {

    // Declare some frequently used constants
    public static final double M_PI = Math.PI;
    public static final double M_1_PI = 1.0 / Math.PI;
    public static final double M_PI_2 = Math.PI / 2.0;
    public static final double M_PI_4 = Math.PI / 4.0;
    public static final double M_SQRT1_2 = Math.sqrt(0.5);
    public static final double M_SQRT2 = Math.sqrt(2);
    public static final double M_SQRT3 = Math.sqrt(3);
    public static final double M_E = Math.E;

    // Together these flags define a cell orientation. If SWAP_MASK
    // is true, then canonical traversal order is flipped around the
    // diagonal (i.e. i and j are swapped with each other). If
    // INVERT_MASK is true, then the traversal order is rotated by 180
    // degrees (i.e. the bits of i and j are inverted, or equivalently,
    // the axis directions are reversed).
    public static final int SWAP_MASK = 0x01;
    public static final int INVERT_MASK = 0x02;

    // Number of bits in the mantissa of a double.
    private static final int EXPONENT_SHIFT = 52;
    // Mask to extract the exponent from a double.
    private static final long EXPONENT_MASK = 0x7ff0000000000000L;

    private static double getDblEpsilon() {
        double d = 1.0;
        while (1.0 + d / 2 != 1.0) {
            d /= 2;
        }
        return d;
    }

    public static double DBL_EPSILON = getDblEpsilon();



    /**
     * WARNING! This requires arbitrary precision arithmetic to be truly robust.
     * This means that for nearly colinear AB and AC, this function may return the
     * wrong answer.
     *
     * <p>
     * Like SimpleCCW(), but returns +1 if the points are counterclockwise and -1
     * if the points are clockwise. It satisfies the following conditions:
     * <p>
     * (1) RobustCCW(a,b,c) == 0 if and only if a == b, b == c, or c == a (2)
     * RobustCCW(b,c,a) == RobustCCW(a,b,c) for all a,b,c (3) RobustCCW(c,b,a)
     * ==-RobustCCW(a,b,c) for all a,b,c
     * <p>
     * In other words:
     * <p>
     * (1) The result is zero if and only if two points are the same. (2)
     * Rotating the order of the arguments does not affect the result. (3)
     * Exchanging any two arguments inverts the result.
     * <p>
     * This function is essentially like taking the sign of the determinant of
     * a,b,c, except that it has additional logic to make sure that the above
     * properties hold even when the three points are coplanar, and to deal with
     * the limitations of floating-point arithmetic.
     * <p>
     * Note: a, b and c are expected to be of unit length. Otherwise, the results
     * are undefined.
     */
    public static int robustCCW(S2Point a, S2Point b, S2Point c) {
        return robustCCW(a, b, c, S2Point.crossProd(a, b));
    }

    /**
     * A more efficient version of RobustCCW that allows the precomputed
     * cross-product of A and B to be specified.
     * <p>
     * Note: a, b and c are expected to be of unit length. Otherwise, the results
     * are undefined
     */
    public static int robustCCW(S2Point a, S2Point b, S2Point c, S2Point aCrossB) {
        // assert (isUnitLength(a) && isUnitLength(b) && isUnitLength(c));

        // There are 14 multiplications and additions to compute the determinant
        // below. Since all three points are normalized, it is possible to show
        // that the average rounding error per operation does not exceed 2**-54,
        // the maximum rounding error for an operation whose result magnitude is in
        // the range [0.5,1). Therefore, if the absolute value of the determinant
        // is greater than 2*14*(2**-54), the determinant will have the same sign
        // even if the arguments are rotated (which produces a mathematically
        // equivalent result but with potentially different rounding errors).
        final double kMinAbsValue = 1.6e-15; // 2 * 14 * 2**-54

        double det = aCrossB.dotProd(c);

        // Double-check borderline cases in debug mode.
        // assert ((Math.abs(det) < kMinAbsValue) || (Math.abs(det) > 1000 * kMinAbsValue)
        //    || (det * expensiveCCW(a, b, c) > 0));

        if (det > kMinAbsValue) {
            return 1;
        }

        if (det < -kMinAbsValue) {
            return -1;
        }

        return expensiveCCW(a, b, c);
    }

    /**
     * A relatively expensive calculation invoked by RobustCCW() if the sign of
     * the determinant is uncertain.
     */
    private static int expensiveCCW(S2Point a, S2Point b, S2Point c) {
        // Return zero if and only if two points are the same. This ensures (1).
        if (a.equals(b) || b.equals(c) || c.equals(a)) {
            return 0;
        }

        // Now compute the determinant in a stable way. Since all three points are
        // unit length and we know that the determinant is very close to zero, this
        // means that points are very nearly colinear. Furthermore, the most common
        // situation is where two points are nearly identical or nearly antipodal.
        // To get the best accuracy in this situation, it is important to
        // immediately reduce the magnitude of the arguments by computing either
        // A+B or A-B for each pair of points. Note that even if A and B differ
        // only in their low bits, A-B can be computed very accurately. On the
        // other hand we can't accurately represent an arbitrary linear combination
        // of two vectors as would be required for Gaussian elimination. The code
        // below chooses the vertex opposite the longest edge as the "origin" for
        // the calculation, and computes the different vectors to the other two
        // vertices. This minimizes the sum of the lengths of these vectors.
        //
        // This implementation is very stable numerically, but it still does not
        // return consistent results in all cases. For example, if three points are
        // spaced far apart from each other along a great circle, the sign of the
        // result will basically be random (although it will still satisfy the
        // conditions documented in the header file). The only way to return
        // consistent results in all cases is to compute the result using
        // arbitrary-precision arithmetic. I considered using the Gnu MP library,
        // but this would be very expensive (up to 2000 bits of precision may be
        // needed to store the intermediate results) and seems like overkill for
        // this problem. The MP library is apparently also quite particular about
        // compilers and compilation options and would be a pain to maintain.

        // We want to handle the case of nearby points and nearly antipodal points
        // accurately, so determine whether A+B or A-B is smaller in each case.
        double sab = (a.dotProd(b) > 0) ? -1 : 1;
        double sbc = (b.dotProd(c) > 0) ? -1 : 1;
        double sca = (c.dotProd(a) > 0) ? -1 : 1;
        S2Point vab = S2Point.plus(a, S2Point.times(b, sab));
        S2Point vbc = S2Point.plus(b, S2Point.times(c, sbc));
        S2Point vca = S2Point.plus(c, S2Point.times(a, sca));
        double dab = vab.norm2();
        double dbc = vbc.norm2();
        double dca = vca.norm2();

        // Sort the difference vectors to find the longest edge, and use the
        // opposite vertex as the origin. If two difference vectors are the same
        // length, we break ties deterministically to ensure that the symmetry
        // properties guaranteed in the header file will be true.
        double sign;
        if (dca < dbc || (dca == dbc && a.compareTo(b) < 0)) {
            if (dab < dbc || (dab == dbc && a.compareTo(c) < 0)) {
                // The "sab" factor converts A +/- B into B +/- A.
                sign = S2Point.crossProd(vab, vca).dotProd(a) * sab; // BC is longest
                // edge
            } else {
                sign = S2Point.crossProd(vca, vbc).dotProd(c) * sca; // AB is longest
                // edge
            }
        } else {
            if (dab < dca || (dab == dca && b.compareTo(c) < 0)) {
                sign = S2Point.crossProd(vbc, vab).dotProd(b) * sbc; // CA is longest
                // edge
            } else {
                sign = S2Point.crossProd(vca, vbc).dotProd(c) * sca; // AB is longest
                // edge
            }
        }
        if (sign > 0) {
            return 1;
        }
        if (sign < 0) {
            return -1;
        }

        // The points A, B, and C are numerically indistinguishable from coplanar.
        // This may be due to roundoff error, or the points may in fact be exactly
        // coplanar. We handle this situation by perturbing all of the points by a
        // vector (eps, eps**2, eps**3) where "eps" is an infinitesmally small
        // positive number (e.g. 1 divided by a googolplex). The perturbation is
        // done symbolically, i.e. we compute what would happen if the points were
        // perturbed by this amount. It turns out that this is equivalent to
        // checking whether the points are ordered CCW around the origin first in
        // the Y-Z plane, then in the Z-X plane, and then in the X-Y plane.

        int ccw =
                planarOrderedCCW(new R2Vector(a.y(), a.z()), new R2Vector(b.y(), b.z()), new R2Vector(c.y(), c.z()));
        if (ccw == 0) {
            ccw =
                    planarOrderedCCW(new R2Vector(a.z(), a.x()), new R2Vector(b.z(), b.x()), new R2Vector(c.z(), c.x()));
            if (ccw == 0) {
                ccw = planarOrderedCCW(
                        new R2Vector(a.x(), a.y()), new R2Vector(b.x(), b.y()), new R2Vector(c.x(), c.y()));
                // assert (ccw != 0);
            }
        }
        return ccw;
    }


    public static int planarCCW(R2Vector a, R2Vector b) {
        // Return +1 if the edge AB is CCW around the origin, etc.
        double sab = (a.dotProd(b) > 0) ? -1 : 1;
        R2Vector vab = a.plus(b.times(sab));
        double da = a.norm2();
        double db = b.norm2();
        double sign;
        if (da < db || (da == db && a.compareTo(b) < 0)) {
            sign = a.crossProd(vab) * sab;
        } else {
            sign = vab.crossProd(b);
        }
        if (sign > 0) {
            return 1;
        }
        if (sign < 0) {
            return -1;
        }
        return 0;
    }

    public static int planarOrderedCCW(R2Vector a, R2Vector b, R2Vector c) {
        int sum = 0;
        sum += planarCCW(a, b);
        sum += planarCCW(b, c);
        sum += planarCCW(c, a);
        if (sum > 0) {
            return 1;
        }
        if (sum < 0) {
            return -1;
        }
        return 0;
    }

    public static boolean approxEquals(double a, double b, double maxError) {
        return Math.abs(a - b) <= maxError;
    }

    public static boolean approxEquals(double a, double b) {
        return approxEquals(a, b, 1e-15);
    }

    // Don't instantiate
    private S2() {
    }
}
