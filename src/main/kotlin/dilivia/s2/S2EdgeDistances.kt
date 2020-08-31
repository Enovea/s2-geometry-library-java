package dilivia.s2

import com.google.common.geometry.S2
import kotlin.math.cos
import kotlin.math.sin

object S2EdgeDistances {

    // Like Interpolate(), except that the parameter "ax" represents the desired
    // distance from A to the result X rather than a fraction between 0 and 1.
    fun interpolateAtDistance(ax_angle: S1Angle, a: S2Point, b: S2Point): S2Point {
        val ax = ax_angle.radians

        Assertions.assertPointIsUnitLength(a)
        Assertions.assertPointIsUnitLength(b)

        // Use RobustCrossProd() to compute the tangent vector at A towards B.  The
        // result is always perpendicular to A, even if A=B or A=-B, but it is not
        // necessarily unit length.  (We effectively normalize it below.)
        val normal = S2.robustCrossProd(a, b)
        val tangent = normal.crossProd(a)
        assert(tangent != S2Point(0, 0, 0))

        // Now compute the appropriate linear combination of A and "tangent".  With
        // infinite precision the result would always be unit length, but we
        // normalize it anyway to ensure that the error is within acceptable bounds.
        // (Otherwise errors can build up when the result of one interpolation is
        // fed into another interpolation.)
        return (cos(ax) * a + (sin(ax) / tangent.norm()) * tangent).normalize()
    }

}
