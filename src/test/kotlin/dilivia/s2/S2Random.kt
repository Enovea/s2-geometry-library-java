package dilivia.s2

import com.google.common.geometry.S2.M_PI
import kotlin.math.*
import kotlin.random.Random
import kotlin.random.nextUInt
import kotlin.random.nextULong

/**
 *
 */
object S2Random {

    private var random = Random(1)

    // Reset the generator state using the given seed.
    @JvmStatic
    fun reset(seed: Int) {
        random = Random(seed)
    }

    // Return a uniformly distributed 64-bit unsigned integer.
    @JvmStatic
    fun randomULong(): ULong = random.nextULong()

    // Return a uniformly distributed 32-bit unsigned integer.
    @JvmStatic
    fun randomUInt(): UInt = random.nextUInt()

    @JvmStatic
    fun randomDouble(): Double = random.nextDouble()

    // Return a uniformly distributed integer in the range [0,n).
    @JvmStatic
    fun randomInt(n: Int): Int = random.nextInt()

    // Return a uniformly distributed "double" in the range [min, limit).
    @JvmStatic
    fun randomDouble(min: Double, limit: Double): Double = random.nextDouble(min, limit)

    // Return true with probability 1 in n.
    @JvmStatic
    fun oneIn(n: Int): Boolean = randomInt(n) == 0

    // Skewed: pick "base" uniformly from range [0,max_log] and then
    // return "base" random bits.  The effect is to pick a number in the
    // range [0,2^max_log-1] with bias towards smaller numbers.

    // Pick "base" uniformly from range [0,maxLog] and then return
    // "base" random bits. The effect is to pick a number in the range
    // [0,2^maxLog-1] with bias towards smaller numbers.
    @JvmStatic
    fun skewed(maxLog: Int): Int {
        val base = Math.abs(random.nextInt()) % (maxLog + 1)
        // if (!base) return 0; // if 0==base, we & with 0 below.
        //
        // this distribution differs slightly from ACMRandom's Skewed,
        // since 0 occurs approximately 3 times more than 1 here, and
        // ACMRandom's Skewed never outputs 0.
        return random.nextInt() and (1 shl base) - 1
    }

    // Return a random unit-length vector.
    @JvmStatic
    fun randomPoint(): S2Point {
        // The order of evaluation of function arguments is unspecified,
        // so we may not just call S2Point with three RandDouble-based args.
        // Use temporaries to induce sequence points between calls.
        // The order of evaluation of function arguments is unspecified,
        // so we may not just call S2Point with three RandDouble-based args.
        // Use temporaries to induce sequence points between calls.
        val x: Double = randomDouble(-1.0, 1.0)
        val y: Double = randomDouble(-1.0, 1.0)
        val z: Double = randomDouble(-1.0, 1.0)
        return S2Point(x, y, z).normalize()
    }

    // Return a right-handed coordinate frame (three orthonormal vectors).
    @JvmStatic
    fun randomFrame(): Triple<S2Point, S2Point, S2Point> {
        return randomFrameAt(randomPoint())
    }

    // Given a unit-length z-axis, compute x- and y-axes such that (x,y,z) is a
    // right-handed coordinate frame (three orthonormal vectors).
    @JvmStatic
    fun randomFrameAt(z: S2Point): Triple<S2Point, S2Point, S2Point> {
        val x = z.crossProd(randomPoint()).normalize()
        val y = z.crossProd(x).normalize()
        return Triple(x, y, z)
    }

    // Return a cap with a random axis such that the log of its area is
    // uniformly distributed between the logs of the two given values.
    // (The log of the cap angle is also approximately uniformly distributed.)
    @JvmStatic
    fun randomCap(min_area: Double, max_area: Double): S2Cap {
        val capArea = max_area * (min_area / max_area).pow(random.nextDouble())
        Assertions.assertGE(capArea, min_area)
        Assertions.assertLE(capArea, max_area)

        // The surface area of a cap is 2*Pi times its height.
        return S2Cap.fromCenterArea(randomPoint(), capArea)
    }

    // Return a point chosen uniformly at random (with respect to area)
    // from the given cap.
    @JvmStatic
    fun samplePoint(cap: S2Cap): S2Point {
        // We consider the cap axis to be the "z" axis.  We choose two other axes to
        // complete the coordinate frame.
        val m = S2Point.getFrame(cap.center)

        // The surface area of a spherical cap is directly proportional to its
        // height.  First we choose a random height, and then we choose a random
        // point along the circle at that height.
        val h = random.nextDouble() * cap.height
        val theta = 2 * M_PI * random.nextDouble()
        val r = sqrt(h * (2 - h))  // Radius of circle.

        // The result should already be very close to unit-length, but we might as
        // well make it accurate as possible.
        return S2Point.fromFrame(m, S2Point(cos(theta) * r, sin(theta) * r, 1 - h))
                .normalize();
    }

    // Return a point chosen uniformly at random (with respect to area on the
    // sphere) from the given latitude-longitude rectangle.
    @JvmStatic
    fun samplePoint(rect: S2LatLngRect): S2Point {
        // First choose a latitude uniformly with respect to area on the sphere.
        val sin_lo = sin(rect.lat.lo)
        val sin_hi = sin(rect.lat.hi)
        val lat = asin(random.nextDouble(sin_lo, sin_hi))

        // Now choose longitude uniformly within the given range.
        val lng = rect.lng.lo + random.nextDouble() * rect.lng.length
        return S2LatLng.fromRadians(lat, lng).normalized().toPoint()
    }

    // Return a random cell id at the given level or at a randomly chosen
    // level.  The distribution is uniform over the space of cell ids,
    // but only approximately uniform over the surface of the sphere.
    @JvmStatic
    fun randomCellId(level: Int): S2CellId {
        val face = random.nextInt(S2CellId.kNumFaces)
        val pos = random.nextLong().toULong() and ((1UL shl S2CellId.kPosBits) - 1UL)
        return S2CellId.fromFacePosLevel(face, pos, level)
    }
    @JvmStatic
    fun randomCellId(): S2CellId {
        return randomCellId(random.nextInt(S2CellId.kMaxLevel + 1))
    }
}