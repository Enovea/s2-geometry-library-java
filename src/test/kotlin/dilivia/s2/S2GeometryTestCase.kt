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

import com.google.common.base.Splitter
import com.google.common.collect.ImmutableList
import com.google.common.collect.Iterables
import com.google.common.collect.Lists
import com.google.common.geometry.S2
import com.google.common.geometry.S2.M_PI
import com.google.common.geometry.S2Loop
import com.google.common.geometry.S2Polygon
import com.google.common.geometry.S2Polyline
import dilivia.s2.S1Angle.Companion.radians
import dilivia.s2.S2Cap.Companion.fromCenterArea
import dilivia.s2.S2CellId.Companion.fromFacePosLevel
import dilivia.s2.S2LatLng.Companion.fromDegrees
import dilivia.s2.S2Point.Companion.crossProd
import dilivia.s2.S2Point.Companion.normalize
import dilivia.s2.S2Point.Companion.plus
import dilivia.s2.S2Point.Companion.times
import junit.framework.TestCase
import java.util.*
import kotlin.math.cos
import kotlin.math.sin
import kotlin.math.sqrt

@Strictfp
abstract class S2GeometryTestCase : TestCase() {

    @JvmField
    var rand: Random? = null
    override fun setUp() {
        rand = Random(123456)
    }

    @JvmOverloads
    fun assertDoubleNear(a: Double, b: Double, error: Double = 1e-9) {
        assertTrue("a (=$a) + error (=$error) is not > b (=$b)", a + error > b)
        assertTrue(a < b + error)
    }


    fun <T: Comparable<T>> assertLessThan(v1: T, v2: T) {
        assertTrue(v1 < v2)
    }

    fun <T: Comparable<T>> assertLessOrEquals(v1: T, v2: T) {
        assertTrue(v1 <= v2)
    }

    fun <T: Comparable<T>> assertGreaterThan(v1: T, v2: T) {
        assertTrue(v1 > v2)
    }

    fun kmToAngle(km: Double): S1Angle {
        return radians(km / kEarthRadiusKm)
    }
    // maybe these should be put in a special testing util class
    /** Return a random unit-length vector.  */
    fun randomPoint(): S2Point {
        return normalize(S2Point(
                2 * rand!!.nextDouble() - 1,
                2 * rand!!.nextDouble() - 1,
                2 * rand!!.nextDouble() - 1))
    }

    fun randomPoint(cap: S2Cap): S2Point {
        // We consider the cap axis to be the "z" axis.  We choose two other axes to
        // complete the coordinate frame.
        val m = S2Point.getFrame(cap.center)

        // The surface area of a spherical cap is directly proportional to its
        // height.  First we choose a random height, and then we choose a random
        // point along the circle at that height.
        val rnd = rand!!
        val h = rnd.nextDouble() * cap.height
        val theta = 2 * M_PI * rnd.nextDouble()
        val r = sqrt(h * (2 - h))  // Radius of circle.

        // The result should already be very close to unit-length, but we might as
        // well make it accurate as possible.
        return S2Point.fromFrame(m, S2Point(cos(theta) * r, sin(theta) * r, 1 - h)).normalize()
    }
    /**
     * Return a right-handed coordinate frame (three orthonormal vectors). Returns
     * an array of three points: x,y,z
     */
    val randomFrame: ImmutableList<S2Point>
        get() {
            val p0 = randomPoint()
            val p1 = normalize(crossProd(p0, randomPoint()))
            val p2 = normalize(crossProd(p0, p1))
            return ImmutableList.of(p0, p1, p2)
        }

    /**
     * Return a random cell id at the given level or at a randomly chosen level.
     * The distribution is uniform over the space of cell ids, but only
     * approximately uniform over the surface of the sphere.
     */
    fun getRandomCellId(level: Int): S2CellId {
        val face = random(S2CellId.kNumFaces)
        val pos = rand!!.nextLong() and (1L shl 2 * S2CellId.kMaxLevel) - 1
        return fromFacePosLevel(face, pos.toULong(), level)
    }

    val randomCellId: S2CellId
        get() = getRandomCellId(random(S2CellId.kMaxLevel + 1))

    protected fun random(n: Int): Int {
        return if (n == 0) {
            0
        } else rand!!.nextInt(n)
    }

    // Pick "base" uniformly from range [0,maxLog] and then return
    // "base" random bits. The effect is to pick a number in the range
    // [0,2^maxLog-1] with bias towards smaller numbers.
    fun skewed(maxLog: Int): Int {
        val base = Math.abs(rand!!.nextInt()) % (maxLog + 1)
        // if (!base) return 0; // if 0==base, we & with 0 below.
        //
        // this distribution differs slightly from ACMRandom's Skewed,
        // since 0 occurs approximately 3 times more than 1 here, and
        // ACMRandom's Skewed never outputs 0.
        return rand!!.nextInt() and (1 shl base) - 1
    }

    /**
     * Checks that "covering" completely covers the given region. If "check_tight"
     * is true, also checks that it does not contain any cells that do not
     * intersect the given region. ("id" is only used internally.)
     */
    fun checkCovering(region: S2Region, covering: S2CellUnion, checkTight: Boolean, id: S2CellId) {
        if (!id.isValid()) {
            for (face in 0..5) {
                checkCovering(region, covering, checkTight, fromFacePosLevel(face, 0UL, 0))
            }
            return
        }
        if (!region.mayIntersect(S2Cell(id))) {
            // If region does not intersect id, then neither should the covering.
            if (checkTight) {
                assertTrue(!covering.intersects(id))
            }
        } else if (!covering.contains(id)) {
            // The region may intersect id, but we can't assert that the covering
            // intersects id because we may discover that the region does not actually
            // intersect upon further subdivision. (MayIntersect is not exact.)
            assertTrue(!region.contains(S2Cell(id)))
            assertTrue(!id.isLeaf())
            val end = id.childEnd()
            var child = id.childBegin()
            while (!child.equals(end)) {
                checkCovering(region, covering, checkTight, child)
                child = child.next()
            }
        }
    }

    fun getRandomCap(minArea: Double, maxArea: Double): S2Cap {
        val capArea = (maxArea
                * Math.pow(minArea / maxArea, rand!!.nextDouble()))
        assertTrue(capArea >= minArea && capArea <= maxArea)

        // The surface area of a cap is 2*Pi times its height.
        return fromCenterArea(randomPoint(), capArea)
    }

    fun samplePoint(cap: S2Cap): S2Point {
        // We consider the cap axis to be the "z" axis. We choose two other axes to
        // complete the coordinate frame.
        val z = cap.center
        val x = z.ortho()
        val y = crossProd(z, x)

        // The surface area of a spherical cap is directly proportional to its
        // height. First we choose a random height, and then we choose a random
        // point along the circle at that height.
        val h = rand!!.nextDouble() * cap.height
        val theta = 2 * S2.M_PI * rand!!.nextDouble()
        val r = Math.sqrt(h * (2 - h)) // Radius of circle.

        // (cos(theta)*r*x + sin(theta)*r*y + (1-h)*z).Normalize()
        return normalize(plus(
                plus(times(x, Math.cos(theta) * r), times(y, Math.sin(theta) * r)),
                times(z, 1 - h)))
    }

    companion object {
        const val kEarthRadiusKm = 6371.01
        fun parseVertices(str: String?, vertices: MutableList<S2Point?>) {
            if (str == null) {
                return
            }
            for (token in Splitter.on(',').split(str)) {
                val colon = token.indexOf(':')
                require(colon != -1) { "Illegal string:$token. Should look like '35:20'" }
                val lat = token.substring(0, colon).toDouble()
                val lng = token.substring(colon + 1).toDouble()
                vertices.add(fromDegrees(lat, lng).toPoint())
            }
        }

        @JvmStatic
        fun makePoint(str: String?): S2Point? {
            val vertices: MutableList<S2Point?> = Lists.newArrayList()
            parseVertices(str, vertices)
            return Iterables.getOnlyElement(vertices)
        }

        fun makeLoop(str: String?): S2Loop {
            val vertices: MutableList<S2Point?> = Lists.newArrayList()
            parseVertices(str, vertices)
            return S2Loop(vertices)
        }

        @JvmStatic
        fun makePolygon(str: String?): S2Polygon {
            val loops: MutableList<S2Loop> = Lists.newArrayList()
            for (token in Splitter.on(';').omitEmptyStrings().split(str)) {
                val loop = makeLoop(token)
                loop.normalize()
                loops.add(loop)
            }
            return S2Polygon(loops)
        }

        @JvmStatic
        fun makePolyline(str: String?): S2Polyline {
            val vertices: MutableList<S2Point?> = Lists.newArrayList()
            parseVertices(str, vertices)
            return S2Polyline(vertices)
        }
    }
}