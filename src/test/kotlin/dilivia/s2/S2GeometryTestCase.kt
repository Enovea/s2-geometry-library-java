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

    @JvmOverloads
    fun assertDoubleNear(a: Double, b: Double, error: Double = 1e-9) {
        assertTrue("a (=$a) + error (=$error) is not > b (=$b)", a + error > b)
        assertTrue("a (=$a) !< b (=$b) + error (=$error) = ${b + error}", a < b + error)
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