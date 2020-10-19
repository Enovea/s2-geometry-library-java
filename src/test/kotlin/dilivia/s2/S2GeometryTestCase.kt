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
import com.google.common.collect.Iterables
import com.google.common.collect.Lists
import dilivia.s2.S1Angle.Companion.radians
import dilivia.s2.S2LatLng.Companion.fromDegrees
import dilivia.s2.collections.isSorted
import dilivia.s2.index.Delta
import dilivia.s2.index.Distance
import dilivia.s2.index.DistanceFactory
import dilivia.s2.region.S2Cell
import dilivia.s2.region.S2CellUnion
import dilivia.s2.region.S2Loop
import dilivia.s2.region.S2Loop.Companion.makeRegularLoop
import dilivia.s2.region.S2Polyline
import dilivia.s2.region.S2Region
import junit.framework.TestCase

@Strictfp
abstract class S2GeometryTestCase : TestCase() {

    @JvmOverloads
    fun assertDoubleNear(a: Double, b: Double, error: Double = 1e-9) {
        assertTrue("a (=$a) + error (=$error) is not > b (=$b)", a + error > b)
        assertTrue("a (=$a) !< b (=$b) + error (=$error) = ${b + error}", a < b + error)
    }


    fun <T : Comparable<T>> assertLessThan(v1: T, v2: T) {
        assertTrue(v1 < v2)
    }

    fun <T : Comparable<T>> assertLessOrEquals(v1: T, v2: T) {
        assertTrue(v1 <= v2)
    }

    fun <T : Comparable<T>> assertGreaterThan(v1: T, v2: T) {
        assertTrue(v1 > v2)
    }

    fun kmToAngle(km: Double): S1Angle {
        return radians(km / kEarthRadiusKm)
    }

    fun checkCovering(region: S2Region, covering: S2CellUnion, check_tight: Boolean, id: S2CellId = S2CellId()) {
        if (!id.isValid()) {
            for (face in 0..5) {
                checkCovering(region, covering, check_tight, S2CellId.fromFace(face))
            }
            return
        }

        if (!region.mayIntersect(S2Cell(id))) {
            // If region does not intersect id, then neither should the covering.
            if (check_tight) assertTrue(!covering.intersects(id))

        } else if (!covering.contains(id)) {
            // The region may intersect id, but we can't assert that the covering
            // intersects id because we may discover that the region does not actually
            // intersect upon further subdivision.  (MayIntersect is not exact.)
            assertTrue(!region.contains(S2Cell(id)))
            assertTrue(!id.isLeaf())
            val end = id.childEnd()
            var child = id.childBegin();
            while (child != end) {
                checkCovering(region, covering, check_tight, child)
                child = child.next()
            }
        }
    }


    companion object {
        const val kEarthRadiusKm = 6371.01


        // Compare two sets of "closest" items, where "expected" is computed via brute
        // force (i.e., considering every possible candidate) and "actual" is computed
        // using a spatial data structure.  Here "max_size" is a bound on the maximum
        // number of items, "max_distance" is a limit on the distance to any item, and
        // "max_error" is the maximum error allowed when selecting which items are
        // closest (see S2ClosestEdgeQuery::Options::max_error).
        fun <T: Distance<T>, D : Comparable<D>> checkDistanceResults(
                distanceFactory: DistanceFactory<T>,
                expected: List<Pair<T, D>>,
                actual: List<Pair<T, D>>,
                max_size: Int,
                max_distance: T,
                max_error: Delta): Boolean {
            // This is a conservative bound on the error in computing the distance from
            // the target geometry to an S2Cell.  Such errors can cause candidates to be
            // pruned from the result set even though they may be slightly closer.
            val kMaxPruningError: Delta = S1ChordAngle.radians(1e-15)
            return (checkResultSet(distanceFactory, actual, expected, max_size, max_distance, max_error, kMaxPruningError, "Missing") and /*not &&*/
                    checkResultSet(distanceFactory, expected, actual, max_size, max_distance, max_error, Delta.zero, "Extra"))
        }


        // Check that result set "x" contains all the expected results from "y", and
        // does not include any duplicate results.
        private fun <T: Distance<T>, D : Comparable<D>> checkResultSet(
                distanceFactory: DistanceFactory<T>,
                x: List<Pair<T, D>>,
                y:  List<Pair<T, D>>,
                max_size: Int,
                max_distance: T,
                max_error: Delta,
                max_pruning_error: Delta,
                label: String): Boolean {
            // Results should be sorted by distance, but not necessarily then by Id.
            assertTrue(x.isSorted { p1, p2 -> p1.first.compareTo(p2.first) })

            // Result set X should contain all the items from Y whose distance is less
            // than "limit" computed below.
            var limit = distanceFactory.zero()
            if (x.size < max_size) {
                // Result set X was not limited by "max_size", so it should contain all
                // the items up to "max_distance", except that a few items right near the
                // distance limit may be missed because the distance measurements used for
                // pruning S2Cells are not conservative.
                if (max_distance == distanceFactory.infinity()) {
                    limit = max_distance;
                } else {
                    limit = max_distance - max_pruning_error;
                }
            } else if (!x.isEmpty()) {
                // Result set X contains only the closest "max_size" items, to within a
                // tolerance of "max_error + max_pruning_error".
                limit = (x.last().first - max_error) - max_pruning_error;
            }

            var result = true
            for (yp in y) {
                // Note that this test also catches duplicate values.
                val count = x.count { xp -> xp.second == yp.second }
                if (yp.first < limit && count != 1) {
                    result = false
                    println((if (count > 1) "Duplicate" else label) + " distance = ${yp.first}, id = ${yp.second}\n")
                }
            }

            return result
        }

        fun parseVertices(str: String): List<S2Point> {
            val vertices: MutableList<S2Point> = mutableListOf()
            if (str != "") {
                for (token in Splitter.on(',').split(str)) {
                    val colon = token.indexOf(':')
                    require(colon != -1) { "Illegal string:$token. Should look like '35:20'" }
                    val lat = token.substring(0, colon).toDouble()
                    val lng = token.substring(colon + 1).toDouble()
                    vertices.add(fromDegrees(lat, lng).toPoint())
                }
            }
            return vertices
        }

        @JvmStatic
        fun makePoint(str: String): S2Point {
            val vertices = parseVertices(str)
            return Iterables.getOnlyElement(vertices)
        }
        @JvmStatic
        fun makeLoop(str: String): S2Loop {
            val vertices = parseVertices(str)
            return S2Loop(vertices)
        }

        @JvmStatic
        @JvmOverloads
        fun makePolyline(str: String, check: Boolean = true): S2Polyline {
            val vertices = parseVertices(str)
            return S2Polyline(vertices.toMutableList(), check)
        }

        @JvmStatic
        fun makeRegularPoints(center: S2Point, radius: S1Angle, num_vertices: Int): List<S2Point> {
            val loop = makeRegularLoop(center, radius, num_vertices)
            val points = ArrayList<S2Point>(loop.numVertices())
            for (i in 0 until loop.numVertices()) {
                points.add(loop.vertex(i))
            }
            return points;
        }


        fun KmToAngle(km: Double): S1Angle = S1Angle.radians(km / kEarthRadiusKm)
    }
}
