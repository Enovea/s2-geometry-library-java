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
package dilivia.s2.index.point

import dilivia.s2.S1Angle
import dilivia.s2.S1ChordAngle
import dilivia.s2.S2CellId
import dilivia.s2.S2CellMetrics
import dilivia.s2.S2GeometryTestCase
import dilivia.s2.S2LatLng
import dilivia.s2.S2LatLngRect
import dilivia.s2.S2Point
import dilivia.s2.S2Random
import dilivia.s2.index.shape.MutableS2ShapeIndex
import dilivia.s2.index.S2MinDistance
import dilivia.s2.index.S2MinDistanceFactory
import dilivia.s2.index.S2MinDistanceTarget
import dilivia.s2.index.point.S2ClosestPointQuery.S2ClosestPointQueryCellTarget
import dilivia.s2.index.point.S2ClosestPointQuery.S2ClosestPointQueryEdgeTarget
import dilivia.s2.index.point.S2ClosestPointQuery.S2ClosestPointQueryPointTarget
import dilivia.s2.index.point.S2ClosestPointQuery.S2ClosestPointQueryShapeIndexTarget
import dilivia.s2.index.shape.FractalLoopShapeIndexFactory
import dilivia.s2.region.S2Cap
import dilivia.s2.region.S2Cell
import java.lang.Math.pow
import kotlin.math.min
import kotlin.math.pow
import dilivia.s2.index.point.S2ClosestPointQuery as TestQuery
import dilivia.s2.index.point.S2PointIndex as TestIndex


class S2ClosestPointQueryTest : S2GeometryTestCase() {

    fun testNoPoints() {
        val index = TestIndex<Int>()
        val query = TestQuery(index)
        val target = S2ClosestPointQueryPointTarget(S2Point(1, 0, 0))
        val results = query.findClosestPoints(target)
        assertEquals(0, results.size)
    }

    fun testManyDuplicatePoints() {
        val kNumPoints = 10000
        val kTestPoint = S2Point(1, 0, 0)
        val index = TestIndex<Int>()
        repeat(kNumPoints) { i ->
            index.add(kTestPoint, i)
        }
        val query = TestQuery(index)
        val target = S2ClosestPointQueryPointTarget(kTestPoint)
        val results = query.findClosestPoints(target)
        assertEquals(kNumPoints, results.size)
    }

    fun testEmptyTargetOptimized() {
        // Ensure that the optimized algorithm handles empty targets when a distance
        // limit is specified.
        val index = TestIndex<Int>()
        repeat(1000) { i ->
            index.add(S2Random.randomPoint(), i)
        }
        val query = TestQuery(index, TestQuery.Options())
        query.options().setMaxDistance(S1Angle.radians(1e-5))
        val targetIndex = MutableS2ShapeIndex()
        val target = S2ClosestPointQueryShapeIndexTarget(targetIndex)
        assertEquals(0, query.findClosestPoints(target).size)
    }

    fun testCirclePoints() {
        fail("Infinite loop")
        testWithIndexFactory(CirclePointIndexFactory(), kNumIndexes, kNumPoints, kNumQueries)
    }

    fun testFractalPoints() {
        testWithIndexFactory(FractalPointIndexFactory(), kNumIndexes, kNumPoints, kNumQueries)
    }

    fun testGridPoints() {
        fail("Infinite loop")
        testWithIndexFactory(GridPointIndexFactory(), kNumIndexes, kNumPoints, kNumQueries)
    }

    fun testConservativeCellDistanceIsUsed() {
        val savedSeed = randomSeed
        // These specific test cases happen to fail if max_error() is not properly
        // taken into account when measuring distances to S2PointIndex cells.  They
        // all involve S2ShapeIndexTarget, which takes advantage of max_error() to
        // optimize its distance calculation.
        for (seed in listOf(16, 586, 589, 822, 1959, 2298, 3155, 3490, 3723, 4953)) {
            randomSeed = seed
            testWithIndexFactory(FractalPointIndexFactory(), 5, 100, 10)
        }
        randomSeed = savedSeed
    }

    companion object {

        var randomSeed = 1

        // The approximate radius of S2Cap from which query points are chosen.
        val kTestCapRadius: S1Angle = KmToAngle(10.0)


        const val kNumIndexes = 10
        const val kNumPoints = 1000
        const val kNumQueries = 50

        // Use "query" to find the closest point(s) to the given target, and extract
        // the query results into the given vector.  Also verify that the results
        // satisfy the search criteria.
        fun getClosestPoints(target: S2MinDistanceTarget, query: TestQuery<Int>, results: MutableList<Pair<S2MinDistance, Int>>) {
            val queryResults = query.findClosestPoints(target)
            assertTrue(queryResults.size <= query.options().getMaxResult())
            val region = query.options().region
            if (region != null && query.options().maxDistance.value == S1ChordAngle.infinity) {
                // We can predict exactly how many points should be returned.
                assertEquals(min(query.options().getMaxResult(), query.index().numPoints()), queryResults.size)
            }
            for (result in queryResults) {
                // Check that the point satisfies the region() condition.
                if (region != null) assertTrue(region.contains(result.point()))

                // Check that it satisfies the max_distance() condition.
                assertTrue(result.distance < query.options().maxDistance)
                results.add(Pair(result.distance, result.data()))
            }
        }

        fun testFindClosestPoints(target: S2MinDistanceTarget, query: TestQuery<Int>) {
            val expected = mutableListOf<Pair<S2MinDistance, Int>>()
            val actual = mutableListOf<Pair<S2MinDistance, Int>>()
            query.options().useBruteForce = true
            getClosestPoints(target, query, expected)
            query.options().useBruteForce = false
            getClosestPoints(target, query, actual)
            assertTrue(
                    """
               max_results  = ${query.options().getMaxResult()},
               max_distance = ${query.options().maxDistance},
               max_error    = ${query.options().maxError}
            """.trimIndent(),
                    checkDistanceResults(S2MinDistanceFactory, expected, actual, query.options().getMaxResult(), query.options().maxDistance, query.options().maxError))


            if (expected.isEmpty()) return

            // Note that when options.max_error() > 0, expected[0].distance may not be
            // the minimum distance.  It is never larger by more than max_error(), but
            // the actual value also depends on max_results().
            //
            // Here we verify that GetDistance() and IsDistanceLess() return results
            // that are consistent with the max_error() setting.
            val maxError = query.options().maxError
            val minDistance = expected[0].first
            assertTrue(query.getDistance(target) <= minDistance.value + maxError)

            // Test IsDistanceLess().
            assertFalse(query.isDistanceLess(target, minDistance.value - maxError))
            assertTrue(query.isConservativeDistanceLessOrEqual(target, minDistance.value))
        }

        // (Note that every query is checked using the brute force algorithm.)
        fun testWithIndexFactory(factory: PointIndexFactory, num_indexes: Int, num_points: Int, num_queries: Int) {
            // Build a set of S2PointIndexes containing the desired geometry.
            val index_caps = mutableListOf<S2Cap>()
            val indexes = mutableListOf<TestIndex<Int>>()
            for (i in 0 until num_indexes) {
                S2Random.reset(randomSeed + i)
                index_caps.add(S2Cap.fromCenterAngle(S2Random.randomPoint(), kTestCapRadius))
                indexes.add(TestIndex())
                factory.addPoints(index_caps.last(), num_points, indexes.last())
            }
            for (i in 0 until num_queries) {
                S2Random.reset(randomSeed + i)
                val iIndex = S2Random.randomInt(num_indexes)
                val indexCap = index_caps[iIndex]

                // Choose query points from an area approximately 4x larger than the
                // geometry being tested.
                val queryRadius = indexCap.radius() * 2.0
                val queryCap = S2Cap.fromCenterAngle(indexCap.center, queryRadius)
                val query = TestQuery(indexes[iIndex])

                // Occasionally we don't set any limit on the number of result points.
                // (This may return all points if we also don't set a distance limit.)
                if (!S2Random.oneIn(5)) {
                    query.options().setMaxResult(1 + S2Random.randomInt(10))
                }
                // We set a distance limit 2/3 of the time.
                if (!S2Random.oneIn(3)) {
                    query.options().setMaxDistance(queryRadius * S2Random.randomDouble())
                }
                if (S2Random.oneIn(2)) {
                    // Choose a maximum error whose logarithm is uniformly distributed over
                    // a reasonable range, except that it is sometimes zero.
                    query.options().setMaxError(S1Angle.radians(pow(1e-4, S2Random.randomDouble()) * queryRadius.radians))
                }
                val filter_rect = S2LatLngRect.fromCenterSize(
                        S2LatLng.fromPoint(S2Random.samplePoint(queryCap)),
                        S2LatLng.fromLatLng(kTestCapRadius * S2Random.randomDouble(), kTestCapRadius * S2Random.randomDouble()))

                if (S2Random.oneIn(5)) {
                    query.options().region = filter_rect
                }
                when (val targetType = S2Random.randomInt(4)) {
                    0 -> {
                        // Find the points closest to a given point.
                        val point = S2Random.samplePoint(queryCap)
                        val target = S2ClosestPointQueryPointTarget(point)
                        testFindClosestPoints(target, query)
                    }
                    1 -> {
                        // Find the points closest to a given edge.
                        val a = S2Random.samplePoint(queryCap)
                        val b = S2Random.samplePoint(S2Cap.fromCenterAngle(a, queryRadius * 1e-4.pow(S2Random.randomDouble())))
                        val target = S2ClosestPointQueryEdgeTarget(a, b)
                        testFindClosestPoints(target, query)
                    }
                    2 -> {
                        // Find the points closest to a given cell.
                        val minLevel = S2CellMetrics.kMaxDiag.getLevelForMaxValue(queryRadius.radians)
                        val level = minLevel + S2Random.randomInt(S2CellId.kMaxLevel - minLevel + 1)
                        val a = S2Random.samplePoint(queryCap)
                        val cell = S2Cell(S2CellId.fromPoint(a).parent(level))
                        val target = S2ClosestPointQueryCellTarget(cell)
                        testFindClosestPoints(target, query)
                    }
                    else -> {
                        assert(3 == targetType)
                        val targetIndex = MutableS2ShapeIndex()
                        FractalLoopShapeIndexFactory().addEdges(indexCap, 100, targetIndex)
                        val target = S2ClosestPointQueryShapeIndexTarget(targetIndex)
                        target.setIncludeInteriors(S2Random.oneIn(2))
                        testFindClosestPoints(target, query)
                    }
                }
            }
        }

    }


}
