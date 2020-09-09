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
package com.google.common.geometry

import com.google.common.collect.Lists
import com.google.common.collect.Maps
import com.google.common.collect.Sets
import com.google.common.geometry.S2Loop
import dilivia.s2.*
import dilivia.s2.S2CellId.Companion.begin
import dilivia.s2.S2CellId.Companion.end
import dilivia.s2.S2LatLng.Companion.fromDegrees
import dilivia.s2.S2LatLng.Companion.fromPoint
import dilivia.s2.S2Point.Companion.crossProd
import dilivia.s2.S2Point.Companion.minus
import dilivia.s2.S2Point.Companion.normalize
import dilivia.s2.S2Point.Companion.origin
import dilivia.s2.S2Point.Companion.ortho
import dilivia.s2.S2Point.Companion.plus
import dilivia.s2.S2Point.Companion.times
import junit.framework.TestCase
import java.util.logging.Logger

/**
 * Tests for [S2Loop].
 *
 * Note that testLoopRelations2() is suppressed because it fails in corner
 * cases due to a problem with S2.robustCCW().
 *
 */
@Strictfp
class S2LoopTest : GeometryTestCase() {
    // A stripe that slightly over-wraps the equator.
    private val candyCane = makeLoop("-20:150, -20:-70, 0:70, 10:-150, 10:70, -10:-70")

    // A small clockwise loop in the northern & eastern hemisperes.
    private val smallNeCw = makeLoop("35:20, 45:20, 40:25")

    // Loop around the north pole at 80 degrees.
    private val arctic80 = makeLoop("80:-150, 80:-30, 80:90")

    // Loop around the south pole at 80 degrees.
    private val antarctic80 = makeLoop("-80:120, -80:0, -80:-120")

    // The northern hemisphere, defined using two pairs of antipodal points.
    private var northHemi = makeLoop("0:-180, 0:-90, 0:0, 0:90")

    // The northern hemisphere, defined using three points 120 degrees apart.
    private val northHemi3 = makeLoop("0:-180, 0:-60, 0:60")

    // The western hemisphere, defined using two pairs of antipodal points.
    private var westHemi = makeLoop("0:-180, -90:0, 0:0, 90:0")

    // The "near" hemisphere, defined using two pairs of antipodal points.
    private val nearHemi = makeLoop("0:-90, -90:0, 0:90, 90:0")

    // A diamond-shaped loop around the point 0:180.
    private val loopA = makeLoop("0:178, -1:180, 0:-179, 1:-180")

    // Another diamond-shaped loop around the point 0:180.
    private val loopB = makeLoop("0:179, -1:180, 0:-178, 1:-180")

    // The intersection of A and B.
    private val aIntersectB = makeLoop("0:179, -1:180, 0:-179, 1:-180")

    // The union of A and B.
    private val aUnionB = makeLoop("0:178, -1:180, 0:-178, 1:-180")

    // A minus B (concave)
    private val aMinusB = makeLoop("0:178, -1:180, 0:179, 1:-180")

    // B minus A (concave)
    private val bMinusA = makeLoop("0:-179, -1:180, 0:-178, 1:-180")

    // A self-crossing loop with a duplicated vertex
    private val bowtie = makeLoop("0:0, 2:0, 1:1, 0:2, 2:2, 1:1")

    // Initialized below.
    private lateinit var southHemi: S2Loop
    private lateinit var eastHemi: S2Loop
    private lateinit var farHemi: S2Loop
    public override fun setUp() {
        super.setUp()
        southHemi = S2Loop(northHemi)
        southHemi.invert()
        eastHemi = S2Loop(westHemi)
        eastHemi.invert()
        farHemi = S2Loop(nearHemi)
        farHemi.invert()
    }

    fun testBounds() {
        assertTrue(candyCane.rectBound.lng.isFull)
        assertTrue(candyCane.rectBound.latLo().degrees() < -20)
        assertTrue(candyCane.rectBound.latHi().degrees() > 10)
        assertTrue(smallNeCw.rectBound.isFull)
        assertEquals(arctic80.rectBound,
                S2LatLngRect(fromDegrees(80, -180), fromDegrees(90, 180)))
        assertEquals(antarctic80.rectBound,
                S2LatLngRect(fromDegrees(-90, -180), fromDegrees(-80, 180)))
        arctic80.invert()
        // The highest latitude of each edge is attained at its midpoint.
        val mid = times(plus(arctic80.vertex(0), arctic80.vertex(1)), 0.5)
        assertDoubleNear(arctic80.rectBound.latHi().radians, fromPoint(mid).lat().radians)
        arctic80.invert()
        assertTrue(southHemi.rectBound.lng.isFull)
        TestCase.assertEquals(southHemi.rectBound.lat, R1Interval(-S2.M_PI_2, 0.0))
    }

    fun testAreaCentroid() {
        assertDoubleNear(northHemi.area, 2 * S2.M_PI)
        assertDoubleNear(eastHemi.area, 2 * S2.M_PI)

        // Construct spherical caps of random height, and approximate their boundary
        // with closely spaces vertices. Then check that the area and centroid are
        // correct.
        for (i in 0..99) {
            // Choose a coordinate frame for the spherical cap.
            val x = randomPoint()
            val y = normalize(crossProd(x, randomPoint()))
            val z = normalize(crossProd(x, y))

            // Given two points at latitude phi and whose longitudes differ by dtheta,
            // the geodesic between the two points has a maximum latitude of
            // atan(tan(phi) / cos(dtheta/2)). This can be derived by positioning
            // the two points at (-dtheta/2, phi) and (dtheta/2, phi).
            //
            // We want to position the vertices close enough together so that their
            // maximum distance from the boundary of the spherical cap is kMaxDist.
            // Thus we want fabs(atan(tan(phi) / cos(dtheta/2)) - phi) <= kMaxDist.
            val kMaxDist = 1e-6
            val height = 2 * rand!!.nextDouble()
            val phi = Math.asin(1 - height)
            var maxDtheta = 2 * Math.acos(Math.tan(Math.abs(phi)) / Math.tan(Math.abs(phi) + kMaxDist))
            maxDtheta = Math.min(S2.M_PI, maxDtheta) // At least 3 vertices.
            val vertices: MutableList<S2Point> = Lists.newArrayList()
            var theta = 0.0
            while (theta < 2 * S2.M_PI) {
                val xCosThetaCosPhi = times(x, Math.cos(theta) * Math.cos(phi))
                val ySinThetaCosPhi = times(y, Math.sin(theta) * Math.cos(phi))
                val zSinPhi = times(z, Math.sin(phi))
                val sum = plus(plus(xCosThetaCosPhi, ySinThetaCosPhi), zSinPhi)
                vertices.add(sum)
                theta += rand!!.nextDouble() * maxDtheta
            }
            val loop = S2Loop(vertices)
            val areaCentroid = loop.areaAndCentroid
            val area = loop.area
            val centroid = loop.centroid
            val expectedArea = 2 * S2.M_PI * height
            assertTrue(areaCentroid.area == area)
            assertTrue(centroid.equals(areaCentroid.centroid))
            assertTrue(Math.abs(area - expectedArea) <= 2 * S2.M_PI * kMaxDist)

            // high probability
            assertTrue(Math.abs(area - expectedArea) >= 0.01 * kMaxDist)
            val expectedCentroid = times(z, expectedArea * (1 - 0.5 * height))
            assertTrue(minus(centroid, expectedCentroid).norm() <= 2 * kMaxDist)
        }
    }

    private fun rotate(loop: S2Loop?): S2Loop {
        val vertices: MutableList<S2Point> = Lists.newArrayList()
        for (i in 1..loop!!.numVertices()) {
            vertices.add(loop.vertex(i))
        }
        return S2Loop(vertices)
    }

    fun testContains() {
        assertTrue(candyCane.contains(fromDegrees(5, 71).toPoint()))
        for (i in 0..3) {
            assertTrue(northHemi.contains(S2Point(0, 0, 1)))
            assertTrue(!northHemi.contains(S2Point(0, 0, -1)))
            assertTrue(!southHemi.contains(S2Point(0, 0, 1)))
            assertTrue(southHemi.contains(S2Point(0, 0, -1)))
            assertTrue(!westHemi.contains(S2Point(0, 1, 0)))
            assertTrue(westHemi.contains(S2Point(0, -1, 0)))
            assertTrue(eastHemi.contains(S2Point(0, 1, 0)))
            assertTrue(!eastHemi.contains(S2Point(0, -1, 0)))
            northHemi = rotate(northHemi)
            southHemi = rotate(southHemi)
            eastHemi = rotate(eastHemi)
            westHemi = rotate(westHemi)
        }

        // This code checks each cell vertex is contained by exactly one of
        // the adjacent cells.
        for (level in 0..2) {
            val loops: MutableList<S2Loop> = Lists.newArrayList()
            val loopVertices: MutableList<S2Point> = Lists.newArrayList()
            val points: MutableSet<S2Point> = Sets.newHashSet()
            var id = begin(level)
            while (!id.equals(end(level))) {
                val cell = S2Cell(id)
                points.add(cell.getCenter())
                for (k in 0..3) {
                    loopVertices.add(cell.getVertex(k))
                    points.add(cell.getVertex(k))
                }
                loops.add(S2Loop(loopVertices))
                loopVertices.clear()
                id = id.next()
            }
            for (point in points) {
                var count = 0
                for (j in loops.indices) {
                    if (loops[j].contains(point)) {
                        ++count
                    }
                }
                assertEquals(count, 1)
            }
        }
    }

    private fun advance(id: S2CellId, n: Int): S2CellId {
        var id = id
        var n = n
        while (id.isValid() && --n >= 0) {
            id = id.next()
        }
        return id
    }

    private fun makeCellLoop(begin: S2CellId, end: S2CellId): S2Loop {
        // Construct a CCW polygon whose boundary is the union of the cell ids
        // in the range [begin, end). We add the edges one by one, removing
        // any edges that are already present in the opposite direction.
        val edges: MutableMap<S2Point, MutableSet<S2Point>?> = Maps.newHashMap()
        var id = begin
        while (!id.equals(end)) {
            val cell = S2Cell(id)
            for (k in 0..3) {
                val a = cell.getVertex(k)
                val b = cell.getVertex(k + 1 and 3)
                if (edges[b] == null) {
                    edges[b] = Sets.newHashSet()
                }
                // if a is in b's set, remove it and remove b's set if it's empty
                // otherwise, add b to a's set
                if (!edges[b]!!.remove(a)) {
                    if (edges[a] == null) {
                        edges[a] = Sets.newHashSet()
                    }
                    edges[a]!!.add(b)
                } else if (edges[b]!!.isEmpty()) {
                    edges.remove(b)
                }
            }
            id = id.next()
        }

        // The remaining edges form a single loop. We simply follow it starting
        // at an arbitrary vertex and build up a list of vertices.
        val vertices: MutableList<S2Point> = Lists.newArrayList()
        var p = edges.keys.iterator().next()
        while (!edges.isEmpty()) {
            assertEquals(1, edges[p]!!.size)
            val next = edges[p]!!.iterator().next()
            vertices.add(p)
            edges.remove(p)
            p = next
        }
        return S2Loop(vertices)
    }

    private fun assertRelation(
            a: S2Loop?, b: S2Loop?, containsOrCrosses: Int, intersects: Boolean, nestable: Boolean) {
        assertEquals(a!!.contains(b), containsOrCrosses == 1)
        assertEquals(a.intersects(b), intersects)
        if (nestable) {
            assertEquals(a.containsNested(b), a.contains(b))
        }
        if (containsOrCrosses >= -1) {
            assertEquals(a.containsOrCrosses(b), containsOrCrosses)
        }
    }

    fun testLoopRelations() {
        assertRelation(northHemi, northHemi, 1, true, false)
        assertRelation(northHemi, southHemi, 0, false, false)
        assertRelation(northHemi, eastHemi, -1, true, false)
        assertRelation(northHemi, arctic80, 1, true, true)
        assertRelation(northHemi, antarctic80, 0, false, true)
        assertRelation(northHemi, candyCane, -1, true, false)

        // We can't compare northHemi3 vs. northHemi or southHemi.
        assertRelation(northHemi3, northHemi3, 1, true, false)
        assertRelation(northHemi3, eastHemi, -1, true, false)
        assertRelation(northHemi3, arctic80, 1, true, true)
        assertRelation(northHemi3, antarctic80, 0, false, true)
        assertRelation(northHemi3, candyCane, -1, true, false)
        assertRelation(southHemi, northHemi, 0, false, false)
        assertRelation(southHemi, southHemi, 1, true, false)
        assertRelation(southHemi, farHemi, -1, true, false)
        assertRelation(southHemi, arctic80, 0, false, true)
        assertRelation(southHemi, antarctic80, 1, true, true)
        assertRelation(southHemi, candyCane, -1, true, false)
        assertRelation(candyCane, northHemi, -1, true, false)
        assertRelation(candyCane, southHemi, -1, true, false)
        assertRelation(candyCane, arctic80, 0, false, true)
        assertRelation(candyCane, antarctic80, 0, false, true)
        assertRelation(candyCane, candyCane, 1, true, false)
        assertRelation(nearHemi, westHemi, -1, true, false)
        assertRelation(smallNeCw, southHemi, 1, true, false)
        assertRelation(smallNeCw, westHemi, 1, true, false)
        assertRelation(smallNeCw, northHemi, -2, true, false)
        assertRelation(smallNeCw, eastHemi, -2, true, false)
        assertRelation(loopA, loopA, 1, true, false)
        assertRelation(loopA, loopB, -1, true, false)
        assertRelation(loopA, aIntersectB, 1, true, false)
        assertRelation(loopA, aUnionB, 0, true, false)
        assertRelation(loopA, aMinusB, 1, true, false)
        assertRelation(loopA, bMinusA, 0, false, false)
        assertRelation(loopB, loopA, -1, true, false)
        assertRelation(loopB, loopB, 1, true, false)
        assertRelation(loopB, aIntersectB, 1, true, false)
        assertRelation(loopB, aUnionB, 0, true, false)
        assertRelation(loopB, aMinusB, 0, false, false)
        assertRelation(loopB, bMinusA, 1, true, false)
        assertRelation(aIntersectB, loopA, 0, true, false)
        assertRelation(aIntersectB, loopB, 0, true, false)
        assertRelation(aIntersectB, aIntersectB, 1, true, false)
        assertRelation(aIntersectB, aUnionB, 0, true, true)
        assertRelation(aIntersectB, aMinusB, 0, false, false)
        assertRelation(aIntersectB, bMinusA, 0, false, false)
        assertRelation(aUnionB, loopA, 1, true, false)
        assertRelation(aUnionB, loopB, 1, true, false)
        assertRelation(aUnionB, aIntersectB, 1, true, true)
        assertRelation(aUnionB, aUnionB, 1, true, false)
        assertRelation(aUnionB, aMinusB, 1, true, false)
        assertRelation(aUnionB, bMinusA, 1, true, false)
        assertRelation(aMinusB, loopA, 0, true, false)
        assertRelation(aMinusB, loopB, 0, false, false)
        assertRelation(aMinusB, aIntersectB, 0, false, false)
        assertRelation(aMinusB, aUnionB, 0, true, false)
        assertRelation(aMinusB, aMinusB, 1, true, false)
        assertRelation(aMinusB, bMinusA, 0, false, true)
        assertRelation(bMinusA, loopA, 0, false, false)
        assertRelation(bMinusA, loopB, 0, true, false)
        assertRelation(bMinusA, aIntersectB, 0, false, false)
        assertRelation(bMinusA, aUnionB, 0, true, false)
        assertRelation(bMinusA, aMinusB, 0, false, true)
        assertRelation(bMinusA, bMinusA, 1, true, false)
    }

    /**
     * TODO(user, ericv) Fix this test. It fails sporadically.
     *
     *
     * The problem is not in this test, it is that
     * [S2.robustCCW] currently requires
     * arbitrary-precision arithmetic to be truly robust. That means it can give
     * the wrong answers in cases where we are trying to determine edge
     * intersections.
     *
     *
     * It seems the strictfp modifier here in java (required for correctness in
     * other areas of the library) restricts the size of temporary registers,
     * causing us to lose some of the precision that the C++ version gets.
     *
     *
     * This test fails when it randomly chooses a cell loop with nearly colinear
     * edges. That's where S2.robustCCW provides the wrong answer. Note that there
     * is an attempted workaround in [S2Loop.isValid], but it
     * does not cover all cases.
     */
    fun suppressedTestLoopRelations2() {
        // Construct polygons consisting of a sequence of adjacent cell ids
        // at some fixed level. Comparing two polygons at the same level
        // ensures that there are no T-vertices.
        for (iter in 0..999) {
            val num = rand!!.nextLong()
            var begin = S2CellId((num or 1).toULong())
            if (!begin.isValid()) {
                continue
            }
            begin = begin.parent(Math.round(rand!!.nextDouble() * S2CellId.kMaxLevel).toInt())
            val aBegin = advance(begin, skewed(6))
            val aEnd = advance(aBegin, skewed(6) + 1)
            val bBegin = advance(begin, skewed(6))
            val bEnd = advance(bBegin, skewed(6) + 1)
            if (!aEnd.isValid() || !bEnd.isValid()) {
                continue
            }
            val a = makeCellLoop(aBegin, aEnd)
            val b = makeCellLoop(bBegin, bEnd)
            val contained = aBegin.lessOrEquals(bBegin) && bEnd.lessOrEquals(aEnd)
            val intersects = aBegin.lessThan(bEnd) && bBegin.lessThan(aEnd)
            log.info(
                    "Checking " + a.numVertices() + " vs. " + b.numVertices() + ", contained = " + contained
                            + ", intersects = " + intersects)
            assertEquals(contained, a.contains(b))
            assertEquals(intersects, a.intersects(b))
        }
    }

    /**
     * Tests that nearly colinear points pass S2Loop.isValid()
     */
    fun testRoundingError() {
        val a = S2Point(-0.9190364081111774, 0.17231932652084575, 0.35451111445694833)
        val b = S2Point(-0.92130667053206, 0.17274500072476123, 0.3483578383756171)
        val c = S2Point(-0.9257244057938284, 0.17357332608634282, 0.3360158106235289)
        val d = S2Point(-0.9278712595449962, 0.17397586116468677, 0.32982923679138537)
        assertTrue(S2Loop.isValid(Lists.newArrayList(a, b, c, d)))
    }

    /**
     * Tests [S2Loop.isValid].
     */
    fun testIsValid() {
        assertTrue(loopA.isValid)
        assertTrue(loopB.isValid)
        assertFalse(bowtie.isValid)
    }

    /**
     * Tests [S2Loop.compareTo].
     */
    fun testComparisons() {
        val abc = makeLoop("0:1, 0:2, 1:2")
        val abcd = makeLoop("0:1, 0:2, 1:2, 1:1")
        val abcde = makeLoop("0:1, 0:2, 1:2, 1:1, 1:0")
        assertTrue(abc.compareTo(abcd) < 0)
        assertTrue(abc.compareTo(abcde) < 0)
        assertTrue(abcd.compareTo(abcde) < 0)
        assertTrue(abcd.compareTo(abc) > 0)
        assertTrue(abcde.compareTo(abc) > 0)
        assertTrue(abcde.compareTo(abcd) > 0)
        val bcda = makeLoop("0:2, 1:2, 1:1, 0:1")
        assertEquals(0, abcd.compareTo(bcda))
        assertEquals(0, bcda.compareTo(abcd))
        val wxyz = makeLoop("10:11, 10:12, 11:12, 11:11")
        assertTrue(abcd.compareTo(wxyz) > 0)
        assertTrue(wxyz.compareTo(abcd) < 0)
    }

    fun testGetDistance() {
        // Error margin since we're doing numerical computations
        val epsilon = 1e-15

        // A square with (lat,lng) vertices (0,1), (1,1), (1,2) and (0,2)
        // Tests the case where the shortest distance is along a normal to an edge,
        // onto a vertex
        val s1 = makeLoop("0:1, 1:1, 1:2, 0:2")

        // A square with (lat,lng) vertices (-1,1), (1,1), (1,2) and (-1,2)
        // Tests the case where the shortest distance is along a normal to an edge,
        // not onto a vertex
        val s2 = makeLoop("-1:1, 1:1, 1:2, -1:2")

        // A diamond with (lat,lng) vertices (1,0), (2,1), (3,0) and (2,-1)
        // Test the case where the shortest distance is NOT along a normal to an
        // edge
        val s3 = makeLoop("1:0, 2:1, 3:0, 2:-1")

        // All the vertices should be distance 0
        for (i in 0 until s1.numVertices()) {
            assertEquals(0.0, s1.getDistance(s1.vertex(i)).radians, epsilon)
        }

        // A point on one of the edges should be distance 0
        assertEquals(0.0, s1.getDistance(fromDegrees(0.5, 1.0).toPoint()).radians, epsilon)

        // In all three cases, the closest point to the origin is (0,1), which is at
        // a distance of 1 degree.
        // Note: all of these are intentionally distances measured along the
        // equator, since that makes the math significantly simpler. Otherwise, the
        // distance wouldn't actually be 1 degree.
        val origin = fromDegrees(0, 0).toPoint()
        assertEquals(1.0, s1.getDistance(origin).degrees(), epsilon)
        assertEquals(1.0, s2.getDistance(origin).degrees(), epsilon)
        assertEquals(1.0, s3.getDistance(origin).degrees(), epsilon)
    }

    /**
     * This function is useful for debugging.
     */
    private fun dumpCrossings(loop: S2Loop) {
        println("Ortho(v1): " + ortho(loop.vertex(1)))
        System.out.printf("Contains(kOrigin): %b\n", loop.contains(origin()))
        for (i in 1..loop.numVertices()) {
            val a = ortho(loop.vertex(i))
            val b = loop.vertex(i - 1)
            val c = loop.vertex(i + 1)
            val o = loop.vertex(i)
            System.out.printf("""
    Vertex %d: [%.17g, %.17g, %.17g], %d%dR=%d, %d%d%d=%d, R%d%d=%d, inside: %b
    
    """.trimIndent(),
                    i,
                    loop.vertex(i).x(),
                    loop.vertex(i).y(),
                    loop.vertex(i).z(),
                    i - 1,
                    i,
                    S2.robustCCW(b, o, a),
                    i + 1,
                    i,
                    i - 1,
                    S2.robustCCW(c, o, b),
                    i,
                    i + 1,
                    S2.robustCCW(a, o, c),
                    S2.orderedCCW(a, b, c, o))
        }
        for (i in 0 until loop.numVertices() + 2) {
            var orig: S2Point? = origin()
            var dest: S2Point?
            if (i < loop.numVertices()) {
                dest = loop.vertex(i)
                System.out.printf("Origin->%d crosses:", i)
            } else {
                dest = S2Point(0, 0, 1)
                if (i == loop.numVertices() + 1) {
                    orig = loop.vertex(1)
                }
                System.out.printf("Case %d:", i)
            }
            for (j in 0 until loop.numVertices()) {
                println(
                        " " + S2EdgeUtil.edgeOrVertexCrossing(orig, dest, loop.vertex(j), loop.vertex(j + 1)))
            }
            println()
        }
        var i = 0
        while (i <= 2) {
            System.out.printf("Origin->v1 crossing v%d->v1: ", i)
            val a = ortho(loop.vertex(1))
            val b = loop.vertex(i)
            val c = origin()
            val o = loop.vertex(1)
            System.out.printf("%d1R=%d, M1%d=%d, R1M=%d, crosses: %b\n",
                    i,
                    S2.robustCCW(b, o, a),
                    i,
                    S2.robustCCW(c, o, b),
                    S2.robustCCW(a, o, c),
                    S2EdgeUtil.edgeOrVertexCrossing(c, o, b, a))
            i += 2
        }
    }

    companion object {
        private val log = Logger.getLogger(S2LoopTest::class.java.canonicalName)
    }
}