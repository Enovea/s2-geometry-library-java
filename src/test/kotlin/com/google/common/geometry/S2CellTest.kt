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

import com.google.common.geometry.S2.Metric
import dilivia.s2.GeometryTestCase
import dilivia.s2.S2Cell
import dilivia.s2.S2Cell.Companion.averageArea
import dilivia.s2.S2Cell.Companion.fromFacePosLevel
import dilivia.s2.S2CellId
import dilivia.s2.S2LatLng.Companion.fromPoint
import dilivia.s2.S2Point
import dilivia.s2.S2Point.Companion.crossProd
import dilivia.s2.S2Point.Companion.normalize
import dilivia.s2.S2Point.Companion.plus
import dilivia.s2.S2Point.Companion.unaryMinus
import java.util.*

@Strictfp
class S2CellTest : GeometryTestCase() {
    fun testFaces() {
        val edgeCounts: MutableMap<S2Point, Int> = HashMap()
        val vertexCounts: MutableMap<S2Point, Int> = HashMap()
        for (face in 0..5) {
            val id = S2CellId.fromFacePosLevel(face, 0UL, 0)
            val cell = S2Cell(id)
            assertEquals(cell.id(), id)
            assertEquals(cell.face(), face)
            assertEquals(cell.level().toInt(), 0)
            // Top-level faces have alternating orientations to get RHS coordinates.
            assertEquals(cell.orientation().toInt(), face and S2.SWAP_MASK)
            assertFalse(cell.isLeaf)
            for (k in 0..3) {
                val edgeRaw = cell.getEdgeRaw(k)
                println(edgeCounts.containsKey(edgeRaw).toString() + " Edge raw: " + edgeRaw + " in: " + edgeCounts)
                if (edgeCounts.containsKey(edgeRaw)) {
                    edgeCounts[edgeRaw] = edgeCounts[edgeRaw]!! + 1
                } else {
                    edgeCounts[edgeRaw] = 1
                }
                if (vertexCounts.containsKey(cell.getVertexRaw(k))) {
                    vertexCounts[cell.getVertexRaw(k)] = vertexCounts[cell
                            .getVertexRaw(k)]!! + 1
                } else {
                    vertexCounts[cell.getVertexRaw(k)] = 1
                }
                assertDoubleNear(cell.getVertexRaw(k).dotProd(edgeRaw), 0.0)
                assertDoubleNear(cell.getVertexRaw(k + 1 and 3).dotProd(
                        edgeRaw), 0.0)
                assertDoubleNear(normalize(
                        crossProd(cell.getVertexRaw(k), cell
                                .getVertexRaw(k + 1 and 3))).dotProd(cell.getEdge(k)), 1.0)
            }
        }
        // Check that edges have multiplicity 2 and vertices have multiplicity 3.
        for (i in edgeCounts.values) {
            assertEquals(2, i)
        }
        for (i in vertexCounts.values) {
            assertEquals(i, 3)
        }
    }

    class LevelStats {
        var count = 0.0
        var minArea = 100.0
        var maxArea = 0.0
        var avgArea = 0.0
        var minWidth = 100.0
        var maxWidth = 0.0
        var avgWidth = 0.0
        var minEdge = 100.0
        var maxEdge = 0.0
        var avgEdge = 0.0
        var maxEdgeAspect = 0.0
        var minDiag = 100.0
        var maxDiag = 0.0
        var avgDiag = 0.0
        var maxDiagAspect = 0.0
        var minAngleSpan = 100.0
        var maxAngleSpan = 0.0
        var avgAngleSpan = 0.0
        var minApproxRatio = 100.0
        var maxApproxRatio = 0.0
    }

    companion object {
        const val DEBUG_MODE = true
        var levelStats: MutableList<LevelStats> = ArrayList(
                S2CellId.kMaxLevel + 1)

        fun gatherStats(cell: S2Cell?) {
            val s = levelStats[cell!!.level().toInt()]
            val exactArea = cell.exactArea()
            val approxArea = cell.approxArea()
            var minEdge = 100.0
            var maxEdge = 0.0
            var avgEdge = 0.0
            var minDiag = 100.0
            var maxDiag = 0.0
            var minWidth = 100.0
            var maxWidth = 0.0
            var minAngleSpan = 100.0
            var maxAngleSpan = 0.0
            for (i in 0..3) {
                val edge = cell.getVertexRaw(i).angle(cell.getVertexRaw(i + 1 and 3))
                minEdge = Math.min(edge, minEdge)
                maxEdge = Math.max(edge, maxEdge)
                avgEdge += 0.25 * edge
                val mid = plus(cell.getVertexRaw(i), cell
                        .getVertexRaw(i + 1 and 3))
                val width = S2.M_PI_2 - mid.angle(cell.getEdgeRaw(i xor 2))
                minWidth = Math.min(width, minWidth)
                maxWidth = Math.max(width, maxWidth)
                if (i < 2) {
                    val diag = cell.getVertexRaw(i).angle(cell.getVertexRaw(i xor 2))
                    minDiag = Math.min(diag, minDiag)
                    maxDiag = Math.max(diag, maxDiag)
                    val angleSpan = cell.getEdgeRaw(i).angle(
                            unaryMinus(cell.getEdgeRaw(i xor 2)))
                    minAngleSpan = Math.min(angleSpan, minAngleSpan)
                    maxAngleSpan = Math.max(angleSpan, maxAngleSpan)
                }
            }
            s.count += 1.0
            s.minArea = Math.min(exactArea, s.minArea)
            s.maxArea = Math.max(exactArea, s.maxArea)
            s.avgArea += exactArea
            s.minWidth = Math.min(minWidth, s.minWidth)
            s.maxWidth = Math.max(maxWidth, s.maxWidth)
            s.avgWidth += 0.5 * (minWidth + maxWidth)
            s.minEdge = Math.min(minEdge, s.minEdge)
            s.maxEdge = Math.max(maxEdge, s.maxEdge)
            s.avgEdge += avgEdge
            s.maxEdgeAspect = Math.max(maxEdge / minEdge, s.maxEdgeAspect)
            s.minDiag = Math.min(minDiag, s.minDiag)
            s.maxDiag = Math.max(maxDiag, s.maxDiag)
            s.avgDiag += 0.5 * (minDiag + maxDiag)
            s.maxDiagAspect = Math.max(maxDiag / minDiag, s.maxDiagAspect)
            s.minAngleSpan = Math.min(minAngleSpan, s.minAngleSpan)
            s.maxAngleSpan = Math.max(maxAngleSpan, s.maxAngleSpan)
            s.avgAngleSpan += 0.5 * (minAngleSpan + maxAngleSpan)
            val approxRatio = approxArea / exactArea
            s.minApproxRatio = Math.min(approxRatio, s.minApproxRatio)
            s.maxApproxRatio = Math.max(approxRatio, s.maxApproxRatio)
        }

        val MAX_LEVEL = if (DEBUG_MODE) 6 else 10

        init {
            for (i in 0 until S2CellId.kMaxLevel + 1) {
                levelStats.add(LevelStats())
            }
        }
    }

    fun testSubdivide(cell: S2Cell?) {
        gatherStats(cell)
        if (cell!!.isLeaf) {
            return
        }
        val children = Array<S2Cell>(4)  { S2Cell() }
        assertTrue(cell.subdivide(children))
        var childId = cell.id().childBegin()
        var exactArea = 0.0
        var approxArea = 0.0
        var averageArea = 0.0
        var i = 0
        while (i < 4) {
            exactArea += children[i].exactArea()
            approxArea += children[i].approxArea()
            averageArea += children[i].averageArea()

            // Check that the child geometry is consistent with its cell id.
            assertEquals(children[i].id(), childId)
            assertTrue(children[i].center.approxEquals(childId.toPoint(), 1e-15))
            val direct = S2Cell(childId)
            assertEquals(children[i].face(), direct.face())
            assertEquals(children[i].level(), direct.level())
            assertEquals(children[i].orientation(), direct.orientation())
            assertEquals(children[i].centerRaw, direct.centerRaw)
            for (k in 0..3) {
                assertEquals(children[i].getVertexRaw(k), direct.getVertexRaw(k))
                assertEquals(children[i].getEdgeRaw(k), direct.getEdgeRaw(k))
            }

            // Test Contains() and MayIntersect().
            assertTrue(cell.contains(children[i]))
            assertTrue(cell.mayIntersect(children[i]))
            assertTrue(!children[i].contains(cell))
            assertTrue(cell.contains(children[i].centerRaw))
            for (j in 0..3) {
                assertTrue(cell.contains(children[i].getVertexRaw(j)))
                if (j != i) {
                    assertTrue(!children[i].contains(children[j].centerRaw))
                    assertTrue(!children[i].mayIntersect(children[j]))
                }
            }

            // Test GetCapBound and GetRectBound.
            val parentCap = cell.capBound
            val parentRect = cell.rectBound
            if (cell.contains(S2Point(0, 0, 1))
                    || cell.contains(S2Point(0, 0, -1))) {
                assertTrue(parentRect.lng.isFull)
            }
            val childCap = children[i].capBound
            val childRect = children[i].rectBound
            assertTrue(childCap.contains(children[i].center))
            assertTrue(childRect.contains(children[i].centerRaw))
            assertTrue(parentCap.contains(children[i].center))
            assertTrue(parentRect.contains(children[i]!!.centerRaw))
            for (j in 0..3) {
                assertTrue(childCap.contains(children[i]!!.getVertex(j)))
                assertTrue(childRect.contains(children[i]!!.getVertex(j)))
                assertTrue(childRect.contains(children[i]!!.getVertexRaw(j)))
                assertTrue(parentCap.contains(children[i]!!.getVertex(j)))
                if (!parentRect.contains(children[i]!!.getVertex(j))) {
                    println("cell: $cell i: $i j: $j")
                    println("Children " + i + ": " + children[i])
                    println("Parent rect: $parentRect")
                    println("Vertex raw(j) " + children[i]!!.getVertex(j))
                    println("Latlng of vertex: " + fromPoint(children[i]!!.getVertex(j)))
                    cell.rectBound
                }
                assertTrue(parentRect.contains(children[i]!!.getVertex(j)))
                if (!parentRect.contains(children[i]!!.getVertexRaw(j))) {
                    println("cell: $cell i: $i j: $j")
                    println("Children " + i + ": " + children[i])
                    println("Parent rect: $parentRect")
                    println("Vertex raw(j) " + children[i]!!.getVertexRaw(j))
                    println("Latlng of vertex: " + fromPoint(children[i]!!.getVertexRaw(j)))
                    cell.rectBound
                }
                assertTrue(parentRect.contains(children[i]!!.getVertexRaw(j)))
                if (j != i) {
                    // The bounding caps and rectangles should be tight enough so that
                    // they exclude at least two vertices of each adjacent cell.
                    var capCount = 0
                    var rectCount = 0
                    for (k in 0..3) {
                        if (childCap.contains(children[j]!!.getVertex(k))) {
                            ++capCount
                        }
                        if (childRect.contains(children[j]!!.getVertexRaw(k))) {
                            ++rectCount
                        }
                    }
                    assertTrue(capCount <= 2)
                    if (childRect.latLo().radians > -S2.M_PI_2
                            && childRect.latHi().radians < S2.M_PI_2) {
                        // Bounding rectangles may be too large at the poles because the
                        // pole itself has an arbitrary fixed longitude.
                        assertTrue(rectCount <= 2)
                    }
                }
            }

            // Check all children for the first few levels, and then sample randomly.
            // Also subdivide one corner cell, one edge cell, and one center cell
            // so that we have a better chance of sample the minimum metric values.
            var forceSubdivide = false
            val center = S2Projections.getNorm(children[i]!!.face())
            val edge = plus(center, S2Projections.getUAxis(children[i]!!.face()))
            val corner = plus(edge, S2Projections.getVAxis(children[i]!!.face()))
            for (j in 0..3) {
                val p = children[i]!!.getVertexRaw(j)
                if (p.equals(center) || p.equals(edge) || p.equals(corner)) {
                    forceSubdivide = true
                }
            }
            if (forceSubdivide || cell.level() < (if (DEBUG_MODE) 5 else 6) || random(if (DEBUG_MODE) 10 else 4) == 0) {
                testSubdivide(children[i])
            }
            ++i
            childId = childId.next()
        }

        // Check sum of child areas equals parent area.
        //
        // For ExactArea(), the best relative error we can expect is about 1e-6
        // because the precision of the unit vector coordinates is only about 1e-15
        // and the edge length of a leaf cell is about 1e-9.
        //
        // For ApproxArea(), the areas are accurate to within a few percent.
        //
        // For AverageArea(), the areas themselves are not very accurate, but
        // the average area of a parent is exactly 4 times the area of a child.
        assertTrue(Math.abs(Math.log(exactArea / cell.exactArea())) <= Math
                .abs(Math.log(1 + 1e-6)))
        assertTrue(Math.abs(Math.log(approxArea / cell.approxArea())) <= Math
                .abs(Math.log(1.03)))
        assertTrue(Math.abs(Math.log(averageArea / cell.averageArea())) <= Math
                .abs(Math.log(1 + 1e-15)))
    }

    fun testMinMaxAvg(label: String?, level: Int, count: Double,
                      absError: Double, minValue: Double, maxValue: Double, avgValue: Double,
                      minMetric: Metric, maxMetric: Metric, avgMetric: Metric) {

        // All metrics are minimums, maximums, or averages of differential
        // quantities, and therefore will not be exact for cells at any finite
        // level. The differential minimum is always a lower bound, and the maximum
        // is always an upper bound, but these minimums and maximums may not be
        // achieved for two different reasons. First, the cells at each level are
        // sampled and we may miss the most extreme examples. Second, the actual
        // metric for a cell is obtained by integrating the differential quantity,
        // which is not constant across the cell. Therefore cells at low levels
        // (bigger cells) have smaller variations.
        //
        // The "tolerance" below is an attempt to model both of these effects.
        // At low levels, error is dominated by the variation of differential
        // quantities across the cells, while at high levels error is dominated by
        // the effects of random sampling.
        var tolerance = maxMetric.getValue(level) - minMetric.getValue(level) / Math.sqrt(Math.min(count, 0.5 * (1L shl level))) * 10
        if (tolerance == 0.0) {
            tolerance = absError
        }
        val minError = minValue - minMetric.getValue(level)
        val maxError = maxMetric.getValue(level) - maxValue
        val avgError = Math.abs(avgMetric.getValue(level) - avgValue)
        System.out.printf("""
    %-10s (%6.0f samples, tolerance %8.3g) - min (%9.3g : %9.3g) max (%9.3g : %9.3g), avg (%9.3g : %9.3g)
    
    """.trimIndent(), label, count,
                tolerance, minError / minValue, minError / tolerance, maxError
                / maxValue, maxError / tolerance, avgError / avgValue, avgError
                / tolerance)
        assertTrue(minMetric.getValue(level) <= minValue + absError)
        assertTrue(minMetric.getValue(level) >= minValue - tolerance)
        println("Level: " + maxMetric.getValue(level) + " max " + (maxValue + tolerance))
        assertTrue(maxMetric.getValue(level) <= maxValue + tolerance)
        assertTrue(maxMetric.getValue(level) >= maxValue - absError)
        assertDoubleNear(avgMetric.getValue(level), avgValue, 10 * tolerance)
    }

    fun testSubdivide() {
        for (face in 0..5) {
            testSubdivide(fromFacePosLevel(face, 0.toByte(), 0))
        }

        // The maximum edge *ratio* is the ratio of the longest edge of any cell to
        // the shortest edge of any cell at the same level (and similarly for the
        // maximum diagonal ratio).
        //
        // The maximum edge *aspect* is the maximum ratio of the longest edge of a
        // cell to the shortest edge of that same cell (and similarly for the
        // maximum diagonal aspect).
        System.out
                .printf("Level    Area      Edge          Diag          Approx       Average\n")
        System.out
                .printf("        Ratio  Ratio Aspect  Ratio Aspect    Min    Max    Min    Max\n")
        for (i in 0..S2CellId.kMaxLevel) {
            val s = levelStats[i]
            if (s.count > 0) {
                s.avgArea /= s.count
                s.avgWidth /= s.count
                s.avgEdge /= s.count
                s.avgDiag /= s.count
                s.avgAngleSpan /= s.count
            }
            System.out.printf(
                    "%5d  %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f\n", i,
                    s.maxArea / s.minArea, s.maxEdge / s.minEdge, s.maxEdgeAspect,
                    s.maxDiag / s.minDiag, s.maxDiagAspect, s.minApproxRatio,
                    s.maxApproxRatio, averageArea(i) / s.maxArea, averageArea(i)
                    / s.minArea)
        }

        // Now check the validity of the S2 length and area metrics.
        for (i in 0..S2CellId.kMaxLevel) {
            val s = levelStats[i]
            if (s.count == 0.0) {
                continue
            }
            System.out.printf(
                    "Level %2d - metric (error/actual : error/tolerance)\n", i)

            // The various length calculations are only accurate to 1e-15 or so,
            // so we need to allow for this amount of discrepancy with the theoretical
            // minimums and maximums. The area calculation is accurate to about 1e-15
            // times the cell width.
            testMinMaxAvg("area", i, s.count, 1e-15 * s.minWidth, s.minArea,
                    s.maxArea, s.avgArea, S2Projections.MIN_AREA, S2Projections.MAX_AREA,
                    S2Projections.AVG_AREA)
            testMinMaxAvg("width", i, s.count, 1e-15, s.minWidth, s.maxWidth,
                    s.avgWidth, S2Projections.MIN_WIDTH, S2Projections.MAX_WIDTH,
                    S2Projections.AVG_WIDTH)
            testMinMaxAvg("edge", i, s.count, 1e-15, s.minEdge, s.maxEdge,
                    s.avgEdge, S2Projections.MIN_EDGE, S2Projections.MAX_EDGE,
                    S2Projections.AVG_EDGE)
            testMinMaxAvg("diagonal", i, s.count, 1e-15, s.minDiag, s.maxDiag,
                    s.avgDiag, S2Projections.MIN_DIAG, S2Projections.MAX_DIAG,
                    S2Projections.AVG_DIAG)
            testMinMaxAvg("angle span", i, s.count, 1e-15, s.minAngleSpan,
                    s.maxAngleSpan, s.avgAngleSpan, S2Projections.MIN_ANGLE_SPAN,
                    S2Projections.MAX_ANGLE_SPAN, S2Projections.AVG_ANGLE_SPAN)

            // The aspect ratio calculations are ratios of lengths and are therefore
            // less accurate at higher subdivision levels.
            assertTrue(s.maxEdgeAspect <= S2Projections.MAX_EDGE_ASPECT + 1e-15
                    * (1 shl i))
            assertTrue(s.maxDiagAspect <= S2Projections.MAX_DIAG_ASPECT + 1e-15
                    * (1 shl i))
        }
    }

    fun expandChildren1(cell: S2Cell?) {
        val children = Array<S2Cell>(4) { S2Cell() }
        assertTrue(cell!!.subdivide(children))
        if (children[0].level() < MAX_LEVEL) {
            for (pos in 0..3) {
                expandChildren1(children[pos])
            }
        }
    }

    fun expandChildren2(cell: S2Cell) {
        var id = cell.id().childBegin()
        var pos = 0
        while (pos < 4) {
            val child = S2Cell(id)
            if (child.level() < MAX_LEVEL) {
                expandChildren2(child)
            }
            ++pos
            id = id.next()
        }
    }
}