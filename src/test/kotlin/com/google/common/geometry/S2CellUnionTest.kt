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
import dilivia.s2.*
import dilivia.s2.S1Angle.Companion.radians
import dilivia.s2.S2Cap.Companion.fromCenterHeight
import dilivia.s2.S2Cell.Companion.averageArea
import dilivia.s2.S2CellId.Companion.fromFacePosLevel
import dilivia.s2.S2CellId.Companion.none
import java.util.*
import java.util.logging.Logger

@Strictfp
class S2CellUnionTest : GeometryTestCase() {
    fun testBasic() {
        logger.info("TestBasic")
        val empty = S2CellUnion()
        val ids = Lists.newArrayList<S2CellId>()
        empty.initFromCellIds(ids)
        assertEquals(0, empty.size())
        val face1Id = fromFacePosLevel(1, 0UL, 0)
        val face1Union = S2CellUnion()
        ids.add(face1Id)
        face1Union.initFromCellIds(ids)
        assertEquals(1, face1Union.size())
        assertEquals(face1Id, face1Union.cellId(0))
        val face2Id = fromFacePosLevel(2, 0UL, 0)
        val face2Union = S2CellUnion()
        val cellids = Lists.newArrayList<Long>()
        cellids.add(face2Id.id.toLong())
        face2Union.initFromIds(cellids)
        assertEquals(1, face2Union.size())
        assertEquals(face2Id, face2Union.cellId(0))
        val face1Cell = S2Cell(face1Id)
        val face2Cell = S2Cell(face2Id)
        assertTrue(face1Union.contains(face1Cell))
        assertTrue(!face1Union.contains(face2Cell))
    }

    fun testContainsCellUnion() {
        logger.info("TestContainsCellUnion")
        val randomCells: MutableSet<S2CellId> = HashSet()
        for (i in 0..99) {
            randomCells.add(getRandomCellId(S2CellId.kMaxLevel))
        }
        val union = S2CellUnion()
        union.initFromCellIds(Lists.newArrayList(randomCells))

        // Add one more
        while (!randomCells.add(getRandomCellId(S2CellId.kMaxLevel))) {
        }
        val unionPlusOne = S2CellUnion()
        unionPlusOne.initFromCellIds(Lists.newArrayList(randomCells))
        assertTrue(unionPlusOne.contains(union))
        assertFalse(union.contains(unionPlusOne))

        // Build the set of parent cells and check containment
        val parents: MutableSet<S2CellId> = HashSet()
        for (cellId in union) {
            parents.add(cellId.parent())
        }
        val parentUnion = S2CellUnion()
        parentUnion.initFromCellIds(Lists.newArrayList(parents))
        assertTrue(parentUnion.contains(union))
        assertFalse(union.contains(parentUnion))
    }

    private fun addCells(id: S2CellId, selected: Boolean, input: MutableList<S2CellId>,
                         expected: ArrayList<S2CellId>) {
        // Decides whether to add "id" and/or some of its descendants to the
        // test case. If "selected" is true, then the region covered by "id"
        // *must* be added to the test case (either by adding "id" itself, or
        // some combination of its descendants, or both). If cell ids are to
        // the test case "input", then the corresponding expected result after
        // simplification is added to "expected".
        var selected = selected
        if (id.equals(none())) {
            // Initial call: decide whether to add cell(s) from each face.
            for (face in 0..5) {
                addCells(fromFacePosLevel(face, 0UL, 0), false, input, expected)
            }
            return
        }
        if (id.isLeaf()) {
            // The rnd.OneIn() call below ensures that the parent of a leaf cell
            // will always be selected (if we make it that far down the hierarchy).
            assertTrue(selected)
            input.add(id)
            return
        }
        // The following code ensures that the probability of selecting a cell
        // at each level is approximately the same, i.e. we test normalization
        // of cells at all levels.
        if (!selected && random(S2CellId.kMaxLevel - id.level()) != 0) {
            // Once a cell has been selected, the expected output is predetermined.
            // We then make sure that cells are selected that will normalize to
            // the desired output.
            expected.add(id)
            selected = true
        }

        // With the rnd.OneIn() constants below, this function adds an average
        // of 5/6 * (kMaxLevel - level) cells to "input" where "level" is the
        // level at which the cell was first selected (level 15 on average).
        // Therefore the average number of input cells in a test case is about
        // (5/6 * 15 * 6) = 75. The average number of output cells is about 6.

        // If a cell is selected, we add it to "input" with probability 5/6.
        var added = false
        if (selected && random(6) != 0) {
            input.add(id)
            added = true
        }
        var numChildren = 0
        var child = id.childBegin()
        var pos = 0
        while (pos < 4) {

            // If the cell is selected, on average we recurse on 4/12 = 1/3 child.
            // This intentionally may result in a cell and some of its children
            // being included in the test case.
            //
            // If the cell is not selected, on average we recurse on one child.
            // We also make sure that we do not recurse on all 4 children, since
            // then we might include all 4 children in the input case by accident
            // (in which case the expected output would not be correct).
            if (random(if (selected) 12 else 4) == 0 && numChildren < 3) {
                addCells(child, selected, input, expected)
                ++numChildren
            }
            // If this cell was selected but the cell itself was not added, we
            // must ensure that all 4 children (or some combination of their
            // descendents) are added.
            if (selected && !added) {
                addCells(child, selected, input, expected)
            }
            ++pos
            child = child.next()
        }
    }

    fun testNormalize() {
        logger.info("TestNormalize")

        // Try a bunch of random test cases, and keep track of average
        // statistics for normalization (to see if they agree with the
        // analysis above).
        val cellunion = S2CellUnion()
        var inSum = 0.0
        var outSum = 0.0
        val kIters = 2000
        for (i in 0 until kIters) {
            val input = Lists.newArrayList<S2CellId>()
            val expected = Lists.newArrayList<S2CellId>()
            addCells(none(), false, input, expected)
            inSum += input.size.toDouble()
            outSum += expected.size.toDouble()
            cellunion.initFromCellIds(input)
            assertEquals(cellunion.size(), expected.size)
            assertEquals(expected, cellunion.cellIds())

            // Test GetCapBound().
            val cap = cellunion.capBound
            for (k in 0 until cellunion.size()) {
                assertTrue(cap.contains(S2Cell(cellunion.cellId(k))))
            }

            // Test Contains(S2CellId) and Intersects(S2CellId).
            for (j in input.indices) {
                assertTrue(cellunion.contains(input[j]))
                assertTrue(cellunion.intersects(input[j]))
                if (!input[j].isFace()) {
                    assertTrue(cellunion.intersects(input[j].parent()))
                    if (input[j].level() > 1) {
                        assertTrue(cellunion.intersects(input[j].parent().parent()))
                        assertTrue(cellunion.intersects(input[j].parent(0)))
                    }
                }
                if (!input[j].isLeaf()) {
                    assertTrue(cellunion.contains(input[j].childBegin()))
                    assertTrue(cellunion.intersects(input[j].childBegin()))
                    assertTrue(cellunion.contains(input[j].childEnd().previous()))
                    assertTrue(cellunion.intersects(input[j].childEnd().previous()))
                    assertTrue(cellunion.contains(input[j].childBegin(S2CellId.kMaxLevel)))
                    assertTrue(cellunion.intersects(input[j].childBegin(S2CellId.kMaxLevel)))
                }
            }
            for (j in expected.indices) {
                if (!expected[j].isFace()) {
                    assertTrue(!cellunion.contains(expected[j].parent()))
                    assertTrue(!cellunion.contains(expected[j].parent(0)))
                }
            }

            // Test contains(S2CellUnion) and intersects(S2CellUnion)
            val x = Lists.newArrayList<S2CellId>()
            val y = Lists.newArrayList<S2CellId>()
            val xOrY = Lists.newArrayList<S2CellId>()
            val xAndY = Lists.newArrayList<S2CellId>()
            for (j in input.indices) {
                val inX = random(2) == 0
                val inY = random(2) == 0
                if (inX) {
                    x.add(input[j])
                }
                if (inY) {
                    y.add(input[j])
                }
                if (inX || inY) {
                    xOrY.add(input[j])
                }
            }
            val xCells = S2CellUnion()
            val yCells = S2CellUnion()
            val xOrYExpected = S2CellUnion()
            val xAndYExpected = S2CellUnion()
            xCells.initFromCellIds(x)
            yCells.initFromCellIds(y)
            xOrYExpected.initFromCellIds(xOrY)
            val xOrYCells = S2CellUnion()
            xOrYCells.getUnion(xCells, yCells)
            assertEquals(xOrYExpected, xOrYCells)

            // Compute the intersection of "x" with each cell of "y",
            // check that this intersection is correct, and append the
            // results to xAndYExpected.
            for (j in 0 until yCells.size()) {
                val yId = yCells.cellId(j)
                val u = S2CellUnion()
                u.getIntersection(xCells, yId)
                for (k in 0 until xCells.size()) {
                    val xId = xCells.cellId(k)
                    if (xId.contains(yId)) {
                        assertEquals(1, u.size())
                        assertEquals(yId, u.cellId(0))
                    } else if (yId.contains(xId)) {
                        if (!u.contains(xId)) {
                            u.getIntersection(xCells, yId)
                        }
                        assertTrue(u.contains(xId))
                    }
                }
                for (k in 0 until u.size()) {
                    assertTrue(xCells.contains(u.cellId(k)))
                    assertTrue(yId.contains(u.cellId(k)))
                }
                xAndY.addAll(u.cellIds())
            }
            xAndYExpected.initFromCellIds(xAndY)
            val xAndYCells = S2CellUnion()
            xAndYCells.getIntersection(xCells, yCells)
            assertEquals(xAndYExpected, xAndYCells)
            val test = Lists.newArrayList<S2CellId>()
            val dummy = Lists.newArrayList<S2CellId>()
            addCells(none(), false, test, dummy)
            for (j in test.indices) {
                var contains = false
                var intersects = false
                for (k in expected.indices) {
                    if (expected[k].contains(test[j])) {
                        contains = true
                    }
                    if (expected[k].intersects(test[j])) {
                        intersects = true
                    }
                }
                assertEquals(cellunion.contains(test[j]), contains)
                assertEquals(cellunion.intersects(test[j]), intersects)
            }
        }
    }

    fun getMaxAngle(covering: S2CellUnion, axis: S2Point): Double {
        var maxAngle = 0.0
        for (i in 0 until covering.size()) {
            val cell = S2Cell(covering.cellId(i))
            val cellCap = cell.capBound
            val angle = axis.angle(cellCap.center) + cellCap.radius().radians
            maxAngle = Math.max(maxAngle, angle)
        }
        return maxAngle
    }

    fun testExpand() {
        logger.info("TestExpand")

        // This test generates coverings for caps of random sizes, and expands
        // the coverings by a random radius, and then make sure that the new
        // covering covers the expanded cap. It also makes sure that the
        // new covering is not too much larger than expected.
        val coverer = S2RegionCoverer()
        for (i in 0..999) {
            val cap = getRandomCap(averageArea(S2CellId.kMaxLevel), 4 * S2.M_PI)

            // Expand the cap by a random factor whose log is uniformly distributed
            // between 0 and log(1e2).
            val expandedCap = fromCenterHeight(cap.center, Math.min(2.0, Math.pow(1e2, rand!!.nextDouble())
                    * cap.height))
            val radius = expandedCap.radius().radians - cap.radius().radians
            val maxLevelDiff = random(8)
            val covering = S2CellUnion()
            coverer.setMaxCells(1 + skewed(10))
            coverer.getCovering(cap, covering)
            checkCovering(cap, covering, true, S2CellId())
            val maxAngle = getMaxAngle(covering, cap.center)
            var minLevel = S2CellId.kMaxLevel
            for (j in 0 until covering.size()) {
                minLevel = Math.min(minLevel, covering.cellId(j).level())
            }
            covering.expand(radians(radius), maxLevelDiff)
            checkCovering(expandedCap, covering, false, S2CellId())
            val expandLevel = Math.min(minLevel + maxLevelDiff, S2Projections.MIN_WIDTH.getMaxLevel(radius))
            val expandedMaxAngle = getMaxAngle(covering, cap.center)

            // If the covering includes a tiny cell along the boundary, in theory the
            // maximum angle of the covering from the cap axis can increase by up to
            // twice the maximum length of a cell diagonal. We allow for an increase
            // of slightly more than this because cell bounding caps are not exact.
            assertTrue(expandedMaxAngle - maxAngle <= 2.01 * S2Projections.MAX_DIAG
                    .getValue(expandLevel))
        }
    }

    fun testLeafCellsCovered() {
        val cellUnion = S2CellUnion()

        // empty union
        assertEquals(0, cellUnion.leafCellsCovered())
        val ids = Lists.newArrayList<S2CellId>()
        ids.add(fromFacePosLevel(
                0, (1L shl (S2CellId.kMaxLevel shl 1) - 1).toULong(), S2CellId.kMaxLevel))

        // One leaf on face 0.
        cellUnion.initFromCellIds(ids)
        assertEquals(1L, cellUnion.leafCellsCovered())

        // Face 0.
        ids.add(fromFacePosLevel(0, 0UL, 0))
        cellUnion.initFromCellIds(ids)
        assertEquals(1L shl 60, cellUnion.leafCellsCovered())

        // Five faces.
        cellUnion.expand(0)
        assertEquals(5L shl 60, cellUnion.leafCellsCovered())

        // Whole world.
        cellUnion.expand(0)
        assertEquals(6L shl 60, cellUnion.leafCellsCovered())

        // Add some disjoint cells.
        ids.add(fromFacePosLevel(1, 0UL, 1))
        ids.add(fromFacePosLevel(2, 0UL, 2))
        ids.add(fromFacePosLevel(2, (1L shl 60).toULong(), 2))
        ids.add(fromFacePosLevel(3, 0UL, 14))
        ids.add(fromFacePosLevel(4, (1L shl 60).toULong(), 15))
        ids.add(fromFacePosLevel(4, 0UL, 27))
        ids.add(fromFacePosLevel(5, 0UL, 30))
        cellUnion.initFromCellIds(ids)
        val expected = 1L + (1L shl 6) + (1L shl 30) + (1L shl 32) + (2L shl 56) + (1L shl 58) + (1L shl 60)
        assertEquals(expected, cellUnion.leafCellsCovered())
    }

    fun testAverageBasedArea() {
        val cellUnion = S2CellUnion()

        // empty union
        assertEquals(0.0, cellUnion.averageBasedArea())
        val ids = Lists.newArrayList<S2CellId>()
        ids.add(fromFacePosLevel(1, 0UL, 1))
        ids.add(fromFacePosLevel(5, 0UL, 30))
        cellUnion.initFromCellIds(ids)
        val expected = averageArea(S2CellId.kMaxLevel) * (1L + (1L shl 58))
        assertEquals(expected, cellUnion.averageBasedArea())
    }

    fun testApproxArea() {
        val cellUnion = S2CellUnion()

        // empty union
        assertEquals(0.0, cellUnion.approxArea())
        val ids = Lists.newArrayList<S2CellId>()
        ids.add(fromFacePosLevel(1, 0UL, 1))
        ids.add(fromFacePosLevel(5, 0UL, 30))
        cellUnion.initFromCellIds(ids)
        val expected = S2Cell(ids[0]).approxArea() + S2Cell(ids[1]).approxArea()
        assertEquals(expected, cellUnion.approxArea())
    }

    fun testExactArea() {
        val cellUnion = S2CellUnion()

        // empty union
        assertEquals(0.0, cellUnion.exactArea())
        val ids = Lists.newArrayList<S2CellId>()
        ids.add(fromFacePosLevel(1, 0UL, 1))
        ids.add(fromFacePosLevel(5, 0UL, 30))
        cellUnion.initFromCellIds(ids)
        val expected = S2Cell(ids[0]).exactArea() + S2Cell(ids[1]).exactArea()
        assertEquals(expected, cellUnion.averageBasedArea())
    }

    companion object {
        var logger = Logger.getLogger(S2CellUnionTest::class.java.name)
    }
}