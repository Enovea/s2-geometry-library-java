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
package dilivia.s2.region

import dilivia.s2.region.S2CellUnion.S2CellUnionTestPeer.fromVerbatimNoChecks
import dilivia.s2.S2.M_PI
import dilivia.s2.S2.M_PI_2
import dilivia.s2.*
import dilivia.s2.S2Random.oneIn
import dilivia.s2.S2Random.randomCap
import dilivia.s2.S2Random.randomCellId
import dilivia.s2.S2Random.randomDouble
import dilivia.s2.S2Random.randomInt
import dilivia.s2.S2Random.skewed
import mu.KotlinLogging
import java.lang.Math.pow
import kotlin.math.max
import kotlin.math.min

class S2CellUnionTest : S2GeometryTestCase() {

    private val logger = KotlinLogging.logger {  }

    /*
    class S2CellUnionTestPeer {
        public:
        // Creates a possibly invalid S2CellUnion without any checks.
        static S2CellUnion FromVerbatimNoChecks(vector<S2CellId> cell_ids) {
            return S2CellUnion(std::move(cell_ids), S2CellUnion::VERBATIM)
        }
    };*/

    fun testDefaultConstructor() {
        val empty = S2CellUnion(emptyList())
        assertTrue(empty.isEmpty())
    }

    fun testS2CellIdConstructor() {
        val face1_id = S2CellId.fromFace(1)
        val face1_union = S2CellUnion(listOf(face1_id))
        assertEquals(1, face1_union.numCells())
        assertEquals(face1_id, face1_union.cellId(0))
    }

    fun testWholeSphere() {
        val whole_sphere = S2CellUnion.wholeSphere()
        assertEquals(whole_sphere.leafCellsCovered(), 6UL * (1UL shl 60))
        whole_sphere.expand(0)
        assertEquals(whole_sphere, S2CellUnion.wholeSphere())
    }

    fun testDuplicateCellsNotValid() {
        val id = S2CellId.fromPoint(S2Point(1, 0, 0))
        val cell_union = fromVerbatimNoChecks(listOf(id, id))
        assertFalse(cell_union.isValid())
    }

    fun testUnsortedCellsNotValid() {
        val id = S2CellId.fromPoint(S2Point(1, 0, 0)).parent(10)
        val cell_union = fromVerbatimNoChecks(listOf(id, id.previous()))
        assertFalse(cell_union.isValid())
    }

    fun testInvalidCellIdNotValid() {
        assertFalse(S2CellId.none().isValid())
        val cell_union = fromVerbatimNoChecks(listOf(S2CellId.none()))
        assertFalse(cell_union.isValid())
    }

    fun testInvalidCellIdNotValidWithDebugFlag() {
        // Manually save and restore flag, to preserve test state in opensource
        // without gflags.
        assertFalse(S2CellId.none().isValid())
        val cell_union = fromVerbatimNoChecks(listOf(S2CellId.none()))
        assertFalse(cell_union.isValid())
    }

    fun testIsNormalized() {
        val id = S2CellId.fromPoint(S2Point(1, 0, 0)).parent(10)
        val cell_union = S2CellUnion.fromVerbatim((0..3).map { i -> id.child(i) })
        assertTrue(cell_union.isValid())
        assertFalse(cell_union.isNormalized())
    }

    fun addCells(id: S2CellId, selected: Boolean, input: MutableList<S2CellId>, expected: MutableList<S2CellId>) {
        var selected = selected
        // Decides whether to add "id" and/or some of its descendants to the
        // test case.  If "selected" is true, then the region covered by "id"
        // *must* be added to the test case (either by adding "id" itself, or
        // some combination of its descendants, or both).  If cell ids are to
        // the test case "input", then the corresponding expected result after
        // simplification is added to "expected".

        if (id == S2CellId.none()) {
            // Initial call: decide whether to add cell(s) from each face.
            for (face in 0..5) {
                addCells(S2CellId.fromFace(face), false, input, expected)
            }
            return
        }
        if (id.isLeaf()) {
            // The rnd.OneIn() call below ensures that the parent of a leaf cell
            // will always be selected (if we make it that far down the hierarchy).
            assert(selected)
            input.add(id)
            return
        }
        // The following code ensures that the probability of selecting a cell
        // at each level is approximately the same, i.e. we test normalization
        // of cells at all levels.
        if (!selected && oneIn(S2CellId.kMaxLevel - id.level())) {
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
        // (5/6 * 15 * 6) = 75.  The average number of output cells is about 6.

        // If a cell is selected, we add it to "input" with probability 5/6.
        var added = false
        if (selected && !oneIn(6)) {
            input.add(id)
            added = true
        }
        var num_children = 0
        var child = id.childBegin()
        for (pos in 0..3) {
            // If the cell is selected, on average we recurse on 4/12 = 1/3 child.
            // This intentionally may result in a cell and some of its children
            // being included in the test case.
            //
            // If the cell is not selected, on average we recurse on one child.
            // We also make sure that we do not recurse on all 4 children, since
            // then we might include all 4 children in the input case by accident
            // (in which case the expected output would not be correct).
            if (oneIn(if (selected) 12 else 4) && num_children < 3) {
                addCells(child, selected, input, expected)
                ++num_children
            }
            // If this cell was selected but the cell itself was not added, we
            // must ensure that all 4 children (or some combination of their
            // descendants) are added.
            if (selected && !added) addCells(child, selected, input, expected)

            child = child.next()
        }
    }

    fun testNormalize() {
        // Try a bunch of random test cases, and keep track of average
        // statistics for normalization (to see if they agree with the
        // analysis above).
        var in_sum = 0
        var out_sum = 0
        val kIters = 2000
        for (i in 0 until kIters) {
            S2Random.reset(i)
            val input = mutableListOf<S2CellId>()
            val expected = mutableListOf<S2CellId>()
            addCells(S2CellId.none(), false, input, expected)
            in_sum += input.size
            out_sum += expected.size
            val cellUnion = S2CellUnion(input)

            logger.trace { """
                | Test Normalize: Iteration $i
                | ----------------------------------
                | union: $cellUnion
                | expected: $expected
            """.trimMargin() }

            assertEquals(expected.size, cellUnion.numCells())
            for (t in 0 until expected.size) {
                assertEquals(expected[t], cellUnion[t])
            }

            // Test GetCapBound().
            val cap = cellUnion.capBound
            for (id in cellUnion) {
                assertTrue("Cell $id is not in cap bound $cap", cap.contains(S2Cell(id)))
            }

            // Test Contains(S2CellId) and Intersects(S2CellId).
            for (input_id in input) {
                assertTrue(cellUnion.contains(input_id))
                assertTrue(cellUnion.contains(input_id.toPoint()))
                assertTrue(cellUnion.intersects(input_id))
                if (!input_id.isFace()) {
                    assertTrue(cellUnion.intersects(input_id.parent()))
                    if (input_id.level() > 1) {
                        assertTrue(cellUnion.intersects(input_id.parent().parent()))
                        assertTrue(cellUnion.intersects(input_id.parent(0)))
                    }
                }
                if (!input_id.isLeaf()) {
                    assertTrue(cellUnion.contains(input_id.childBegin()))
                    assertTrue(cellUnion.intersects(input_id.childBegin()))
                    assertTrue(cellUnion.contains(input_id.childEnd().previous()))
                    assertTrue(cellUnion.intersects(input_id.childEnd().previous()))
                    assertTrue(cellUnion.contains(input_id.childBegin(S2CellId.kMaxLevel)))
                    assertTrue(cellUnion.intersects(input_id.childBegin(S2CellId.kMaxLevel)))
                }
            }
            for (expected_id in expected) {
                if (!expected_id.isFace()) {
                    assertTrue(!cellUnion.contains(expected_id.parent()))
                    assertTrue(!cellUnion.contains(expected_id.parent(0)))
                }
            }

            // Test Contains(S2CellUnion), Intersects(S2CellUnion), Union(),
            // Intersection(), and Difference().
            val x = mutableListOf<S2CellId>()
            val y = mutableListOf<S2CellId>()
            val x_or_y = mutableListOf<S2CellId>()
            val x_and_y = mutableListOf<S2CellId>()
            for (input_id in input) {
                val in_x = oneIn(2)
                val in_y = oneIn(2)
                if (in_x) x.add(input_id)
                if (in_y) y.add(input_id)
                if (in_x || in_y) x_or_y.add(input_id)
            }
            val xcells = S2CellUnion(x)
            val ycells = S2CellUnion(y)
            val x_or_y_expected = S2CellUnion(x_or_y)
            val x_or_y_cells = xcells.union(ycells)
            assertTrue(x_or_y_cells == x_or_y_expected)

            // Compute the intersection of "x" with each cell of "y",
            // check that this intersection is correct, and append the
            // results to x_and_y_expected.
            for (yid in ycells) {
                val ucells = xcells.intersection(yid)
                for (xid in xcells) {
                    if (xid.contains(yid)) {
                        assertTrue(ucells.numCells() == 1 && ucells[0] == yid)
                    } else if (yid.contains(xid)) {
                        assertTrue(ucells.contains(xid))
                    }
                }
                for (uid in ucells) {
                    assertTrue(xcells.contains(uid))
                    assertTrue(yid.contains(uid))
                }
                x_and_y.addAll(ucells)
            }
            val x_and_y_expected = S2CellUnion(x_and_y)
            val x_and_y_cells = xcells.intersection(ycells)
            assertTrue(x_and_y_cells == x_and_y_expected)

            val x_minus_y_cells = xcells.difference(ycells)
            val y_minus_x_cells = ycells.difference(xcells)
            assertTrue(xcells.contains(x_minus_y_cells))
            assertTrue(!x_minus_y_cells.intersects(ycells))
            assertTrue(ycells.contains(y_minus_x_cells))
            assertTrue(!y_minus_x_cells.intersects(xcells))
            assertTrue(!x_minus_y_cells.intersects(y_minus_x_cells))

            val diff_intersection_union = x_minus_y_cells.union(y_minus_x_cells).union(x_and_y_cells)
            assertTrue(diff_intersection_union == x_or_y_cells)

            val test = mutableListOf<S2CellId>()
            val dummy = mutableListOf<S2CellId>()
            addCells(S2CellId.none(), false, test, dummy)
            for (test_id in test) {
                var contains = false
                var intersects = false
                for (expected_id in expected) {
                    if (expected_id.contains(test_id)) contains = true
                    if (expected_id.intersects(test_id)) intersects = true
                }
                assertEquals(contains, cellUnion.contains(test_id))
                assertEquals(intersects, cellUnion.intersects(test_id))
            }
        }
        println(String.format("avg in %.2f, avg out %.2f\n", in_sum.toDouble() / kIters, out_sum.toDouble() / kIters))
    }

    // Return the maximum geodesic distance from "axis" to any point of
// "covering".
    fun getRadius(covering: S2CellUnion, axis: S2Point): Double {
        var max_dist = 0.0
        for (id in covering) {
            val cell = S2Cell(id)
            for (j in 0..3) {
                val a = cell.getVertex(j)
                val b = cell.getVertex(j + 1)
                var dist = 0.0
                // The maximum distance is not always attained at a cell vertex: if at
                // least one vertex is in the opposite hemisphere from "axis" then the
                // maximum may be attained along an edge.  We solve this by computing
                // the minimum distance from the edge to (-axis) instead.  We can't
                // simply do this all the time because S2::GetDistance() has
                // poor accuracy when the result is close to Pi.
                //
                // TODO(ericv): Improve S2::GetDistance() accuracy near Pi.
                if (a.angle(axis) > M_PI_2 || b.angle(axis) > M_PI_2) {
                    dist = M_PI - S2EdgeDistances.getDistance(-axis, a, b).radians
                } else {
                    dist = a.angle(axis)
                }
                max_dist = max(max_dist, dist)
            }
        }
        return max_dist
    }

    fun testExpand() {
        // This test generates coverings for caps of random sizes, expands
        // the coverings by a random radius, and then make sure that the new
        // covering covers the expanded cap.  It also makes sure that the
        // new covering is not too much larger than expected.

        val coverer = S2RegionCoverer()
        for (i in 0 until 1000) {
            val cap = randomCap(S2Cell.averageArea(S2CellId.kMaxLevel), 4 * M_PI)

            // Expand the cap area by a random factor whose log is uniformly
            // distributed between 0 and log(1e2).
            val expanded_cap = S2Cap.fromCenterHeight(cap.center, min(2.0, pow(1e2, randomDouble()) * cap.height))

            val radius = (expanded_cap.radius() - cap.radius()).radians
            val max_level_diff = randomInt(8)

            // Generate a covering for the original cap, and measure the maximum
            // distance from the cap center to any point in the covering.
            coverer.setMaxCells(1 + skewed(10))
            val covering = coverer.getCovering(cap)
            checkCovering(cap, covering, true)
            val covering_radius = getRadius(fromVerbatimNoChecks(covering.cellIds()), cap.center)

            // This code duplicates the logic in Expand(min_radius, max_level_diff)
            // that figures out an appropriate cell level to use for the expansion.
            var min_level = S2CellId.kMaxLevel
            for (id in covering) {
                min_level = min(min_level, id.level())
            }
            val expand_level = min(min_level + max_level_diff, S2CellMetrics.kMinWidth.getLevelForMinValue(radius))

            // Generate a covering for the expanded cap, and measure the new maximum
            // distance from the cap center to any point in the covering.
            covering.expand(S1Angle.radians(radius), max_level_diff)
            checkCovering(expanded_cap, covering, false)
            val expanded_covering_radius = getRadius(fromVerbatimNoChecks(covering.cellIds()), cap.center)

            // If the covering includes a tiny cell along the boundary, in theory the
            // maximum angle of the covering from the cap center can increase by up to
            // twice the maximum length of a cell diagonal.
            assertTrue(expanded_covering_radius - covering_radius <= 2 * S2CellMetrics.kMaxDiag.getValue(expand_level))
        }
    }

    fun testFromMinMax(min_id: S2CellId, max_id: S2CellId) {
        val cell_union = S2CellUnion.fromMinMax(min_id, max_id)
        val cell_ids = cell_union.cellIds()

        assertTrue(cell_ids.size > 0)
        assertEquals(min_id, cell_ids.first().rangeMin())
        assertEquals(max_id, cell_ids.last().rangeMax())
        for (i in cell_ids.indices) {
            if (i > 0) {
                assertEquals(cell_ids[i].rangeMin(), cell_ids[i - 1].rangeMax().next())
            }
        }
        assertTrue(cell_union.isNormalized())
    }

    fun testFromMinMax() {
        // Check the very first leaf cell and face cell.
        val face1_id = S2CellId.fromFace(0)
        testFromMinMax(face1_id.rangeMin(), face1_id.rangeMin())
        testFromMinMax(face1_id.rangeMin(), face1_id.rangeMax())

        // Check the very last leaf cell and face cell.
        val face5_id = S2CellId.fromFace(5)
        testFromMinMax(face5_id.rangeMin(), face5_id.rangeMax())
        testFromMinMax(face5_id.rangeMax(), face5_id.rangeMax())

        // Check random ranges of leaf cells.
        for (iter in 0 until 100) {
            var x = randomCellId(S2CellId.kMaxLevel)
            var y = randomCellId(S2CellId.kMaxLevel)
            if (x > y) {
                val temp = x
                x = y
                y = temp
            }
            testFromMinMax(x, y)
        }
    }

    fun testFromBeginEnd() {
        // Since FromMinMax() is implemented in terms of FromBeginEnd(), we
        // focus on test cases that generate an empty range.
        val id_begin = S2CellId.begin(S2CellId.kMaxLevel)
        var cell_union = S2CellUnion.fromBeginEnd(id_begin, id_begin)
        assertTrue(cell_union.isEmpty())

        // Test the full sphere.
        val id_end = S2CellId.end(S2CellId.kMaxLevel)
        cell_union = S2CellUnion.fromBeginEnd(id_begin, id_end)
        assertEquals(6, cell_union.numCells())
        for (id in cell_union) {
            assertTrue(id.isFace())
        }
    }

    fun testEmpty() {
        val empty_cell_union = S2CellUnion()
        val face1_id = S2CellId.fromFace(1)

        // Normalize()
        empty_cell_union.normalize()
        assertTrue(empty_cell_union.isEmpty())

        // Denormalize(...)
        val output = empty_cell_union.denormalize(0, 2)
        assertTrue(empty_cell_union.isEmpty())

        // Contains(...)
        assertFalse(empty_cell_union.contains(face1_id))
        assertTrue(empty_cell_union.contains(empty_cell_union))

        // Intersects(...)
        assertFalse(empty_cell_union.intersects(face1_id))
        assertFalse(empty_cell_union.intersects(empty_cell_union))

        // Union(...)
        val cell_union = empty_cell_union.union(empty_cell_union)
        assertTrue(cell_union.isEmpty())

        // Intersection(...)
        var intersection = empty_cell_union.intersection(face1_id)
        assertTrue(intersection.isEmpty())
        intersection = empty_cell_union.intersection(empty_cell_union)
        assertTrue(intersection.isEmpty())

        // Difference(...)
        val difference = empty_cell_union.difference(empty_cell_union)
        assertEquals(0, difference.numCells())

        // Expand(...)
        empty_cell_union.expand(S1Angle.radians(1), 20)
        assertTrue(empty_cell_union.isEmpty())
        empty_cell_union.expand(10)
        assertTrue(empty_cell_union.isEmpty())
    }

    fun testClear() {
        val face1_id = S2CellId.fromFace(1)
        val face1_union = S2CellUnion(face1_id)

        assertEquals(1, face1_union.numCells())
        assertEquals(1, face1_union.cellIds().size)

        face1_union.clear()
        assertEquals(0, face1_union.numCells())
        assertEquals(0, face1_union.cellIds().size)
    }


    fun testLeafCellsCovered() {
        var cell_union = S2CellUnion()
        assertEquals(0UL, cell_union.leafCellsCovered())

        val ids = mutableListOf<S2CellId>()
        // One leaf cell on face 0.
        ids.add(S2CellId.fromFace(0).childBegin(S2CellId.kMaxLevel))
        cell_union = S2CellUnion(ids)
        assertEquals(1UL, cell_union.leafCellsCovered())

        // Face 0 itself (which includes the previous leaf cell).
        ids.add(S2CellId.fromFace(0))
        cell_union = S2CellUnion(ids)
        assertEquals(1UL shl 60, cell_union.leafCellsCovered())
        // Five faces.
        cell_union.expand(0)
        assertEquals(5UL shl 60, cell_union.leafCellsCovered())
        // Whole world.
        cell_union.expand(0)
        assertEquals(6UL shl 60, cell_union.leafCellsCovered())

        // Add some disjoint cells.
        ids.add(S2CellId.fromFace(1).childBegin(1))
        ids.add(S2CellId.fromFace(2).childBegin(2))
        ids.add(S2CellId.fromFace(2).childEnd(2).previous())
        ids.add(S2CellId.fromFace(3).childBegin(14))
        ids.add(S2CellId.fromFace(4).childBegin(27))
        ids.add(S2CellId.fromFace(4).childEnd(15).previous())
        ids.add(S2CellId.fromFace(5).childBegin(30))
        cell_union = S2CellUnion(ids)
        val expected = 1UL + (1UL shl 6) + (1UL shl 30) + (1UL shl 32) + (2UL shl 56) + (1UL shl 58) + (1UL shl 60)
        assertEquals(expected, cell_union.leafCellsCovered())
    }

    fun testCellUnion() {
        val x = listOf("5/1322210002222", "5/132221001", "5/132221002", "5/132221003")
                .map { S2CellId.fromDebugString(it) }
        val y = listOf("5/13222100232000", "5/13222100232003", "5/13222100233110", "5/13222100233111")
                .map { S2CellId.fromDebugString(it) }
        val out = mutableListOf<S2CellId>()
        S2CellUnion.getIntersection(x, y, out)
        println(out)
        assertEquals(y, out)
    }

}
