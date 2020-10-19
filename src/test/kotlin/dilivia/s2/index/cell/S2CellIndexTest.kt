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
package dilivia.s2.index.cell

import dilivia.s2.S2CellId
import dilivia.s2.S2GeometryTestCase
import dilivia.s2.S2Random
import dilivia.s2.index.Label
import dilivia.s2.index.S2CellIndex
import dilivia.s2.index.S2CellIndex.CellIterator
import dilivia.s2.index.S2CellIndex.ContentsIterator
import dilivia.s2.index.S2CellIndex.LabelledCell
import dilivia.s2.region.S2CellUnion
import mu.KotlinLogging
import java.util.*


class S2CellIndexTest : S2GeometryTestCase() {

    private val logger = KotlinLogging.logger {  }

    private val index = S2CellIndex()
    private val contents = mutableListOf<LabelledCell>()


    fun testEmpty() {
        quadraticValidate()
        // Deltas: [(cellId = Invalid: 0000000000000000, label-1) ; (cellId = Invalid: 0000000000000000, label-1) ; ]
    }

    fun testOneFaceCell() {
        add("0/", 0)
        quadraticValidate()
    }

    fun testOneLeafCell() {
        add("1/012301230123012301230123012301", 12)
        quadraticValidate()
    }

    fun testDuplicateValues() {
        add("0/", 0)
        add("0/", 0)
        add("0/", 1)
        add("0/", 17)
        quadraticValidate()
    }

    fun testDisjointCells() {
        add("0/", 0)
        add("3/", 0)
        quadraticValidate()
    }

    fun testNestedCells() {
        // Tests nested cells, including cases where several cells have the same
        // range_min() or range_max() and with randomly ordered labels.
        add("1/", 3)
        add("1/0", 15)
        add("1/000", 9)
        add("1/00000", 11)
        add("1/012", 6)
        add("1/01212", 5)
        add("1/312", 17)
        add("1/31200", 4)
        add("1/3120000", 10)
        add("1/333", 20)
        add("1/333333", 18)
        add("5/", 3)
        add("5/3", 31)
        add("5/3333", 27)
        quadraticValidate()
    }

    fun testRandomCellUnions() {
        // Construct cell unions from random S2CellIds at random levels.  Note that
        // because the cell level is chosen uniformly, there is a very high
        // likelihood that the cell unions will overlap.
        repeat(100) { i ->
            add(getRandomCellUnion(), i)
        }
        quadraticValidate()
    }

    fun testContentsIteratorSuppressesDuplicates() {
        // Checks that ContentsIterator stops reporting values once it reaches a
        // node of the cell tree that was visited by the previous call to Begin().
        add("2/1", 1)
        add("2/1", 2)
        add("2/10", 3)
        add("2/100", 4)
        add("2/102", 5)
        add("2/1023", 6)
        add("2/31", 7)
        add("2/313", 8)
        add("2/3132", 9)
        add("3/1", 10)
        add("3/12", 11)
        add("3/13", 12)
        quadraticValidate()

        val contents = S2CellIndex.ContentsIterator(index)
        expectContents("1/123", contents, emptyList())
        expectContents("2/100123", contents, listOf(Pair("2/1", 1), Pair("2/1", 2), Pair("2/10", 3), Pair("2/100", 4)))
        // Check that a second call with the same key yields no additional results.
        expectContents("2/100123", contents, emptyList())
        // Check that seeking to a different branch yields only the new values.
        expectContents("2/10232", contents, listOf(Pair("2/102", 5), Pair("2/1023", 6)))
        // Seek to a node with a different root.
        expectContents("2/313", contents, listOf(Pair("2/31", 7), Pair("2/313", 8)))
        // Seek to a descendant of the previous node.
        expectContents("2/3132333", contents, listOf(Pair("2/3132", 9)))
        // Seek to an ancestor of the previous node.
        expectContents("2/213", contents, emptyList())
        // A few more tests of incremental reporting.
        expectContents("3/1232", contents, listOf(Pair("3/1", 10), Pair("3/12", 11)))
        expectContents("3/133210", contents, listOf(Pair("3/13", 12)))
        expectContents("3/133210", contents, emptyList())
        expectContents("5/0", contents, emptyList())

        // Now try moving backwards, which is expected to yield values that were
        // already reported above.
        expectContents("3/13221", contents, listOf(Pair("3/1", 10), Pair("3/13", 12)))
        expectContents("2/31112", contents, listOf(Pair("2/31", 7)))
    }

    fun testIntersectionOptimization() {
        // Tests various corner cases for the binary search optimization in
        // VisitIntersectingCells.

        add("1/001", 1)
        add("1/333", 2)
        add("2/00", 3)
        add("2/0232", 4)
        build()
        testIntersection(makeCellUnion(listOf("1/010", "1/3")))
        testIntersection(makeCellUnion(listOf("2/010", "2/011", "2/02")))
    }

    fun testIntersectionRandomUnions() {
        // Construct cell unions from random S2CellIds at random levels.  Note that
        // because the cell level is chosen uniformly, there is a very high
        // likelihood that the cell unions will overlap.
        repeat(100) { i ->
            add(getRandomCellUnion(), i)
        }
        build()
        // Now repeatedly query a cell union constructed in the same way.
        repeat(200) {
            testIntersection(getRandomCellUnion())
        }
    }

    fun testIntersectionSemiRandomUnions() {
        // This test also uses random S2CellUnions, but the unions are specially
        // constructed so that interesting cases are more likely to arise.
        repeat(200) {
            S2Random.reset(it)
            index.clear()
            var id = S2CellId.fromDebugString("1/0123012301230123")
            val target = mutableListOf<S2CellId>()
            repeat(100) { i ->
                if (S2Random.oneIn(10)) add(id, i)
                if (S2Random.oneIn(4)) target.add(id)
                if (S2Random.oneIn(2)) id = id.nextWrap()
                if (S2Random.oneIn(6) && !id.isFace()) id = id.parent()
                if (S2Random.oneIn(6) && !id.isLeaf()) id = id.childBegin()
            }
            build()
            testIntersection(S2CellUnion(target))
        }
    }

    fun testIntersection() {
        add(cellStr = "1/012301230123020021", label = 18)
        add(cellStr = "1/012301230123020110", label = 40)
        add(cellStr = "1/012301230123020112", label = 43)
        add(cellStr = "1/01230123012302012", label = 56)
        add(cellStr = "1/0123012301230211", label = 75)
        add(cellStr = "1/0123012301230212", label = 76)
        add(cellStr = "1/012301230123022", label = 85)
        add(cellStr = "1/01230123012310", label = 93)
        build()

        val union = S2CellUnion(
               S2CellId.fromDebugString("1/0123012301230123"),
               S2CellId.fromDebugString("1/012301230123013"),
               S2CellId.fromDebugString("1/0123012301230200"),
               S2CellId.fromDebugString("1/012301230123020100"),
               S2CellId.fromDebugString("1/012301230123020101"),
               S2CellId.fromDebugString("1/012301230123020102"),
               S2CellId.fromDebugString("1/01230123012302011"),
               S2CellId.fromDebugString("1/01230123012302012"),
               S2CellId.fromDebugString("1/0123012301230210"),
               S2CellId.fromDebugString("1/0123012301230211"),
               S2CellId.fromDebugString("1/01230123012310"),
               S2CellId.fromDebugString("1/01230123012311"),
        )

        testIntersection(union)
    }

    // Adds the (cell_id, label) pair to index_ and also contents_ (which is
    // used for independent validation).
    private fun add(cellId: S2CellId, label: Label) {
        index.add(cellId, label)
        contents.add(LabelledCell(cellId, label))
    }

    private fun add(cellStr: String, label: Label) {
        add(S2CellId.fromDebugString(cellStr), label)
    }

    private fun add(cellUnion: S2CellUnion, label: Label) {
        index.add(cellUnion, label)
        for (cell_id in cellUnion) {
            contents.add(LabelledCell(cell_id, label))
        }
    }

    private fun build() {
        index.build()
    }

    // Verifies that the index computes the correct set of (cell_id, label) pairs
    // for every possible leaf cell.  The running time of this function is
    // quadratic in the size of the index.
    private fun quadraticValidate() {
        build()
        verifyCellIterator()
        verifyIndexContents()
        verifyRangeIterators()
    }

    // Verifies that S2CellIndex::CellIterator visits each (cell_id, label) pair
    // exactly once.
    private fun verifyCellIterator() {
        val actual = mutableListOf<LabelledCell>()
        val iterator = CellIterator(index)
        while (!iterator.done()) {
            actual.add(LabelledCell(iterator.cellId(), iterator.label()))
            iterator.next()
        }
        expectEqual(contents, actual)
    }

    private fun verifyRangeIterators() {
        // Test Finish(), which is not otherwise tested below.
        val iterator = S2CellIndex.RangeIterator(index)
        iterator.begin()
        iterator.finish()
        assertTrue(iterator.done())

        // And also for non-empty ranges.
        val nonEmpty = S2CellIndex.NonEmptyRangeIterator(index)
        nonEmpty.begin()
        nonEmpty.finish()
        assertTrue(nonEmpty.done())

        // Iterate through all the ranges in the index.  We simultaneously iterate
        // through the non-empty ranges and check that the correct ranges are found.
        var prevStart = S2CellId.none()
        var nonEmptyPrevStart = S2CellId.none()
        iterator.begin()
        nonEmpty.begin()
        while (!iterator.done()) {
            // Check that seeking in the current range takes us to this range.
            val rangeIterator = S2CellIndex.RangeIterator(index)
            val start = iterator.startId()
            rangeIterator.seek(iterator.startId())
            assertEquals(start, rangeIterator.startId())
            rangeIterator.seek(iterator.limitId().previous())
            assertEquals(start, rangeIterator.startId())

            // And also for non-empty ranges.
            val nonEmpty2 = S2CellIndex.NonEmptyRangeIterator(index)
            val nonEmptyStart = nonEmpty.startId()
            nonEmpty2.seek(iterator.startId())
            assertEquals(nonEmptyStart, nonEmpty2.startId())
            nonEmpty2.seek(iterator.limitId().previous())
            assertEquals(nonEmptyStart, nonEmpty2.startId())

            // Test Prev() and Next().
            if (rangeIterator.prev()) {
                assertEquals(prevStart, rangeIterator.startId())
                rangeIterator.next()
                assertEquals(start, rangeIterator.startId())
            } else {
                assertEquals(start, rangeIterator.startId())
                assertEquals(S2CellId.none(), prevStart)
            }

            // And also for non-empty ranges.
            if (nonEmpty2.prev()) {
                assertEquals(nonEmptyPrevStart, nonEmpty2.startId())
                nonEmpty2.next()
                assertEquals(nonEmptyStart, nonEmpty2.startId())
            } else {
                assertEquals(nonEmptyStart, nonEmpty2.startId())
                assertEquals(S2CellId.none(), nonEmptyPrevStart)
            }

            // Keep the non-empty iterator synchronized with the regular one.
            if (!iterator.isEmpty()) {
                assertEquals(iterator.startId(), nonEmpty.startId())
                assertEquals(iterator.limitId(), nonEmpty.limitId())
                assertFalse(nonEmpty.done())
                nonEmptyPrevStart = nonEmptyStart
                nonEmpty.next()
            }
            prevStart = start
            iterator.next()
        }
        // Verify that the NonEmptyRangeIterator is also finished.
        assertTrue(nonEmpty.done())
    }

    // Verifies that RangeIterator and ContentsIterator can be used to determine
    // the exact set of (s2cell_id, label) pairs that contain any leaf cell.
    private fun verifyIndexContents() {
        // "min_cellid" is the first S2CellId that has not been validated yet.
        var minCellId = S2CellId.begin(S2CellId.kMaxLevel)
        val range = S2CellIndex.RangeIterator(index)
        range.begin()
        while (!range.done()) {
            logger.trace { "VerifyIndexContents: process range (startId = ${range.startId()}, contents = ${range.contents()}, limitId = ${range.limitId()} )" }
            assertEquals(minCellId, range.startId())
            assertTrue(minCellId < range.limitId())
            assertTrue(range.limitId().isLeaf())
            minCellId = range.limitId()

            logger.trace { "minCellId: $minCellId" }

            // Build a list of expected (cell_id, label) pairs for this range.
            val expected = mutableListOf<LabelledCell>()
            for (x in contents) {
                if (x.cellId.rangeMin() <= range.startId() && x.cellId.rangeMax().next() >= range.limitId()) {
                    // The cell contains the entire range.
                    expected.add(x)
                } else {
                    // Verify that the cell does not intersect the range.
                    assertFalse(x.cellId.rangeMin() <= range.limitId().previous() && x.cellId.rangeMax() >= range.startId())
                }
            }
            val actual = mutableListOf<LabelledCell>()
            val contents = ContentsIterator(index)
            contents.startUnion(range)
            while (!contents.done()) {
                actual.add(LabelledCell(contents.cellId(), contents.label()))
                contents.next()
            }
            expectEqual(expected, actual)
            range.next()
        }
        assertEquals(S2CellId.end(S2CellId.kMaxLevel), minCellId)
    }

    // Tests that VisitIntersectingCells() and GetIntersectingLabels() return
    // correct results for the given target.
    private fun testIntersection(target: S2CellUnion) {
        val expected = mutableListOf<LabelledCell>()
        val actual = mutableListOf<LabelledCell>()
        val expectedLabels = TreeSet<Label>()
        val it = CellIterator(index)
        while (!it.done()) {
            if (target.intersects(it.cellId())) {
                expected.add(LabelledCell(it.cellId(), it.label()))
                expectedLabels.add(it.label())
            }
            it.next()
        }
        index.visitIntersectingCells(target, object : S2CellIndex.CellVisitor {
            override fun visit(cellId: S2CellId, label: Label): Boolean {
                actual.add(LabelledCell(cellId, label))
                return true
            }
        })
        expectEqual(expected, actual)
        var actualLabels = index.getIntersectingLabels(target)
        actualLabels = actualLabels.distinct().sorted()
        assertEquals(expectedLabels.toList(), actualLabels)
    }

    // Given an S2CellId "target_str" in human-readable form, expects that the
    // first leaf cell contained by this target will intersect the exact set of
    // (cell_id, label) pairs given by "expected_strs".
    private fun expectContents(targetStr: String, contents: ContentsIterator, expectedStrs: List<Pair<String, Label>>) {
        val range = S2CellIndex.RangeIterator(index)
        range.seek(S2CellId.fromDebugString(targetStr).rangeMin())
        val expected = mutableListOf<LabelledCell>()
        val actual = mutableListOf<LabelledCell>()
        for (p in expectedStrs) {
            expected.add(LabelledCell(S2CellId.fromDebugString(p.first), p.second))
        }
        contents.startUnion(range)
        while (!contents.done()) {
            actual.add(LabelledCell(contents.cellId(), contents.label()))
            contents.next()
        }
        expectEqual(expected, actual)
    }

    // Verifies that "expected" and "actual" have the same contents.  Note that
    // duplicate values are allowed.
    private fun expectEqual(expected: List<LabelledCell>, actual: List<LabelledCell>) {
        assertEquals(expected.sorted(), actual.sorted())
    }

    // Creates a cell union from a small number of random cells at random levels.
    private fun getRandomCellUnion(): S2CellUnion {
        val ids = mutableListOf<S2CellId>()
        repeat(10) {
            ids.add(S2Random.randomCellId())
        }
        return S2CellUnion(ids)
    }

    private fun makeCellUnion(strs: List<String>): S2CellUnion {
        val ids = mutableListOf<S2CellId>()
        for (str in strs) {
            ids.add(S2CellId.fromDebugString(str))
        }
        return S2CellUnion(ids)
    }

}


