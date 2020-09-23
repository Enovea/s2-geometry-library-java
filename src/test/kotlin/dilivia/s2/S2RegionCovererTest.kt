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

import com.google.common.geometry.S2.M_PI
import mu.KotlinLogging
import java.util.*
import kotlin.math.max
import kotlin.math.min

/*
S2_DEFINE_string(max_cells, "4,8",
              "Comma-separated list of values to use for 'max_cells'")

S2_DEFINE_int32(iters, google::DEBUG_MODE ? 1000 : 100000,
             "Number of random caps to try for each max_cells value")
*/

class S2RegionCovererTest : S2GeometryTestCase() {

    private val logger = KotlinLogging.logger { }

    companion object {
        val kMaxCells = listOf(4, 8)
        val kIters = 100000
    }

    fun testRandomCells() {
        val coverer = S2RegionCoverer(maxCells = 1)

        // Test random cell ids at all levels.
        repeat(10000) { i ->
            val id = S2Random.randomCellId()
            logger.trace { "Iteration $i, cell ID token ${id.toToken()}" }
            val covering = coverer.getCovering(S2Cell(id))
            assertEquals(1, covering.numCells())
            assertEquals(id, covering[0])
        }
    }

    fun checkCovering(options: S2RegionCoverer, region: S2Region, covering: List<S2CellId>, interior: Boolean) {
        // Keep track of how many cells have the same options.min_level() ancestor.
        val min_level_cells = mutableMapOf<S2CellId, Int>()
        for (cell_id in covering) {
            val level = cell_id.level()
            assertTrue(level >= options.getMinLevel())
            assertTrue(level <= options.getMaxLevel())
            assertEquals((level - options.getMinLevel()) % options.getLevelMod(), 0)
            min_level_cells.compute(cell_id.parent(options.getMinLevel())) { _, count -> if (count == null) 1 else count + 1 }
        }
        if (covering.size > options.getMaxCells()) {
            // If the covering has more than the requested number of cells, then check
            // that the cell count cannot be reduced by using the parent of some cell.
            min_level_cells.forEach { (_, count) -> assertEquals(1, count) }
        }
        if (interior) {
            covering.forEach { cellId ->
                assertTrue(region.contains(S2Cell(cellId)))
            }
        } else {
            val cellUnion = S2CellUnion(covering)
            checkCovering(region, cellUnion, true, S2CellId())
        }
    }

    fun testRandomCaps() {
        logger.info { "Test random caps" }
        val kMaxLevel = S2CellId.kMaxLevel
        var minLevel = 0
        var maxLevel = 0
        var maxCells = 0
        var levelMod = 1
        repeat(1000) { i ->
            if (i % 100 == 0 || logger.isTraceEnabled) logger.info { "Iteration $i" }
            maxLevel = S2Random.randomInt(kMaxLevel + 1)
            minLevel =if (maxLevel == 0) 0 else S2Random.randomInt(maxLevel)
            maxCells = S2Random.skewed(10)
            levelMod = S2Random.randomInt(1, 4)
            val max_area = min(4 * M_PI, (3 * maxCells + 1) * S2Cell.averageArea(minLevel))
            val cap = S2Random.randomCap(0.1 * S2Cell.averageArea(kMaxLevel), max_area)
            val coverer = S2RegionCoverer(maxCells, minLevel, maxLevel, levelMod)
            val covering = mutableListOf<S2CellId>()
            val interior = mutableListOf<S2CellId>()

            logger.trace { """
                | Test random cap:
                | ============================================
                | maxCells: $maxCells
                | minLevel: $minLevel
                | maxLevel: $maxLevel
                | levelMod: $levelMod
                | cap: $cap
                | -------------------------------------------
            """.trimMargin() }

            coverer.getCovering(cap, covering)
            logger.trace { "Covering = $covering" }
            checkCovering(coverer, cap, covering, false)
            coverer.getInteriorCovering(cap, interior)
            logger.trace { "Interior covering: $interior" }
            checkCovering(coverer, cap, interior, true)

            // Check that GetCovering is deterministic.
            val covering2 = mutableListOf<S2CellId>()
            coverer.getCovering(cap, covering2)
            assertEquals(covering, covering2)

            // Also check S2CellUnion::Denormalize().  The denormalized covering
            // may still be different and smaller than "covering" because
            // S2RegionCoverer does not guarantee that it will not output all four
            // children of the same parent.
            val cells = S2CellUnion(covering)
            val denormalized = cells.denormalize(minLevel, levelMod)
            logger.trace { "Denormalized: $denormalized" }
            checkCovering(coverer, cap, denormalized, false)
        }
    }

    fun testSimpleCoverings() {
        val kMaxLevel = S2CellId.kMaxLevel
        val maxCells = Int.MAX_VALUE
        repeat(1000) { i ->
            val level = S2Random.randomInt(kMaxLevel + 1)
            val maxArea = min(4 * M_PI, 1000 * S2Cell.averageArea(level))
            val cap = S2Random.randomCap(0.1 * S2Cell.averageArea(kMaxLevel), maxArea)
            val covering = mutableListOf<S2CellId>()
            S2RegionCoverer.getSimpleCovering(cap, cap.center, level, covering)
            checkCovering(S2RegionCoverer(maxCells, level, level), cap, covering, false)
        }
    }

    // We keep a priority queue of the caps that had the worst approximation
    // ratios so that we can print them at the end.
    data class WorstCap(val ratio: Double, val cap: S2Cap, val num_cells: Int) : Comparable<WorstCap> {
        override fun compareTo(other: WorstCap): Int = -ratio.compareTo(other.ratio)
    }

    fun testAccuracy(max_cells: Int) {
        logger.info { "testAccuracy: $max_cells cells" }

        val kNumMethods = 1
        // This code is designed to evaluate several approximation algorithms and
        // figure out which one works better.  The way to do this is to hack the
        // S2RegionCoverer interface to add a global variable to control which
        // algorithm (or variant of an algorithm) is selected, and then assign to
        // this variable in the "method" loop below.  The code below will then
        // collect statistics on all methods, including how often each one wins in
        // terms of cell count and approximation area.

        val coverer = S2RegionCoverer()
        coverer.setMaxCells(max_cells)

        val ratio_total = Array<Double>(kNumMethods) { 0.0 }
        val min_ratio = Array<Double>(kNumMethods) { 1e20 }  // initialized in loop below
        val max_ratio = Array<Double>(kNumMethods) { 0.0 }
        val ratios = Array<MutableList<Double>>(kNumMethods) { mutableListOf() }
        val cell_total = IntArray(kNumMethods)
        val area_winner_tally = IntArray(kNumMethods)
        val cell_winner_tally = IntArray(kNumMethods)
        val kMaxWorstCaps = 10
        val worst_caps = Array<PriorityQueue<WorstCap>>(kNumMethods) { PriorityQueue() }

        repeat(kIters) { i ->
            // Choose the log of the cap area to be uniformly distributed over
            // the allowable range.  Don't try to approximate regions that are so
            // small they can't use the given maximum number of cells efficiently.
            val min_cap_area = S2Cell.averageArea(S2CellId.kMaxLevel) * max_cells * max_cells
            // Coverings for huge caps are not interesting, so limit the max area too.
            val cap = S2Random.randomCap(min_cap_area, 0.1 * M_PI)
            val cap_area = cap.area

            var min_area = 1e30
            var min_cells = 1 shl 30
            val area = DoubleArray(kNumMethods)
            val cells = IntArray(kNumMethods)
            for (method in 0 until kNumMethods) {
                // If you want to play with different methods, do this:
                // S2RegionCoverer::method_number = method
                val covering = mutableListOf<S2CellId>()
                coverer.getCovering(cap, covering)

                var union_area = 0.0
                for (cell_id in covering) {
                    union_area += S2Cell(cell_id).exactArea()
                }
                cells[method] = covering.size
                min_cells = min(cells[method], min_cells)
                area[method] = union_area
                min_area = min(area[method], min_area)
                cell_total[method] += cells[method]
                val ratio = area[method] / cap_area
                ratio_total[method] += ratio
                min_ratio[method] = min(ratio, min_ratio[method])
                max_ratio[method] = max(ratio, max_ratio[method])
                ratios[method].add(ratio)
                if (worst_caps[method].size < kMaxWorstCaps) {
                    worst_caps[method].add(WorstCap(ratio, cap, cells[method]))
                } else if (ratio > worst_caps[method].peek().ratio) {
                    worst_caps[method].poll()
                    worst_caps[method].add(WorstCap(ratio, cap, cells[method]))
                }
            }
            for (method in 0 until kNumMethods) {
                if (area[method] == min_area) ++area_winner_tally[method]
                if (cells[method] == min_cells) ++cell_winner_tally[method]
            }
        }
        for (method in 0 until kNumMethods) {
            println("\nMax cells %d, method %d:".format(max_cells, method))
            println("  Average cells: %.4f".format(cell_total[method] / kIters.toDouble()))
            println("  Average area ratio: %.4f".format(ratio_total[method] / kIters.toDouble()))
            val mratios = ratios[method]
            mratios.sort()
            println("  Median ratio: %.4f".format(mratios[mratios.size / 2]))
            println("  Max ratio: %.4f".format(max_ratio[method]))
            println("  Min ratio: %.4f".format(min_ratio[method]))
            if (kNumMethods > 1) {
                println("  Cell winner probability: %.4f".format(cell_winner_tally[method] / kIters.toDouble()))
                println("  Area winner probability: %.4f".format(area_winner_tally[method] / kIters.toDouble()))
            }
            println("  Caps with the worst approximation ratios:")
            while (worst_caps[method].isNotEmpty()) {
                val w = worst_caps[method].poll()
                val ll = S2LatLng.fromPoint(w.cap.center)
                println("    Ratio %.4f, Cells %d, Center (%.8f, %.8f), Km %.6f".format(
                        w.ratio, w.num_cells,
                        ll.lat().degrees(), ll.lng().degrees(),
                        w.cap.radius().radians * 6367.0))
            }
        }
    }

    fun testAccuracy() {
        kMaxCells.forEach { maxCells -> testAccuracy(maxCells) }
    }

    fun testInteriorCovering() {
        logger.info { "Test interior covering" }
        // We construct the region the following way. Start with S2 cell of level l.
        // Remove from it one of its grandchildren (level l+2). If we then set
        //   min_level < l + 1
        //   max_level > l + 2
        //   max_cells = 3
        // the best interior covering should contain 3 children of the initial cell,
        // that were not effected by removal of a grandchild.
        val level = 12
        val small_cell = S2CellId.fromPoint(S2Random.randomPoint()).parent(level + 2)
        val large_cell = small_cell.parent(level)
        val diff = S2CellUnion(large_cell).difference(S2CellUnion(small_cell))
        val coverer = S2RegionCoverer(maxCells = 3, minLevel = level, maxLevel = level + 3)
        val interior = mutableListOf<S2CellId>()
        coverer.getInteriorCovering(diff, interior)
        logger.info { "Small cell = $small_cell" }
        logger.info { "Large cell = $large_cell" }
        logger.info { "Diff = $diff" }
        logger.info { "Interior = $interior" }
        assertEquals(interior.size, 3)
        for (i in 0..2) {
            assertEquals(interior[i].level(), level + 1)
        }
    }

    fun testGetFastCoveringHugeFixedLevelCovering() {
        // Test a "fast covering" with a huge number of cells due to min_level().
        val coverer = S2RegionCoverer(minLevel = 10)
        val covering = mutableListOf<S2CellId>()
        val region = S2Cell(S2CellId.fromDebugString("1/23"))
        coverer.getFastCovering(region, covering)
        assertTrue(covering.size >= 1 shl 16)
    }

    fun isCanonical(input_str: List<String>, coverer: S2RegionCoverer): Boolean =
            coverer.isCanonical(input_str.map { S2CellId.fromDebugString(it) })

    fun testIsCanonicalInvalidS2CellId() {
        assertTrue(isCanonical(listOf("1/"), S2RegionCoverer()))
        assertFalse(isCanonical(listOf("invalid"), S2RegionCoverer()))
    }

    fun testIsCanonicalUnsorted() {
        assertTrue(isCanonical(listOf("1/1", "1/3"), S2RegionCoverer()))
        assertFalse(isCanonical(listOf("1/3", "1/1"), S2RegionCoverer()))
    }

    fun testIsCanonicalOverlapping() {
        assertTrue(isCanonical(listOf("1/2", "1/33"), S2RegionCoverer()))
        assertFalse(isCanonical(listOf("1/3", "1/33"), S2RegionCoverer()))
    }

    fun testIsCanonicalMinLevel() {
        val options = S2RegionCoverer(minLevel = 2)
        assertTrue(isCanonical(listOf("1/31"), options))
        assertFalse(isCanonical(listOf("1/3"), options))
    }

    fun testIsCanonicalMaxLevel() {
        val options = S2RegionCoverer(maxLevel = 2)
        assertTrue(isCanonical(listOf("1/31"), options))
        assertFalse(isCanonical(listOf("1/312"), options))
    }

    fun testIsCanonicalLevelMod() {
        val options = S2RegionCoverer(levelMod = 2)
        assertTrue(isCanonical(listOf("1/31"), options))
        assertFalse(isCanonical(listOf("1/312"), options))
    }

    fun testIsCanonicalMaxCells() {
        val options = S2RegionCoverer(maxCells = 2)
        assertTrue(isCanonical(listOf("1/1", "1/3"), options))
        assertFalse(isCanonical(listOf("1/1", "1/3", "2/"), options))
        assertTrue(isCanonical(listOf("1/123", "2/1", "3/0122"), options))
    }

    fun testIsCanonicalNormalized() {
        // Test that no sequence of cells could be replaced by an ancestor.
        val options = S2RegionCoverer()
        assertTrue(isCanonical(listOf("1/01", "1/02", "1/03", "1/10", "1/11"), options))
        assertFalse(isCanonical(listOf("1/00", "1/01", "1/02", "1/03", "1/10"), options))

        assertTrue(isCanonical(listOf("0/22", "1/01", "1/02", "1/03", "1/10"), options))
        assertFalse(isCanonical(listOf("0/22", "1/00", "1/01", "1/02", "1/03"), options))

        options.setMaxCells(20)
        options.setLevelMod(2)
        assertTrue(isCanonical(listOf(
                "1/1101", "1/1102", "1/1103", "1/1110",
                "1/1111", "1/1112", "1/1113", "1/1120",
                "1/1121", "1/1122", "1/1123", "1/1130",
                "1/1131", "1/1132", "1/1133", "1/1200"
        ), options))
        assertFalse(isCanonical(listOf(
                "1/1100", "1/1101", "1/1102", "1/1103",
                "1/1110", "1/1111", "1/1112", "1/1113",
                "1/1120", "1/1121", "1/1122", "1/1123",
                "1/1130", "1/1131", "1/1132", "1/1133"
        ), options))
    }

    fun testCanonicalizeCovering(input_str: List<String>, expected_str: List<String>, coverer: S2RegionCoverer) {
        val actual = input_str.map { str -> S2CellId.fromDebugString(str) }.toMutableList()
        val expected = expected_str.map { str -> S2CellId.fromDebugString(str) }
        val isCanonical = coverer.isCanonical(actual)
        coverer.canonicalizeCovering(actual)
        val isCanocicalAfterCanolicalization = coverer.isCanonical(actual)

        logger.trace { """
            | 
            | testCanonicalizeCovering()
            | ===========================================
            | input: $input_str
            | expected: $expected_str
            | ------------------------------------------
            | isCanonical: $isCanonical
            | canonicalized: ${actual.joinToString(", ") { it.toString() }}
            | isCanonicalAfter: $isCanocicalAfterCanolicalization
            | """.trimMargin() }
        assertFalse(isCanonical)
        assertTrue(isCanocicalAfterCanolicalization)
        assertEquals(expected, actual)
    }

    fun testCanonicalizeCoveringUnsortedDuplicateCells() {
        val coverer = S2RegionCoverer()
        testCanonicalizeCovering(
                listOf("1/200", "1/13122", "1/20", "1/131", "1/13100"),
                listOf("1/131", "1/20"),
                coverer)
    }

    fun testCanonicalizeCoveringMaxLevelExceeded() {
        val coverer = S2RegionCoverer()
        coverer.setMaxLevel(2)
        testCanonicalizeCovering(listOf("0/3001", "0/3002", "4/012301230123"), listOf("0/30", "4/01"), coverer)
    }

    fun testCanonicalizeCoveringWrongLevelMod() {
        val coverer = S2RegionCoverer()
        coverer.setMinLevel(1)
        coverer.setLevelMod(3)
         testCanonicalizeCovering(listOf("0/0", "1/11", "2/222", "3/3333"), listOf("0/0", "1/1", "2/2", "3/3333"), coverer)
    }

    fun testCanonicalizeCoveringReplacedByParent() {
        // Test that 16 children are replaced by their parent when level_mod == 2.
        val coverer = S2RegionCoverer()
        coverer.setLevelMod(2)
        testCanonicalizeCovering(listOf(
                "0/00", "0/01", "0/02", "0/03", "0/10", "0/11", "0/12", "0/13",
                "0/20", "0/21", "0/22", "0/23", "0/30", "0/31", "0/32", "0/33"
        ), listOf("0/"), coverer)
    }

    fun testCanonicalizeCoveringDenormalizedCellUnion() {
        // Test that all 4 children of a cell may be used when this is necessary to
        // satisfy min_level() or level_mod()
        val coverer = S2RegionCoverer()
        coverer.setMinLevel(1)
        coverer.setLevelMod(2)
        testCanonicalizeCovering(
                listOf("0/", "1/130", "1/131", "1/132", "1/133"),
                listOf("0/0", "0/1", "0/2", "0/3", "1/130", "1/131", "1/132", "1/133"),
                coverer)
    }

    fun testCanonicalizeCoveringMaxCellsMergesSmallest() {
        // When there are too many cells, the smallest cells should be merged first.
        val coverer = S2RegionCoverer()
        coverer.setMaxCells(3)
        testCanonicalizeCovering(listOf("0/", "1/0", "1/1", "2/01300", "2/0131313"), listOf("0/", "1/", "2/013"), coverer)
    }

    fun testCanonicalizeCoveringMaxCellsMergesRepeatedly() {
        // Check that when merging creates a cell when all 4 children are present,
        // those cells are merged into their parent (repeatedly if necessary).
        val coverer = S2RegionCoverer()
        coverer.setMaxCells(8)
        testCanonicalizeCovering(listOf(
                "0/0121", "0/0123", "1/0", "1/1", "1/2", "1/30", "1/32", "1/33",
                "1/311", "1/312", "1/313", "1/3100", "1/3101", "1/3103",
                "1/31021", "1/31023"
        ), listOf("0/0121", "0/0123", "1/"), coverer)

        // 0/0121, 0/0123, 1/0, 1/1, 1/2, 1/30, 1/31, 1/32, 1/33,
    }

    fun testJavaCcConsistencyCheckCovering() {
        val points = listOf(
            S2LatLng.fromDegrees(-33.8663457, 151.1960891).toPoint(),
            S2LatLng.fromDegrees(-33.866094000000004, 151.19517439999998).toPoint()
        )

        val coverer = S2RegionCoverer()
                coverer.setMinLevel(0)
        coverer.setMaxLevel(22)
        coverer.setMaxCells(Int.MAX_VALUE)
        val ret = coverer.getCovering(S2Polyline(points))
        val expected = listOf(
                    "6b12ae36313d", "6b12ae36313f", "6b12ae363141", "6b12ae363143",
                    "6b12ae363145", "6b12ae363159", "6b12ae36315b", "6b12ae363343",
                    "6b12ae363345", "6b12ae36334d", "6b12ae36334f", "6b12ae363369",
                    "6b12ae36336f", "6b12ae363371", "6b12ae363377", "6b12ae363391",
                    "6b12ae363393", "6b12ae36339b", "6b12ae36339d", "6b12ae3633e3",
                    "6b12ae3633e5", "6b12ae3633ed", "6b12ae3633ef", "6b12ae37cc11",
                    "6b12ae37cc13", "6b12ae37cc1b", "6b12ae37cc1d", "6b12ae37cc63",
                    "6b12ae37cc65", "6b12ae37cc6d", "6b12ae37cc6f", "6b12ae37cc89",
                    "6b12ae37cc8f", "6b12ae37cc91", "6b12ae37cc97", "6b12ae37ccb1",
                    "6b12ae37ccb3", "6b12ae37ccbb", "6b12ae37ccbd", "6b12ae37cea5",
                    "6b12ae37cea7", "6b12ae37cebb"
                )
        val actual = ret.map { it.toToken() }
        assertEquals(expected, actual)
    }



    fun notAtestCap() {
        val coverer = S2RegionCoverer(maxCells = 45, minLevel = 8, maxLevel = 12, levelMod = 1)
        val cap = S2Cap(S2Point(-0.9652843302739064, -0.190115886479274, -0.17911479960030016), S1ChordAngle.fromLength2(3.1560079238412953E-4))
        val covering = mutableListOf<S2CellId>()
        coverer.getCovering(cap, covering)
        val union = S2CellUnion(covering)
        assertEquals(cap, union.capBound)

        assertEquals(listOf(
                "3/200222201", "3/200222202", "3/200222212", "3/200222213", "3/20022222", "3/20022223",
                "3/20022230", "3/201332323", "3/20133233", "3/20133300", "3/20133301", "3/20133302", "3/20133303",
                "3/20133310", "3/201333123", "3/20133313", "3/20133320", "3/201333210", "3/201333213", "3/20133322",
                "3/20133323", "3/20133330", "3/20133331", "3/20133332", "3/20133333", "3/20200000", "3/20200001",
                "3/20200002", "3/20200003", "3/20200010", "3/20200013", "3/20200020", "3/2020002322", "3/20200030",
                "3/20200031", "3/20200032", "3/20200033", "3/20200100", "3/20200101000", "3/2031110333", "3/20311110",
                "3/20311111", "3/203111120", "3/203111121", "3/20311113111").map { S2CellId.fromDebugString(it) },
                covering)
    }
}