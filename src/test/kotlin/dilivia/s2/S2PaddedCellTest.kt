package dilivia.s2

import dilivia.s2.math.R2Point
import dilivia.s2.region.S2Cell
import mu.KotlinLogging
import kotlin.math.pow

class S2PaddedCellTest : S2GeometryTestCase() {

    private val logger = KotlinLogging.logger {  }

    fun compareS2CellToPadded(cell: S2Cell, pcell: S2PaddedCell, padding: Double) {
        logger.trace { "Compare S2Cell $cell with padded (padding=$padding) $pcell" }
        assertEquals(cell.id(), pcell.id)
        assertEquals(cell.level(), pcell.level)
        assertEquals(padding, pcell.padding)
        logger.trace { "Cell bound          = ${cell.boundUV()}, ${S2PaddedCell(cell.id(), 0.0).bound}" }
        logger.trace { "Expanded cell bound = ${cell.boundUV().expanded(padding)}" }
        logger.trace { "Padded cell bound   = ${pcell.bound}" }
        assertEquals(cell.boundUV().expanded(padding), pcell.bound)
        val centerUv = cell.id().getCenterUV()
        assertEquals(R2Rect.fromPoint(centerUv).expanded(padding), pcell.middle())
        assertEquals(cell.getCenter(), pcell.getCenter())
    }

    fun testS2CellMethods() {
        // Test the S2PaddedCell methods that have approximate S2Cell equivalents.
        val kIters = 1000
        repeat(kIters) { iter ->
            logger.trace { "iteration ${iter + 1}" }
            val id = S2Random.randomCellId()
            val padding = 1e-15.pow(S2Random.randomDouble())
            val cell = S2Cell(id)
            val pcell = S2PaddedCell(id, padding)
            compareS2CellToPadded(cell, pcell, padding)
            if (!id.isLeaf()) {
                val children = Array(4) { S2Cell() }
                assertTrue(cell.subdivide(children))
                for (pos in 0..3) {
                  val (i, j) = pcell.getChildIJ(pos)
                    compareS2CellToPadded(children[pos], S2PaddedCell(pcell, i, j), padding)
                }
            }
        }
    }

    fun testGetEntryExitVertices() {
        val kIters = 1000
        repeat(kIters) {
            val id = S2Random.randomCellId()
            // Check that entry/exit vertices do not depend on padding.
            assertEquals(S2PaddedCell(id, 0.0).getEntryVertex(), S2PaddedCell(id, 0.5).getEntryVertex())
            assertEquals(S2PaddedCell(id, 0.0).getExitVertex(), S2PaddedCell(id, 0.5).getExitVertex())

            // Check that the exit vertex of one cell is the same as the entry vertex
            // of the immediately following cell.  (This also tests wrapping from the
            // end to the start of the S2CellId curve with high probability.)
            assertEquals(S2PaddedCell(id, 0.0).getExitVertex(), S2PaddedCell(id.nextWrap(), 0.0).getEntryVertex())

            // Check that the entry vertex of a cell is the same as the entry vertex
            // of its first child, and similarly for the exit vertex.
            if (!id.isLeaf()) {
                assertEquals(S2PaddedCell(id, 0.0).getEntryVertex(), S2PaddedCell(id.child(0), 0.0).getEntryVertex())
                assertEquals(S2PaddedCell(id, 0.0).getExitVertex(), S2PaddedCell(id.child(3), 0.0).getExitVertex())
            }
        }
    }

    fun sampleInterval(x: R1Interval): Double {
        return S2Random.randomDouble(x.lo, x.hi)
    }

    fun testShrinkToFit() {
        val kIters = 1000
        repeat(kIters) {
            // Start with the desired result and work backwards.
            val result = S2Random.randomCellId()
            val resultUv = result.getBoundUV()
            val sizeUv = resultUv.size

            // Find the biggest rectangle that fits in "result" after padding.
            // (These calculations ignore numerical errors.)
            val maxPadding = 0.5 * kotlin.math.min(sizeUv[0], sizeUv[1])
            val padding = maxPadding * S2Random.randomDouble()
            val maxRect = resultUv.expanded(-padding)

            // Start with a random subset of the maximum rectangle.
            val a = R2Point(sampleInterval(maxRect[0]), sampleInterval(maxRect[1])).toMutable()
            val b = R2Point(sampleInterval(maxRect[0]), sampleInterval(maxRect[1])).toMutable()
            if (!result.isLeaf()) {
                // If the result is not a leaf cell, we must ensure that no child of
                // "result" also satisfies the conditions of ShrinkToFit().  We do this
                // by ensuring that "rect" intersects at least two children of "result"
                // (after padding).
                val axis = S2Random.randomInt(2)
                val center = result.getCenterUV()[axis]

                // Find the range of coordinates that are shared between child cells
                // along that axis.
                val shared = R1Interval(center - padding, center + padding)
                val mid = sampleInterval(shared.intersection(maxRect[axis]))
                a[axis] = sampleInterval(R1Interval(maxRect[axis].lo, mid))
                b[axis] = sampleInterval(R1Interval(mid, maxRect[axis].hi))
            }
            val rect = R2Rect.fromPointPair(a, b)

            // Choose an arbitrary ancestor as the S2PaddedCell.
            val initialId = result.parent(S2Random.randomInt(result.level() + 1))
            assertEquals(result, S2PaddedCell(initialId, padding).shrinkToFit(rect))
        }
    }

}
