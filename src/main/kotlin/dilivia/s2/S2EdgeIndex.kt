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

import com.google.common.base.Preconditions
import com.google.common.collect.Lists
import com.google.common.collect.Sets
import dilivia.s2.S2CellId.Companion.begin
import dilivia.s2.S2CellId.Companion.end
import dilivia.s2.S2CellId.Companion.fromPoint
import dilivia.s2.S2CellId.Companion.sentinel
import dilivia.s2.S2Point.Companion.crossProd
import dilivia.s2.S2Point.Companion.div
import dilivia.s2.S2Point.Companion.minus
import dilivia.s2.S2Point.Companion.normalize
import dilivia.s2.S2Point.Companion.plus
import dilivia.s2.S2Point.Companion.times
import java.util.*

@Strictfp
abstract class S2EdgeIndex {
    /**
     * The cell containing each edge, as given in the parallel array
     * `edges`.
     */
    private var cells: LongArray? = null

    /**
     * The edge contained by each cell, as given in the parallel array
     * `cells`.
     */
    private var edges: IntArray? = null

    /**
     * No cell strictly below this level appears in mapping. Initially leaf level,
     * that's the minimum level at which we will ever look for test edges.
     */
    private var minimumS2LevelUsed = 0

    /**
     * Has the index been computed already?
     */
    var isIndexComputed = false
        private set

    /**
     * Number of queries so far
     */
    private var queryCount = 0

    /**
     * Empties the index in case it already contained something.
     */
    fun reset() {
        minimumS2LevelUsed = S2CellId.kMaxLevel
        isIndexComputed = false
        queryCount = 0
        cells = null
        edges = null
    }

    /** Computes the index (if it has not been previously done).  */
    fun computeIndex() {
        if (isIndexComputed) {
            return
        }
        val cellList: MutableList<Long> = Lists.newArrayList()
        val edgeList: MutableList<Int> = Lists.newArrayList()
        for (i in 0 until numEdges) {
            val from = edgeFrom(i)
            val to = edgeTo(i)
            val cover = Lists.newArrayList<S2CellId>()
            val level = getCovering(from, to, true, cover)
            minimumS2LevelUsed = Math.min(minimumS2LevelUsed, level)
            for (cellId in cover) {
                cellList.add(cellId.id.toLong())
                edgeList.add(i)
            }
        }
        cells = LongArray(cellList.size)
        edges = IntArray(edgeList.size)
        for (i in cells!!.indices) {
            cells!![i] = cellList[i]
            edges!![i] = edgeList[i]
        }
        sortIndex()
        isIndexComputed = true
    }

    /** Sorts the parallel `cells` and `edges` arrays.  */
    private fun sortIndex() {
        // create an array of indices and sort based on the values in the parallel
        // arrays at each index
        val indices = arrayOfNulls<Int>(cells!!.size)
        for (i in indices.indices) {
            indices[i] = i
        }
        Arrays.sort(indices) { index1, index2 -> compare(cells!![index1!!], edges!![index1], cells!![index2!!], edges!![index2]) }
        // copy the cells and edges in the order given by the sorted list of indices
        val newCells = LongArray(cells!!.size)
        val newEdges = IntArray(edges!!.size)
        for (i in indices.indices) {
            newCells[i] = cells!![indices[i]!!]
            newEdges[i] = edges!![indices[i]!!]
        }
        // replace the cells and edges with the sorted arrays
        cells = newCells
        edges = newEdges
    }

    /**
     * Tell the index that we just received a new request for candidates. Useful
     * to compute when to switch to quad tree.
     */
    protected fun incrementQueryCount() {
        ++queryCount
    }

    /**
     * If the index hasn't been computed yet, looks at how much work has gone into
     * iterating using the brute force method, and how much more work is planned
     * as defined by 'cost'. If it were to have been cheaper to use a quad tree
     * from the beginning, then compute it now. This guarantees that we will never
     * use more than twice the time we would have used had we known in advance
     * exactly how many edges we would have wanted to test. It is the theoretical
     * best.
     *
     * The value 'n' is the number of iterators we expect to request from this
     * edge index.
     *
     * If we have m data edges and n query edges, then the brute force cost is m
     * * n * testCost where testCost is taken to be the cost of
     * EdgeCrosser.robustCrossing, measured to be about 30ns at the time of this
     * writing.
     *
     * If we compute the index, the cost becomes: m * costInsert + n *
     * costFind(m)
     *
     * - costInsert can be expected to be reasonably stable, and was measured at
     * 1200ns with the BM_QuadEdgeInsertionCost benchmark.
     *
     * - costFind depends on the length of the edge . For m=1000 edges, we got
     * timings ranging from 1ms (edge the length of the polygon) to 40ms. The
     * latter is for very long query edges, and needs to be optimized. We will
     * assume for the rest of the discussion that costFind is roughly 3ms.
     *
     * When doing one additional query, the differential cost is m * testCost -
     * costFind(m) With the numbers above, it is better to use the quad tree (if
     * we have it) if m >= 100.
     *
     * If m = 100, 30 queries will give m*n*testCost = m * costInsert = 100ms,
     * while the marginal cost to find is 3ms. Thus, this is a reasonable thing to
     * do.
     */
    fun predictAdditionalCalls(n: Int) {
        if (isIndexComputed) {
            return
        }
        if (numEdges > 100 && queryCount + n > 30) {
            computeIndex()
        }
    }

    /**
     * Overwrite these functions to give access to the underlying data. The
     * function getNumEdges() returns the number of edges in the index, while
     * edgeFrom(index) and edgeTo(index) return the "from" and "to" endpoints of
     * the edge at the given index.
     */
    protected abstract val numEdges: Int
    protected abstract fun edgeFrom(index: Int): S2Point
    protected abstract fun edgeTo(index: Int): S2Point

    /**
     * Appends to "candidateCrossings" all edge references which may cross the
     * given edge. This is done by covering the edge and then finding all
     * references of edges whose coverings overlap this covering. Parent cells are
     * checked level by level. Child cells are checked all at once by taking
     * advantage of the natural ordering of S2CellIds.
     */
    protected fun findCandidateCrossings(a: S2Point, b: S2Point, candidateCrossings: MutableList<Int>) {
        Preconditions.checkState(isIndexComputed)
        val cover = Lists.newArrayList<S2CellId>()
        getCovering(a, b, false, cover)

        // Edge references are inserted into the map once for each covering cell, so
        // absorb duplicates here
        val uniqueSet: MutableSet<Int> = HashSet()
        getEdgesInParentCells(cover, uniqueSet)

        // TODO(user): An important optimization for long query
        // edges (Contains queries): keep a bounding cap and clip the query
        // edge to the cap before starting the descent.
        getEdgesInChildrenCells(a, b, cover, uniqueSet)
        candidateCrossings.clear()
        candidateCrossings.addAll(uniqueSet)
    }

    /**
     * Computes a cell covering of an edge. Clears edgeCovering and returns the
     * level of the s2 cells used in the covering (only one level is ever used for
     * each call).
     *
     * If thickenEdge is true, the edge is thickened and extended by 1% of its
     * length.
     *
     * It is guaranteed that no child of a covering cell will fully contain the
     * covered edge.
     */
    private fun getCovering(
            a: S2Point, b: S2Point, thickenEdge: Boolean, edgeCovering: ArrayList<S2CellId>): Int {
        edgeCovering.clear()

        // Selects the ideal s2 level at which to cover the edge, this will be the
        // level whose S2 cells have a width roughly commensurate to the length of
        // the edge. We multiply the edge length by 2*THICKENING to guarantee the
        // thickening is honored (it's not a big deal if we honor it when we don't
        // request it) when doing the covering-by-cap trick.
        val edgeLength = a.angle(b)
        val idealLevel = S2CellMetrics.kMinWidth.getLevelForMinValue(edgeLength * (1 + 2 * THICKENING))
        val containingCellId: S2CellId
        containingCellId = if (!thickenEdge) {
            containingCell(a, b)
        } else {
            if (idealLevel == S2CellId.kMaxLevel) {
                // If the edge is tiny, instabilities are more likely, so we
                // want to limit the number of operations.
                // We pretend we are in a cell much larger so as to trigger the
                // 'needs covering' case, so we won't try to thicken the edge.
                S2CellId(0xFFF0.toULong()).parent(3)
            } else {
                val pq = times(minus(b, a), THICKENING)
                val ortho = times(normalize(crossProd(pq, a)), edgeLength * THICKENING)
                val p = minus(a, pq)
                val q = plus(b, pq)
                // If p and q were antipodal, the edge wouldn't be lengthened,
                // and it could even flip! This is not a problem because
                // idealLevel != 0 here. The farther p and q can be is roughly
                // a quarter Earth away from each other, so we remain
                // Theta(THICKENING).
                containingCell(minus(p, ortho), plus(p, ortho), minus(q, ortho),
                        plus(q, ortho))
            }
        }

        // Best case: edge is fully contained in a cell that's not too big.
        if (!containingCellId.equals(sentinel())
                && containingCellId.level() >= idealLevel - 2) {
            edgeCovering.add(containingCellId)
            return containingCellId.level()
        }
        if (idealLevel == 0) {
            // Edge is very long, maybe even longer than a face width, so the
            // trick below doesn't work. For now, we will add the whole S2 sphere.
            // TODO(user): Do something a tad smarter (and beware of the
            // antipodal case).
            var cellid = begin(0)
            while (!cellid.equals(end(0))) {
                edgeCovering.add(cellid)
                cellid = cellid.next()
            }
            return 0
        }
        // TODO(user): Check trick below works even when vertex is at
        // interface
        // between three faces.

        // Use trick as in S2PolygonBuilder.PointIndex.findNearbyPoint:
        // Cover the edge by a cap centered at the edge midpoint, then cover
        // the cap by four big-enough cells around the cell vertex closest to the
        // cap center.
        val middle = normalize(div(plus(a, b), 2.0))
        val actualLevel = Math.min(idealLevel, S2CellId.kMaxLevel - 1)
        fromPoint(middle).appendVertexNeighbors(actualLevel, edgeCovering)
        return actualLevel
    }

    /**
     * Filters a list of entries down to the inclusive range defined by the given
     * cells, in `O(log N)` time.
     *
     * @param cell1 One side of the inclusive query range.
     * @param cell2 The other side of the inclusive query range.
     * @return An array of length 2, containing the start/end indices.
     */
    private fun getEdges(cell1: Long, cell2: Long): IntArray {
        // ensure cell1 <= cell2
        var cell1 = cell1
        var cell2 = cell2
        if (cell1 > cell2) {
            val temp = cell1
            cell1 = cell2
            cell2 = temp
        }
        // The binary search returns -N-1 to indicate an insertion point at index N,
        // if an exact match cannot be found. Since the edge indices queried for are
        // not valid edge indices, we will always get -N-1, so we immediately
        // convert to N.
        return intArrayOf(
                -1 - binarySearch(cell1, Int.MIN_VALUE),
                -1 - binarySearch(cell2, Int.MAX_VALUE))
    }

    private fun binarySearch(cell: Long, edge: Int): Int {
        var low = 0
        var high = cells!!.size - 1
        while (low <= high) {
            val mid = low + high ushr 1
            val cmp = compare(cells!![mid], edges!![mid], cell, edge)
            if (cmp < 0) {
                low = mid + 1
            } else if (cmp > 0) {
                high = mid - 1
            } else {
                return mid
            }
        }
        return -(low + 1)
    }

    /**
     * Adds to candidateCrossings all the edges present in any ancestor of any
     * cell of cover, down to minimumS2LevelUsed. The cell->edge map is in the
     * variable mapping.
     */
    private fun getEdgesInParentCells(cover: List<S2CellId>, candidateCrossings: MutableSet<Int>) {
        // Find all parent cells of covering cells.
        val parentCells: MutableSet<S2CellId> = Sets.newHashSet()
        for (coverCell in cover) {
            for (parentLevel in coverCell.level() - 1 downTo minimumS2LevelUsed) {
                if (!parentCells.add(coverCell.parent(parentLevel))) {
                    break // cell is already in => parents are too.
                }
            }
        }

        // Put parent cell edge references into result.
        for (parentCell in parentCells) {
            val bounds = getEdges(parentCell.id.toLong(), parentCell.id.toLong())
            for (i in bounds[0] until bounds[1]) {
                candidateCrossings.add(edges!![i])
            }
        }
    }

    /**
     * Appends to candidateCrossings the edges that are fully contained in an S2
     * covering of edge. The covering of edge used is initially cover, but is
     * refined to eliminate quickly subcells that contain many edges but do not
     * intersect with edge.
     */
    private fun getEdgesInChildrenCells(a: S2Point, b: S2Point, cover: MutableList<S2CellId>,
                                        candidateCrossings: MutableSet<Int>) {
        // Put all edge references of (covering cells + descendant cells) into
        // result.
        // This relies on the natural ordering of S2CellIds.
        var children: Array<S2Cell>? = null
        while (!cover.isEmpty()) {
            val cell = cover.removeAt(cover.size - 1)
            var bounds = getEdges(cell.rangeMin().id.toLong(), cell.rangeMax().id.toLong())
            if (bounds[1] - bounds[0] <= 16) {
                for (i in bounds[0] until bounds[1]) {
                    candidateCrossings.add(edges!![i])
                }
            } else {
                // Add cells at this level
                bounds = getEdges(cell.id.toLong(), cell.id.toLong())
                for (i in bounds[0] until bounds[1]) {
                    candidateCrossings.add(edges!![i])
                }
                // Recurse on the children -- hopefully some will be empty.
                if (children == null) {
                    children = Array(4) { S2Cell() }
                }
                S2Cell(cell).subdivide(children)
                for (child in children) {
                    // TODO(user): Do the check for the four cells at once,
                    // as it is enough to check the four edges between the cells. At
                    // this time, we are checking 16 edges, 4 times too many.
                    //
                    // Note that given the guarantee of AppendCovering, it is enough
                    // to check that the edge intersect with the cell boundary as it
                    // cannot be fully contained in a cell.
                    if (edgeIntersectsCellBoundary(a, b, child)) {
                        cover.add(child.id())
                    }
                }
            }
        }
    }

    /*
   * An iterator on data edges that may cross a query edge (a,b). Create the
   * iterator, call getCandidates(), then hasNext()/next() repeatedly.
   *
   * The current edge in the iteration has index index(), goes between from()
   * and to().
   */
    class DataEdgeIterator(
            /**
             * The structure containing the data edges.
             */
            private val edgeIndex: S2EdgeIndex) {
        /**
         * Tells whether getCandidates() obtained the candidates through brute force
         * iteration or using the quad tree structure.
         */
        private var isBruteForce = false

        /**
         * Index of the current edge and of the edge before the last next() call.
         */
        private var currentIndex = 0

        /**
         * Cache of edgeIndex.getNumEdges() so that hasNext() doesn't make an extra
         * call
         */
        private var numEdges = 0

        /**
         * All the candidates obtained by getCandidates() when we are using a
         * quad-tree (i.e. isBruteForce = false).
         */
        var candidates: ArrayList<Int>

        /**
         * Index within array above. We have: currentIndex =
         * candidates.get(currentIndexInCandidates).
         */
        private var currentIndexInCandidates = 0

        /**
         * Initializes the iterator to iterate over a set of candidates that may
         * cross the edge (a,b).
         */
        fun getCandidates(a: S2Point, b: S2Point) {
            edgeIndex.predictAdditionalCalls(1)
            isBruteForce = !edgeIndex.isIndexComputed
            if (isBruteForce) {
                edgeIndex.incrementQueryCount()
                currentIndex = 0
                numEdges = edgeIndex.numEdges
            } else {
                candidates.clear()
                edgeIndex.findCandidateCrossings(a, b, candidates)
                currentIndexInCandidates = 0
                if (!candidates.isEmpty()) {
                    currentIndex = candidates[0]
                }
            }
        }

        /**
         * Index of the current edge in the iteration.
         */
        fun index(): Int {
            Preconditions.checkState(hasNext())
            return currentIndex
        }

        /**
         * False if there are no more candidates; true otherwise.
         */
        operator fun hasNext(): Boolean {
            return if (isBruteForce) {
                currentIndex < numEdges
            } else {
                currentIndexInCandidates < candidates.size
            }
        }

        /**
         * Iterate to the next available candidate.
         */
        operator fun next() {
            Preconditions.checkState(hasNext())
            if (isBruteForce) {
                ++currentIndex
            } else {
                ++currentIndexInCandidates
                if (currentIndexInCandidates < candidates.size) {
                    currentIndex = candidates[currentIndexInCandidates]
                }
            }
        }

        init {
            candidates = Lists.newArrayList()
        }
    }

    companion object {
        /**
         * Thicken the edge in all directions by roughly 1% of the edge length when
         * thickenEdge is true.
         */
        private const val THICKENING = 0.01

        /**
         * Threshold for small angles, that help lenientCrossing to determine whether
         * two edges are likely to intersect.
         */
        private const val MAX_DET_ERROR = 1e-14

        /**
         * Compares [cell1, edge1] to [cell2, edge2], by cell first and edge second.
         *
         * @return -1 if [cell1, edge1] is less than [cell2, edge2], 1 if [cell1,
         * edge1] is greater than [cell2, edge2], 0 otherwise.
         */
        private fun compare(cell1: Long, edge1: Int, cell2: Long, edge2: Int): Int {
            return if (cell1 < cell2) {
                -1
            } else if (cell1 > cell2) {
                1
            } else if (edge1 < edge2) {
                -1
            } else if (edge1 > edge2) {
                1
            } else {
                0
            }
        }

        /**
         * Returns the smallest cell containing all four points, or
         * [S2CellId.sentinel] if they are not all on the same face. The
         * points don't need to be normalized.
         */
        private fun containingCell(pa: S2Point, pb: S2Point, pc: S2Point, pd: S2Point): S2CellId {
            var a = fromPoint(pa)
            var b = fromPoint(pb)
            var c = fromPoint(pc)
            var d = fromPoint(pd)
            if (a.face() != b.face() || a.face() != c.face() || a.face() != d.face()) {
                return sentinel()
            }
            while (!a.equals(b) || !a.equals(c) || !a.equals(d)) {
                a = a.parent()
                b = b.parent()
                c = c.parent()
                d = d.parent()
            }
            return a
        }

        /**
         * Returns the smallest cell containing both points, or Sentinel if they are
         * not all on the same face. The points don't need to be normalized.
         */
        private fun containingCell(pa: S2Point, pb: S2Point): S2CellId {
            var a = fromPoint(pa)
            var b = fromPoint(pb)
            if (a.face() != b.face()) {
                return sentinel()
            }
            while (!a.equals(b)) {
                a = a.parent()
                b = b.parent()
            }
            return a
        }

        /**
         * Returns true if ab possibly crosses cd, by clipping tiny angles to zero.
         */
        private fun lenientCrossing(a: S2Point, b: S2Point, c: S2Point?, d: S2Point?): Boolean {
            // assert (S2.isUnitLength(a));
            // assert (S2.isUnitLength(b));
            // assert (S2.isUnitLength(c));
            val acb = crossProd(a, c!!).dotProd(b)
            val bda = crossProd(b, d!!).dotProd(a)
            if (Math.abs(acb) < MAX_DET_ERROR || Math.abs(bda) < MAX_DET_ERROR) {
                return true
            }
            if (acb * bda < 0) {
                return false
            }
            val cbd = crossProd(c, b).dotProd(d)
            val dac = crossProd(c, a).dotProd(c)
            return if (Math.abs(cbd) < MAX_DET_ERROR || Math.abs(dac) < MAX_DET_ERROR) {
                true
            } else acb * cbd >= 0 && acb * dac >= 0
        }

        /**
         * Returns true if the edge and the cell (including boundary) intersect.
         */
        private fun edgeIntersectsCellBoundary(a: S2Point, b: S2Point, cell: S2Cell): Boolean {
            val vertices = arrayOfNulls<S2Point>(4)
            for (i in 0..3) {
                vertices[i] = cell.getVertex(i)
            }
            for (i in 0..3) {
                val fromPoint = vertices[i]
                val toPoint = vertices[(i + 1) % 4]
                if (lenientCrossing(a, b, fromPoint, toPoint)) {
                    return true
                }
            }
            return false
        }
    }
}