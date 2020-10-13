package dilivia.s2.index

import com.google.common.collect.ComparisonChain
import dilivia.s2.Assertions
import dilivia.s2.Assertions.assertGE
import dilivia.s2.Assertions.assertTrue
import dilivia.s2.S2CellId
import dilivia.s2.index.S2CellIndex.LabelledCell
import dilivia.s2.index.S2CellIndex.NonEmptyRangeIterator
import dilivia.s2.index.S2CellIndex.RangeIterator
import dilivia.s2.region.S2Cap
import dilivia.s2.region.S2Cell
import dilivia.s2.region.S2CellUnion
import dilivia.s2.region.S2Region
import dilivia.s2.region.S2RegionCoverer
import mu.KotlinLogging
import java.util.*
import kotlin.collections.ArrayList
import kotlin.collections.HashSet

// S2ClosestCellQueryBase is a templatized class for finding the closest
// (cell_id, label) pairs in an S2CellIndex to a given target.  It is not
// intended to be used directly, but rather to serve as the implementation of
// various specialized classes with more convenient APIs (such as
// S2ClosestCellQuery).  It is flexible enough so that it can be adapted to
// compute maximum distances and even potentially Hausdorff distances.
//
// By using the appropriate options, this class can answer questions such as:
//
//  - Find the minimum distance between a cell collection A and a target B.
//  - Find all cells in collection A that are within a distance D of target B.
//  - Find the k cells of collection A that are closest to a given point P.
//
// The target is any class that implements the S2DistanceTarget interface.
// There are predefined targets for points, edges, S2Cells, S2CellUnions, and
// S2ShapeIndexes (arbitrary collctions of points, polylines, and polygons).
//
// The Distance template argument is used to represent distances.  Usually it
// is a thin wrapper around S1ChordAngle, but another distance type may be
// used as long as it implements the Distance concept described in
// s2distance_targets.h.  For example this can be used to measure maximum
// distances, to get more accuracy, or to measure non-spheroidal distances.
class S2ClosestCellQueryBase<T : Distance<T>> {

    // Default constructor; requires Init() to be called.
    constructor(distanceFactory: DistanceFactory<T>) {
        testedCells.add(LabelledCell(S2CellId.none(), -1))
        this.distanceFactory = distanceFactory
    }

    // Convenience constructor that calls Init().
    constructor(distanceFactory: DistanceFactory<T>, index: S2CellIndex) {
        this.distanceFactory = distanceFactory
        init(index)
    }

    val distanceFactory: DistanceFactory<T>
    private lateinit var index: S2CellIndex
    private lateinit var options: Options<T>
    private lateinit var target: S2DistanceTarget<T>


    // True if max_error() must be subtracted from priority queue cell distances
    // in order to ensure that such distances are measured conservatively.  This
    // is true only if the target takes advantage of max_error() in order to
    // return faster results, and 0 < max_error() < distance_limit_.
    private var useConservativeCellDistance: Boolean = false

    // For the optimized algorithm we precompute the top-level S2CellIds that
    // will be added to the priority queue.  There can be at most 6 of these
    // cells.  Essentially this is just a covering of the indexed cells.
    private var indexCovering: ArrayList<S2CellId> = ArrayList(6)

    // The distance beyond which we can safely ignore further candidate cells.
    // (Candidates that are exactly at the limit are ignored; this is more
    // efficient for UpdateMinDistance() and should not affect clients since
    // distance measurements have a small amount of error anyway.)
    //
    // Initially this is the same as the maximum distance specified by the user,
    // but it can also be updated by the algorithm (see MaybeAddResult).
    lateinit var distanceLimit: T

    // The current result set is stored in one of three ways:
    //
    //  - If max_results() == 1, the best result is kept in result_singleton_.
    //
    //  - If max_results() == kMaxMaxResults, results are appended to
    //    result_vector_ and sorted/uniqued at the end.
    //
    //  - Otherwise results are kept in a btree_set so that we can progressively
    //    reduce the distance limit once max_results() results have been found.
    //    (A priority queue is not sufficient because we need to be able to
    //    check whether a candidate cell is already in the result set.)
    //
    // TODO(ericv): Check whether it would be faster to use avoid_duplicates_
    // when result_set_ is used so that we could use a priority queue instead.
    private lateinit var resultSingleton: Result<T>
    private var resultVector: MutableList<Result<T>> = mutableListOf()
    private var resultSet: TreeSet<Result<T>> = TreeSet()

    // Used to iterate over the contents of an S2CellIndex range.  It is defined
    // here to take advantage of the fact that when multiple ranges are visited
    // in increasing order, duplicates can automatically be eliminated.
    private lateinit var contentsIt: S2CellIndex.ContentsIterator

    // When the results are stored in a btree_set (see above), usually
    // duplicates can be removed simply by inserting candidate cells in the
    // current result set.  However this is not true if Options::max_error() > 0
    // and the Target subtype takes advantage of this by returning suboptimal
    // distances.  This is because when UpdateMinDistance() is called with
    // different "min_dist" parameters (i.e., the distance to beat), the
    // implementation may return a different distance for the same cell.  Since
    // the btree_set is keyed by (distance, cell_id, label) this can create
    // duplicate results.
    //
    // The flag below is true when duplicates must be avoided explicitly.  This
    // is achieved by maintaining a separate set keyed by (cell_id, label) only,
    // and checking whether each edge is in that set before computing the
    // distance to it.
    //
    // TODO(ericv): Check whether it is faster to avoid duplicates by default
    // (even when Options::max_results() == 1), rather than just when we need to.
    private var avoidDuplicates: Boolean = true
    private val testedCells: MutableSet<LabelledCell> = HashSet()

    // The algorithm maintains a priority queue of unprocessed S2CellIds, sorted
    // in increasing order of distance from the target.
    private val queue: Deque<QueueEntry<T>> = ArrayDeque(16)

    // Temporaries, defined here to avoid multiple allocations / initializations.
    private val maxDistanceCovering: MutableList<S2CellId> = mutableListOf()
    private val intersectionWithMaxDistance: MutableList<S2CellId> = mutableListOf()

    // Initializes the query.
    // REQUIRES: ReInit() must be called if "index" is modified.
    fun init(index: S2CellIndex) {
        this.index = index
        this.contentsIt = S2CellIndex.ContentsIterator(index)
        reInit()
    }

    // Reinitializes the query.  This method must be called whenever the
    // underlying index is modified.
    fun reInit() {
        indexCovering.clear()
        distanceLimit = distanceFactory.infinity()
    }

    // Return a reference to the underlying S2CellIndex.
    fun index(): S2CellIndex = index

    // Returns the closest (cell_id, label) pairs to the given target that
    // satisfy the given options.  This method may be called multiple times.
    fun findClosestCells(target: S2DistanceTarget<T>, options: Options<T>): List<Result<T>> {
        val results = mutableListOf<Result<T>>()
        findClosestCells(target, options, results)
        return results
    }

    // This version can be more efficient when this method is called many times,
    // since it does not require allocating a new vector on each call.
    fun findClosestCells(target: S2DistanceTarget<T>, options: Options<T>, results: MutableList<Result<T>>) {
        findClosestCellsInternal(target, options)
        results.clear()
        if (options.getMaxResults() == 1) {
            if (!resultSingleton.isEmpty()) {
                results.add(resultSingleton)
            }
        } else if (options.getMaxResults() == Options.kMaxMaxResults) {
            results.addAll(resultVector)
            results.sort()
            resultVector.clear()
        } else {
            results.addAll(resultSet)
            resultSet.clear()
        }
    }

    // Convenience method that returns exactly one (cell_id, label) pair.  If no
    // cells satisfy the given search criteria, then a Result with
    // distance() == Infinity() and is_empty() == true is returned.
    //
    // REQUIRES: options.max_results() == 1
    fun findClosestCell(target: S2DistanceTarget<T>, options: Options<T>): Result<T> {
        Assertions.assertEQ(options.getMaxResults(), 1)
        findClosestCellsInternal(target, options)
        return resultSingleton
    }

    // Options that control the set of cells returned.  Note that by default
    // *all* cells are returned, so you will always want to set either the
    // max_results() option or the max_distance() option (or both).
    //
    // This class is also available as S2ClosestCellQueryBase<Data>::Options.
    //
    // The Distance template argument is described below.
    open class Options<T : Distance<T>>(
          
            val distanceFactory: DistanceFactory<T>,
            
            // Specifies that at most "max_results" cells should be returned.
            //
            // REQUIRES: max_results >= 1
            // DEFAULT: kMaxMaxResults
            private var maxResults: Int = kMaxMaxResults,

            // Specifies that only cells whose distance to the target is less than
            // "max_distance" should be returned.
            //
            // Note that cells whose distance is exactly equal to "max_distance" are
            // not returned.  In most cases this doesn't matter (since distances are
            // not computed exactly in the first place), but if such cells are needed
            // then you can retrieve them by specifying "max_distance" as the next
            // largest representable Distance.  For example, if Distance is an
            // S1ChordAngle then you can specify max_distance.Successor().
            //
            // DEFAULT: Distance::Infinity()
            var maxDistance: T = distanceFactory.infinity(),

            // Specifies that cells up to max_error() further away than the true
            // closest cells may be substituted in the result set, as long as such
            // cells satisfy all the remaining search criteria (such as max_distance).
            // This option only has an effect if max_results() is also specified
            // otherwise all cells closer than max_distance() will always be returned.
            //
            // Note that this does not affect how the distance between cells is
            // computed; it simply gives the algorithm permission to stop the search
            // early as soon as the best possible improvement drops below max_error().
            //
            // This can be used to implement distance predicates efficiently.  For
            // example, to determine whether the minimum distance is less than D, the
            // IsDistanceLess() method sets max_results() == 1 and max_distance() ==
            // max_error() == D.  This causes the algorithm to terminate as soon as it
            // finds any cell whose distance is less than D, rather than continuing to
            // search for a cell that is even closer.
            //
            // DEFAULT: Distance::Delta::Zero()
            var maxError: Delta = Delta.zero,

            // Specifies that cells must intersect the given S2Region.  "region" is
            // owned by the caller and must persist during the lifetime of this
            // object.  The value may be changed between calls to FindClosestPoints(),
            // or reset by calling set_region(nullptr).
            //
            // Note that if you want to set the region to a disc around a target
            // point, it is faster to use a PointTarget with set_max_distance()
            // instead.  You can also call both methods, e.g. to set a maximum
            // distance and also require that cells lie within a given rectangle.
            var region: S2Region? = null,

            // Specifies that distances should be computed by examining every cell
            // rather than using the S2ShapeIndex.  This is useful for testing,
            // benchmarking, and debugging.
            //
            // DEFAULT: false
            var useBruteForce: Boolean = false
    ): Cloneable {

        fun setMaxResults(value: Int) {
            assertGE(value, 1)
            this.maxResults = value
        }

        fun getMaxResults(): Int = maxResults
        
        override fun clone(): Options<T> {
            return Options(
                    distanceFactory, maxResults, maxDistance, maxError, region, useBruteForce
            )
        }

        companion object {

            const val kMaxMaxResults = Int.MAX_VALUE

        }
    }

    // Each "Result" object represents a closest (cell_id, label) pair.
    data class Result<T : Distance<T>>(
            val distance: T,
            val cellId: S2CellId = S2CellId.none(),
            val label: Int = -1
    ) : Comparable<Result<T>> {

        // Returns true if this Result object does not refer to any cell.
        // (The only case where an empty Result is returned is when the
        // FindClosestCell() method does not find any cells that meet the
        // specified criteria.)
        fun isEmpty(): Boolean = cellId == S2CellId.none()

        override fun compareTo(other: Result<T>): Int {
            return ComparisonChain.start()
                    .compare(distance, other.distance)
                    .compare(cellId, other.cellId)
                    .compare(label, other.label)
                    .result()
        }

    }

    data class QueueEntry<T : Distance<T>>(
            // A lower bound on the distance from the target to "id".  This is the key
            // of the priority queue.
            val distance: T,

            // The cell being queued.
            val id: S2CellId
    ) : Comparable<QueueEntry<T>> {

        // The priority queue returns the largest elements first, so we want the
        // "largest" entry to have the smallest distance.
        override fun compareTo(other: QueueEntry<T>): Int = -distance.compareTo(other.distance)

    }

    companion object {

        private val logger = KotlinLogging.logger(S2ClosestCellQueryBase::class.java.name)

        // The minimum number of ranges that a cell must contain to enqueue it
        // rather than processing its contents immediately.
        const val kMinRangesToEnqueue = 6
    }


    private fun findClosestCellsInternal(target: S2DistanceTarget<T>, options: Options<T>) {
        this.target = target
        this.options = options

        testedCells.clear()
        contentsIt.clear()
        distanceLimit = options.maxDistance
        resultSingleton = Result(distanceFactory.infinity())
        assertTrue(resultVector.isEmpty())
        assertTrue(resultSet.isEmpty())
        assertGE(target.maxBruteForceIndexSize(), 0)
        if (distanceLimit == distanceFactory.zero()) return

        if (options.getMaxResults() == Options.kMaxMaxResults &&
                options.maxDistance == distanceFactory.infinity() &&
                options.region == null) {
            logger.warn { "Returning all cells (max_results/max_distance/region not set)" }
        }

        // If max_error() > 0 and the target takes advantage of this, then we may
        // need to adjust the distance estimates to the priority queue cells to
        // ensure that they are always a lower bound on the true distance.  For
        // example, suppose max_distance == 100, max_error == 30, and we compute the
        // distance to the target from some cell C0 as d(C0) == 80.  Then because
        // the target takes advantage of max_error(), the true distance could be as
        // low as 50.  In order not to miss edges contained by such cells, we need
        // to subtract max_error() from the distance estimates.  This behavior is
        // controlled by the use_conservative_cell_distance_ flag.
        //
        // However there is one important case where this adjustment is not
        // necessary, namely when max_distance() < max_error().  This is because
        // max_error() only affects the algorithm once at least max_edges() edges
        // have been found that satisfy the given distance limit.  At that point,
        // max_error() is subtracted from distance_limit_ in order to ensure that
        // any further matches are closer by at least that amount.  But when
        // max_distance() < max_error(), this reduces the distance limit to 0,
        // i.e. all remaining candidate cells and edges can safely be discarded.
        // (Note that this is how IsDistanceLess() and friends are implemented.)
        //
        // Note that Distance::Delta only supports operator==.
        val targetUsesMaxError = (options.maxError != Delta.zero && target.setMaxError(options.maxError))

        // Note that we can't compare max_error() and distance_limit_ directly
        // because one is a Delta and one is a Distance.  Instead we subtract them.
        useConservativeCellDistance = targetUsesMaxError &&
                (distanceLimit == distanceFactory.infinity() || distanceFactory.zero() < distanceLimit - options.maxError)

        // Use the brute force algorithm if the index is small enough.
        if (options.useBruteForce || index.numCells() <= target.maxBruteForceIndexSize()) {
            avoidDuplicates = false
            findClosestCellsBruteForce()
        } else {
            // If the target takes advantage of max_error() then we need to avoid
            // duplicate edges explicitly.  (Otherwise it happens automatically.)
            avoidDuplicates = (targetUsesMaxError && options.getMaxResults() > 1)
            findClosestCellsOptimized()
        }
    }

    private fun findClosestCellsBruteForce() {
        val iter = S2CellIndex.CellIterator(index)
        while (!iter.done()) {
            maybeAddResult(iter.cellId(), iter.label())
            iter.next()
        }
    }

    private fun findClosestCellsOptimized() {
        initQueue()
        while (!queue.isEmpty()) {
            // We need to copy the top entry before removing it, and we need to remove
            // it before adding any new entries to the queue.
            val entry = queue.pollFirst()
            // Work around weird parse error in gcc 4.9 by using a local variable for
            // entry.distance.
            val distance = entry.distance
            if (distance >= distanceLimit) {
                queue.clear()  // Clear any remaining entries.
                break
            }
            var child = entry.id.childBegin()
            // We already know that it has too many cells, so process its children.
            // Each child may either be processed directly or enqueued again.  The
            // loop is optimized so that we don't seek unnecessarily.
            var seek = true
            val range = NonEmptyRangeIterator(index)
            for (i in 0..3) {
                seek = processOrEnqueue(child, range, seek)
                child = child.next()
            }
        }
    }

    private fun initQueue() {
        assertTrue(queue.isEmpty())

        // Optimization: rather than starting with the entire index, see if we can
        // limit the search region to a small disc.  Then we can find a covering for
        // that disc and intersect it with the covering for the index.  This can
        // save a lot of work when the search region is small.
        val cap = target.getCapBound()
        if (cap.isEmpty) return  // Empty target.
        if (options.getMaxResults() == 1) {
            // If the user is searching for just the closest cell, we can compute an
            // upper bound on search radius by seeking to the center of the target's
            // bounding cap and looking at the contents of that leaf cell range.  If
            // the range intersects any cells, then the distance is zero.  Otherwise
            // we can still look at the two neighboring ranges, and use the minimum
            // distance to any cell in those ranges as an upper bound on the search
            // radius.  These cells may wind up being processed twice, but in general
            // this is still faster.
            //
            // First check the range containing or immediately following "center".
            val range = NonEmptyRangeIterator(index)
            val target = S2CellId.fromPoint(cap.center)
            range.seek(target)
            addRange(range)
            if (distanceLimit == distanceFactory.zero()) return

            // If the range immediately follows "center" (rather than containing it),
            // then check the previous non-empty range as well.
            if (range.startId() > target && range.prev()) {
                addRange(range)
                if (distanceLimit == distanceFactory.zero()) return
            }
        }

        // We start with a covering of the set of indexed cells, then intersect it
        // with the maximum search radius disc (if any).
        //
        // Note that unlike S2ClosestPointQuery, we can't also intersect with the
        // given region (if any).  This is because the index cells in the result are
        // only required to intersect the region.  This means that an index cell that
        // intersects the region's covering may be much closer to the target than the
        // covering itself, which means that we cannot use the region's covering to
        // restrict the search.
        //
        // TODO(ericv): If this feature becomes important, this could be fixed by
        // (1) computing a covering of the region, (2) looking up any index cells
        // that contain each covering cell by seeking to covering_cell.range_min(),
        // (3) replacing each covering cell by the largest such cell (if any), and
        // (4) normalizing the result.
        if (indexCovering.isEmpty()) initCovering()
        var initialCells = indexCovering.toList()
        if (distanceLimit < distanceFactory.infinity()) {
            val coverer = S2RegionCoverer(maxCells = 4)
            val radius = cap.radius + distanceLimit.getChordAngleBound()
            val searchCap = S2Cap(cap.center, radius)
            coverer.getFastCovering(searchCap, maxDistanceCovering)
            S2CellUnion.getIntersection(initialCells, maxDistanceCovering, intersectionWithMaxDistance)
            initialCells = intersectionWithMaxDistance.toList()
        }
        val range = NonEmptyRangeIterator(index)
        for (i in initialCells.indices) {
            val id = initialCells[i]
            val seek = (i == 0) || id.rangeMin() >= range.limitId()
            processOrEnqueue(id, range, seek)
            if (range.done()) break
        }
    }

    private fun initCovering() {
        // Compute the "index covering", which is a small number of S2CellIds that
        // cover the indexed cells.  There are two cases:
        //
        //  - If the index spans more than one face, then there is one covering cell
        // per spanned face, just big enough to cover the indexed cells on that face.
        //
        //  - If the index spans only one face, then we find the smallest cell "C"
        // that covers the indexed cells on that face (just like the case above).
        // Then for each of the 4 children of "C", if the child contains any index
        // cells then we create a covering cell that is big enough to just fit
        // those indexed cells (i.e., shrinking the child as much as possible to fit
        // its contents).  This essentially replicates what would happen if we
        // started with "C" as the covering cell, since "C" would immediately be
        // split, except that we take the time to prune the children further since
        // this will save work on every subsequent query.
        indexCovering.ensureCapacity(6)
        val it = NonEmptyRangeIterator(index)
        var last = NonEmptyRangeIterator(index)
        it.begin()
        last.finish()
        if (!last.prev()) return  // Empty index.
        val indexLastId = last.limitId().previous()
        if (it.startId() != last.startId()) {
            // The index contains at least two distinct S2CellIds (because otherwise
            // there would only be one non-empty range).  Choose a level such that the
            // entire index can be spanned with at most 6 cells (if the index spans
            // multiple faces) or 4 cells (it the index spans a single face).
            val level = it.startId().getCommonAncestorLevel(indexLastId) + 1

            // Visit each potential covering cell except the last (handled below).
            val startId = it.startId().parent(level)
            val lastId = indexLastId.parent(level)
            var id = startId
            while (id != lastId) {
                // Skip any covering cells that don't contain an indexed range.
                if (id.rangeMax() < it.startId()) continue

                // Find the indexed range contained by this covering cell and then
                // shrink the cell if necessary so that it just covers this range.
                val cellFirstId = it.startId()
                it.seek(id.rangeMax().next())
                // Find the last leaf cell covered by the previous non-empty range.
                last = it
                last.prev()
                addInitialRange(cellFirstId, last.limitId().previous())
                id = id.next()
            }
        }
        addInitialRange(it.startId(), indexLastId)
    }

    private fun addInitialRange(firstId: S2CellId, lastId: S2CellId) {
        // Add the lowest common ancestor of the given range.
        val level = firstId.getCommonAncestorLevel(lastId)
        assertGE(level, 0)
        indexCovering.add(firstId.parent(level))
    }

    // TODO(ericv): Consider having this method return false when distance_limit_
    // is reduced to zero, and terminating any calling loops early.
    private fun maybeAddResult(cell_id: S2CellId, label: Label) {
        if (avoidDuplicates && !testedCells.add(LabelledCell(cell_id, label))) {
            return
        }

        // TODO(ericv): It may be relatively common to add the same S2CellId
        // multiple times with different labels.  This could be optimized by
        // remembering the last "cell_id" argument and its distance.  However this
        // may not be beneficial when Options::max_results() == 1, for example.
        val cell = S2Cell(cell_id)
        val distance = distanceFactory.distance(distanceLimit)
        if (!target.updateMinDistance(cell, distance)) return

        val region = options.region
        if (region != null && !region.mayIntersect(cell)) return

        val result = Result(distance, cell_id, label)
        if (options.getMaxResults() == 1) {
            // Optimization for the common case where only the closest cell is wanted.
            resultSingleton = result
            distanceLimit = result.distance - options.maxError
        } else if (options.getMaxResults() == Options.kMaxMaxResults) {
            resultVector.add(result)  // Sort/unique at end.
        } else {
            // Add this cell to result_set_.  Note that even if we already have enough
            // edges, we can't erase an element before insertion because the "new"
            // edge might in fact be a duplicate.
            resultSet.add(result)
            val size = resultSet.size
            if (size >= options.getMaxResults()) {
                if (size > options.getMaxResults()) {
                    resultSet.pollLast()
                }
                distanceLimit = resultSet.last().distance - options.maxError
            }
        }
    }

    private fun processOrEnqueue(id: S2CellId, iter: NonEmptyRangeIterator, seek: Boolean): Boolean {
        if (seek) iter.seek(id.rangeMin())
        val last = id.rangeMax()
        if (iter.startId() > last) {
            return false  // No need to seek to next child.
        }
        // If this cell intersects at least "kMinRangesToEnqueue" leaf cell ranges
        // (including ranges whose contents are empty), then enqueue it.  We test
        // this by advancing (n - 1) ranges and checking whether that range also
        // intersects this cell.
        val maxIt = RangeIterator(iter)
        if (maxIt.advance(kMinRangesToEnqueue - 1) && maxIt.startId() <= last) {
            // This cell intersects at least kMinRangesToEnqueue ranges, so enqueue it.
            val cell = S2Cell(id)
            var distance = distanceFactory.distance(distanceLimit)
            // We check "region_" second because it may be relatively expensive.
            val region = options.region
            if (target.updateMinDistance(cell, distance) && (region == null || region.mayIntersect(cell))) {
                if (useConservativeCellDistance) {
                    // Ensure that "distance" is a lower bound on distance to the cell.
                    distance -= options.maxError
                }
                queue.push(QueueEntry(distance, id))
            }
            return true  // Seek to next child.
        }
        // There were few enough ranges that we might as well process them now.
        while (iter.startId() <= last) {
            addRange(RangeIterator(iter))
            iter.next()
        }
        return false  // No need to seek to next child.
    }

    private fun addRange(range: RangeIterator) {
        contentsIt.startUnion(range)
        while (!contentsIt.done()) {
            maybeAddResult(contentsIt.cellId(), contentsIt.label())
            contentsIt.next()
        }
    }

}

