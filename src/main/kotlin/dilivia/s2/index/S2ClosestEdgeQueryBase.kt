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
package dilivia.s2.index

import com.google.common.collect.ComparisonChain
import dilivia.s2.Assertions
import dilivia.s2.S2CellId
import dilivia.s2.S2Point
import dilivia.s2.index.shape.CellRelation
import dilivia.s2.index.shape.InitialPosition
import dilivia.s2.index.shape.IteratorBase
import dilivia.s2.index.shape.S2ShapeIndex
import dilivia.s2.index.shape.S2ShapeIndexCell
import dilivia.s2.index.shape.S2ShapeIndexCellIterator
import dilivia.s2.region.S2Cap
import dilivia.s2.region.S2Cell
import dilivia.s2.region.S2CellUnion
import dilivia.s2.region.S2RegionCoverer
import dilivia.s2.shape.S2Shape
import dilivia.s2.shape.ShapeEdgeId
import mu.KotlinLogging
import java.util.*
import kotlin.collections.ArrayList


// S2ClosestEdgeQueryBase is a templatized class for finding the closest
// edge(s) between two geometries.  It is not intended to be used directly,
// but rather to serve as the implementation of various specialized classes
// with more convenient APIs (such as S2ClosestEdgeQuery).  It is flexible
// enough so that it can be adapted to compute maximum distances and even
// potentially Hausdorff distances.
//
// By using the appropriate options, this class can answer questions such as:
//
//  - Find the minimum distance between two geometries A and B.
//  - Find all edges of geometry A that are within a distance D of geometry B.
//  - Find the k edges of geometry A that are closest to a given point P.
//
// You can also specify whether polygons should include their interiors (i.e.,
// if a point is contained by a polygon, should the distance be zero or should
// it be measured to the polygon boundary?)
//
// The input geometries may consist of any number of points, polylines, and
// polygons (collectively referred to as "shapes").  Shapes do not need to be
// disjoint; they may overlap or intersect arbitrarily.  The implementation is
// designed to be fast for both simple and complex geometries.
//
// The Distance template argument is used to represent distances.  Usually it
// is a thin wrapper around S1ChordAngle, but another distance type may be
// used as long as it implements the Distance concept described in
// s2distance_targets.h.  For example this can be used to measure maximum
// distances, to get more accuracy, or to measure non-spheroidal distances.
class S2ClosestEdgeQueryBase<T : Distance<T>> {

    // Default constructor; requires Init() to be called.
    constructor(distanceFactory: DistanceFactory<T>) {
        this.distanceFactory = distanceFactory
    }

    // Convenience constructor that calls Init().
    constructor(distanceFactory: DistanceFactory<T>, index: S2ShapeIndex) {
        this.distanceFactory = distanceFactory
        init(index)
    }

    val distanceFactory: DistanceFactory<T>

    private lateinit var index: S2ShapeIndex
    private lateinit var options: Options<T>
    private lateinit var target: S2DistanceTarget<T>

    // True if max_error() must be subtracted from priority queue cell distances
    // in order to ensure that such distances are measured conservatively.  This
    // is true only if the target takes advantage of max_error() in order to
    // return faster results, and 0 < max_error() < distance_limit_.
    private var useConservativeCellDistance: Boolean = false

    // For the optimized algorihm we precompute the top-level S2CellIds that
    // will be added to the priority queue.  There can be at most 6 of these
    // cells.  Essentially this is just a covering of the indexed edges, except
    // that we also store pointers to the corresponding S2ShapeIndexCells to
    // reduce the number of index seeks required.
    //
    // The covering needs to be stored in a std::vector so that we can use
    // S2CellUnion::GetIntersection().
    private var indexCovering: ArrayList<S2CellId> = ArrayList(6)
    private var index_cells: MutableList<S2ShapeIndexCell?> = ArrayList(6)

    // The decision about whether to use the brute force algorithm is based on
    // counting the total number of edges in the index.  However if the index
    // contains a large number of shapes, this in itself might take too long.
    // So instead we only count edges up to (max_brute_force_index_size() + 1)
    // for the current target type (stored as index_num_edges_limit_).
    private var index_num_edges: Int = 0
    private var index_num_edges_limit: Int = 0

    // The distance beyond which we can safely ignore further candidate edges.
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
    //  - If max_results() == "infinity", results are appended to result_vector_
    //    and sorted/uniqued at the end.
    //
    //  - Otherwise results are kept in a btree_set so that we can progressively
    //    reduce the distance limit once max_results() results have been found.
    //    (A priority queue is not sufficient because we need to be able to
    //    check whether a candidate edge is already in the result set.)
    //
    // TODO(ericv): Check whether it would be faster to use avoid_duplicates_
    // when result_set_ is used so that we could use a priority queue instead.
    private lateinit var resultSingleton: Result<T>
    private var resultVector: MutableList<Result<T>> = mutableListOf()
    private var resultSet: TreeSet<Result<T>> = TreeSet()

    // When the result edges are stored in a btree_set (see above), usually
    // duplicates can be removed simply by inserting candidate edges in the
    // current set.  However this is not true if Options::max_error() > 0 and
    // the Target subtype takes advantage of this by returning suboptimal
    // distances.  This is because when UpdateMinDistance() is called with
    // different "min_dist" parameters (i.e., the distance to beat), the
    // implementation may return a different distance for the same edge.  Since
    // the btree_set is keyed by (distance, shape_id, edge_id) this can create
    // duplicate edges in the results.
    //
    // The flag below is true when duplicates must be avoided explicitly.  This
    // is achieved by maintaining a separate set keyed by (shape_id, edge_id)
    // only, and checking whether each edge is in that set before computing the
    // distance to it.
    //
    // TODO(ericv): Check whether it is faster to avoid duplicates by default
    // (even when Options::max_results() == 1), rather than just when we need to.
    private var avoidDuplicates: Boolean = true
    private val testedEdges: MutableSet<ShapeEdgeId> = HashSet()

    // The algorithm maintains a priority queue of unprocessed S2CellIds, sorted
    // in increasing order of distance from the target.
    private val queue: Deque<QueueEntry<T>> = ArrayDeque(16)

    // Temporaries, defined here to avoid multiple allocations / initializations.
    private lateinit var iter: S2ShapeIndexCellIterator
    private val maxDistanceCovering: MutableList<S2CellId> = mutableListOf()
    private val initial_cells: MutableList<S2CellId> = mutableListOf()

    // Initializes the query.
    // REQUIRES: ReInit() must be called if "index" is modified.
    fun init(index: S2ShapeIndex) {
        this.index = index
        reInit()
    }

    // Reinitializes the query.  This method must be called whenever the
    // underlying index is modified.
    fun reInit() {
        index_num_edges = 0
        index_num_edges_limit = 0
        indexCovering.clear()
        index_cells.clear()
        // We don't initialize iter_ here to make queries on small indexes a bit
        // faster (i.e., where brute force is used).
    }

    // Returns a reference to the underlying S2ShapeIndex.
    fun index(): S2ShapeIndex = index

    // Returns the closest edges to the given target that satisfy the given
    // options.  This method may be called multiple times.
    //
    // Note that if options().include_interiors() is true, the result vector may
    // include some entries with edge_id == -1.  This indicates that the target
    // intersects the indexed polygon with the given shape_id.
    fun findClosestEdges(target: S2DistanceTarget<T>, options: Options<T>): List<Result<T>> {
        val results = mutableListOf<Result<T>>()
        findClosestEdges(target, options, results)
        return results;
    }

    // This version can be more efficient when this method is called many times,
    // since it does not require allocating a new vector on each call.
    fun findClosestEdges(target: S2DistanceTarget<T>, options: Options<T>, results: MutableList<Result<T>>) {
        findClosestEdgesInternal(target, options)
        results.clear()
        if (options.getMaxResults() == 1) {
            if (resultSingleton.shapeId >= 0) {
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

    // Convenience method that returns exactly one edge.  If no edges satisfy
    // the given search criteria, then a Result with distance == Infinity() and
    // shape_id == edge_id == -1 is returned.
    //
    // Note that if options.include_interiors() is true, edge_id == -1 is also
    // used to indicate that the target intersects an indexed polygon (but in
    // that case distance == Zero() and shape_id >= 0).
    //
    // REQUIRES: options.max_results() == 1
    fun findClosestEdge(target: S2DistanceTarget<T>, options: Options<T>): Result<T> {
        check(options.getMaxResults() == 1)
        findClosestEdgesInternal(target, options)
        return resultSingleton
    }

    // Options that control the set of edges returned.  Note that by default
    // *all* edges are returned, so you will always want to set either the
    // max_results() option or the max_distance() option (or both).
    open class Options<T : Distance<T>>(

            val distanceFactory: DistanceFactory<T>,

            // Specifies that at most "max_results" edges should be returned.
            //
            // REQUIRES: max_results >= 1
            // DEFAULT: kMaxMaxResults
            private var maxResults: Int = S2ClosestCellQueryBase.Options.kMaxMaxResults,

            // Specifies that only edges whose distance to the target is less than
            // "max_distance" should be returned.
            //
            // Note that edges whose distance is exactly equal to "max_distance" are
            // not returned.  In most cases this doesn't matter (since distances are
            // not computed exactly in the first place), but if such edges are needed
            // then you can retrieve them by specifying "max_distance" as the next
            // largest representable Distance.  For example, if Distance is an
            // S1ChordAngle then you can specify max_distance.Successor().
            //
            // DEFAULT: Distance::Infinity()
            var maxDistance: T = distanceFactory.infinity(),

            // Specifies that edges up to max_error() further away than the true
            // closest edges may be substituted in the result set, as long as such
            // edges satisfy all the remaining search criteria (such as max_distance).
            // This option only has an effect if max_results() is also specified;
            // otherwise all edges closer than max_distance() will always be returned.
            //
            // Note that this does not affect how the distance between edges is
            // computed; it simply gives the algorithm permission to stop the search
            // early as soon as the best possible improvement drops below max_error().
            //
            // This can be used to implement distance predicates efficiently.  For
            // example, to determine whether the minimum distance is less than D, set
            // max_results() == 1 and max_distance() == max_error() == D.  This causes
            // the algorithm to terminate as soon as it finds any edge whose distance
            // is less than D, rather than continuing to search for an edge that is
            // even closer.
            //
            // DEFAULT: Distance::Delta::Zero()
            var maxError: Delta = Delta.zero,

            // Specifies that polygon interiors should be included when measuring
            // distances.  In other words, polygons that contain the target should
            // have a distance of zero.  (For targets consisting of multiple connected
            // components, the distance is zero if any component is contained.)  This
            // is indicated in the results by returning a (shape_id, edge_id) pair
            // with edge_id == -1, i.e. this value denotes the polygons's interior.
            //
            // Note that for efficiency, any polygon that intersects the target may or
            // may not have an (edge_id == -1) result.  Such results are optional
            // because in that case the distance to the polygon is already zero.
            //
            // DEFAULT: true
            var includeInteriors: Boolean = true,

            // Specifies that distances should be computed by examining every edge
            // rather than using the S2ShapeIndex.  This is useful for testing,
            // benchmarking, and debugging.
            //
            // DEFAULT: false
            var useBruteForce: Boolean = false

    ): Cloneable {

        fun setMaxResults(value: Int) {
            Assertions.assertGE(value, 1)
            this.maxResults = value
        }

        fun getMaxResults(): Int = maxResults

        override fun clone(): S2ClosestEdgeQueryBase.Options<T> {
            return S2ClosestEdgeQueryBase.Options(
                    distanceFactory, maxResults, maxDistance, maxError, includeInteriors, useBruteForce
            )
        }

        companion object {
            const val kMaxMaxResults = Int.MAX_VALUE
        }

    }

    // Each "Result" object represents a closest edge.  Note the following
    // special cases:
    //
    //  - (shape_id() >= 0) && (edge_id() < 0) represents the interior of a shape.
    //    Such results may be returned when options.include_interiors() is true.
    //    Such results can be identified using the is_interior() method.
    //
    //  - (shape_id() < 0) && (edge_id() < 0) is returned by `FindClosestEdge`
    //    to indicate that no edge satisfies the given query options.  Such
    //    results can be identified using is_empty() method.
    data class Result<T : Distance<T>>(
            val distance: T,           // The distance from the target to this edge.
            val shapeId: Int = -1,     // Identifies an indexed shape.
            val edgeId: Int = -1       // Identifies an edge within the shape.
    ) : Comparable<Result<T>> {

        // Returns true if this Result object represents the interior of a shape.
        // (Such results may be returned when options.include_interiors() is true.)
        fun isInterior(): Boolean = shapeId >= 0 && edgeId < 0

        // Returns true if this Result object indicates that no edge satisfies the
        // given query options.  (This result is only returned in one special
        // case, namely when FindClosestEdge() does not find any suitable edges.
        // It is never returned by methods that return a vector of results.)
        fun isEmpty(): Boolean = shapeId < 0

        // Compares edges first by distance, then by (shape_id, edge_id).
        override fun compareTo(other: Result<T>): Int = ComparisonChain.start()
                .compare(distance, other.distance)
                .compare(shapeId, other.shapeId)
                .compare(edgeId, other.edgeId)
                .result()

    }

    data class QueueEntry<T : Distance<T>>(
            // A lower bound on the distance from the target to "id".  This is the key
            // of the priority queue.
            val distance: T,

            // The cell being queued.
            val id: S2CellId,

            // If "id" belongs to the index, this field stores the corresponding
            // S2ShapeIndexCell.  Otherwise "id" is a proper ancestor of one or more
            // S2ShapeIndexCells and this field stores nullptr.  The purpose of this
            // field is to avoid an extra Seek() when the queue entry is processed.
            val indexCell: S2ShapeIndexCell? = null

    ) : Comparable<QueueEntry<T>> {

        // The priority queue returns the largest elements first, so we want the
        // "largest" entry to have the smallest distance.
        override fun compareTo(other: QueueEntry<T>): Int = -distance.compareTo(other.distance)

    }

    companion object {
        private val logger = KotlinLogging.logger(S2ClosestEdgeQueryBase::class.java.name)

        // Return the number of edges in the given index cell.
        private fun countEdges(cell: S2ShapeIndexCell): Int {
            var count = 0;
            for (s in 0 until cell.numClipped) {
                count += cell.clipped(s).numEdges
            }
            return count
        }
    }


    private fun findClosestEdgesInternal(target: S2DistanceTarget<T>, options: Options<T>) {
        this.target = target
        this.options = options

        testedEdges.clear()
        distanceLimit = options.maxDistance
        resultSingleton = Result(distance = distanceFactory.infinity())
        check(resultVector.isEmpty())
        check(resultSet.isEmpty())
        check(target.maxBruteForceIndexSize() >= 0)
        if (distanceLimit == distanceFactory.zero()) return

        if (options.getMaxResults() == Options.kMaxMaxResults && options.maxDistance == distanceFactory.infinity()) {
            logger.warn { "Returning all edges (max_results/max_distance not set)" }
        }

        if (options.includeInteriors) {
            val shapeIds = mutableSetOf<Int>()
            target.visitContainingShapes(index, object : S2DistanceTarget.ShapeVisitor {

                override fun visit(containing_shape: S2Shape, target_point: S2Point): Boolean {
                    shapeIds.add(containing_shape.id)
                    return shapeIds.size < options.getMaxResults()
                }

            })
            for (shape_id in shapeIds) {
                addResult(Result(distanceFactory.zero(), shape_id, -1));
            }
            if (distanceLimit == distanceFactory.zero()) return
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
        // max_error() only affects the algorithm once at least max_results() edges
        // have been found that satisfy the given distance limit.  At that point,
        // max_error() is subtracted from distance_limit_ in order to ensure that
        // any further matches are closer by at least that amount.  But when
        // max_distance() < max_error(), this reduces the distance limit to 0,
        // i.e. all remaining candidate cells and edges can safely be discarded.
        // (Note that this is how IsDistanceLess() and friends are implemented.)
        //
        // Note that Distance::Delta only supports operator==.
        val target_uses_max_error = (!(options.maxError == Delta.zero) && target.setMaxError(options.maxError))

        // Note that we can't compare max_error() and distance_limit_ directly
        // because one is a Delta and one is a Distance.  Instead we subtract them.
        useConservativeCellDistance = target_uses_max_error &&
                (distanceLimit == distanceFactory.infinity() || distanceFactory.zero() < distanceLimit - options.maxError)

        // Use the brute force algorithm if the index is small enough.  To avoid
        // spending too much time counting edges when there are many shapes, we stop
        // counting once there are too many edges.  We may need to recount the edges
        // if we later see a target with a larger brute force edge threshold.
        val min_optimized_edges = target.maxBruteForceIndexSize() + 1
        if (min_optimized_edges > index_num_edges_limit && index_num_edges >= index_num_edges_limit) {
            index_num_edges = S2CountEdges.countEdgesUpTo(index, min_optimized_edges)
            index_num_edges_limit = min_optimized_edges
        }

        if (options.useBruteForce || index_num_edges < min_optimized_edges) {
            // The brute force algorithm considers each edge exactly once.
            avoidDuplicates = false
            findClosestEdgesBruteForce()
        } else {
            // If the target takes advantage of max_error() then we need to avoid
            // duplicate edges explicitly.  (Otherwise it happens automatically.)
            avoidDuplicates = (target_uses_max_error && options.getMaxResults() > 1)
            findClosestEdgesOptimized();
        }
    }

    private fun findClosestEdgesBruteForce() {
        for (shape in index.shapeIterator()) {
            if (shape == null) continue
            val numEdges = shape.numEdges
            for (e in 0 until numEdges) {
                maybeAddResult(shape, e)
            }
        }
    }

    private fun findClosestEdgesOptimized() {
        initQueue();
        // Repeatedly find the closest S2Cell to "target" and either split it into
        // its four children or process all of its edges.
        while (!queue.isEmpty()) {
            // We need to copy the top entry before removing it, and we need to
            // remove it before adding any new entries to the queue.
            val entry = queue.pollFirst()
            // Work around weird parse error in gcc 4.9 by using a local variable for
            // entry.distance.
            val distance = entry.distance
            if (!(distance < distanceLimit)) {
                queue.clear()  // Clear any remaining entries.
                break
            }
            // If this is already known to be an index cell, just process it.
            if (entry.indexCell != null) {
                processEdges(entry)
                continue
            }
            // Otherwise split the cell into its four children.  Before adding a
            // child back to the queue, we first check whether it is empty.  We do
            // this in two seek operations rather than four by seeking to the key
            // between children 0 and 1 and to the key between children 2 and 3.
            val id = entry.id;
            iter.seek(id.child(1).rangeMin())
            if (!iter.done() && iter.id() <= id.child(1).rangeMax()) {
                processOrEnqueue(id.child(1))
            }
            if (iter.prev() && iter.id() >= id.rangeMin()) {
                processOrEnqueue(id.child(0))
            }
            iter.seek(id.child(3).rangeMin())
            if (!iter.done() && iter.id() <= id.rangeMax()) {
                processOrEnqueue(id.child(3));
            }
            if (iter.prev() && iter.id() >= id.child(2).rangeMin()) {
                processOrEnqueue(id.child(2));
            }
        }
    }

    private fun initQueue() {
        check(queue.isEmpty())
        if (indexCovering.isEmpty()) {
            // We delay iterator initialization until now to make queries on very
            // small indexes a bit faster (i.e., where brute force is used).
            iter = index.cellIterator(InitialPosition.UNPOSITIONED)
        }

        // Optimization: if the user is searching for just the closest edge, and the
        // center of the target's bounding cap happens to intersect an index cell,
        // then we try to limit the search region to a small disc by first
        // processing the edges in that cell.  This sets distance_limit_ based on
        // the closest edge in that cell, which we can then use to limit the search
        // area.  This means that the cell containing "target" will be processed
        // twice, but in general this is still faster.
        //
        // TODO(ericv): Even if the cap center is not contained, we could still
        // process one or both of the adjacent index cells in S2CellId order,
        // provided that those cells are closer than distance_limit_.
        val cap = target.getCapBound()
        if (cap.isEmpty) return  // Empty target.
        if (options.getMaxResults() == 1 && iter.locate(cap.center)) {
            processEdges(QueueEntry(distanceFactory.zero(), iter.id(), iter.cell()))
            // Skip the rest of the algorithm if we found an intersecting edge.
            if (distanceLimit == distanceFactory.zero()) return
        }
        if (indexCovering.isEmpty()) initCovering()
        if (distanceLimit == distanceFactory.infinity()) {
            // Start with the precomputed index covering.
            for (i in 0 until indexCovering.size) {
                processOrEnqueue(indexCovering[i], index_cells[i])
            }
        } else {
            // Compute a covering of the search disc and intersect it with the
            // precomputed index covering.
            val coverer = S2RegionCoverer()
            coverer.setMaxCells(4)
            val radius = cap.radius + distanceLimit.getChordAngleBound()
            val search_cap = S2Cap(cap.center, radius)
            coverer.getFastCovering(search_cap, maxDistanceCovering)
            S2CellUnion.getIntersection(indexCovering, maxDistanceCovering, initial_cells)

            // Now we need to clean up the initial cells to ensure that they all
            // contain at least one cell of the S2ShapeIndex.  (Some may not intersect
            // the index at all, while other may be descendants of an index cell.)
            var i = 0
            var j = 0
            while (i < initial_cells.size) {
                val id_i = initial_cells[i]
                // Find the top-level cell that contains this initial cell.
                while (indexCovering[j].rangeMax() < id_i) ++j
                val id_j = indexCovering[j]
                if (id_i == id_j) {
                    // This initial cell is one of the top-level cells.  Use the
                    // precomputed S2ShapeIndexCell pointer to avoid an index seek.
                    processOrEnqueue(id_j, index_cells[j])
                    ++i; ++j
                } else {
                    // This initial cell is a proper descendant of a top-level cell.
                    // Check how it is related to the cells of the S2ShapeIndex.
                    val r = iter.locate(id_i)
                    if (r == CellRelation.INDEXED) {
                        // This cell is a descendant of an index cell.  Enqueue it and skip
                        // any other initial cells that are also descendants of this cell.
                        processOrEnqueue(iter.id(), iter.cell())
                        val last_id = iter.id().rangeMax()
                        while (++i < initial_cells.size && initial_cells[i] <= last_id)
                            continue;
                    } else {
                        // Enqueue the cell only if it contains at least one index cell.
                        if (r == CellRelation.SUBDIVIDED) processOrEnqueue(id_i, null)
                        ++i;
                    }
                }
            }
        }
    }

    private fun initCovering() {
        // Find the range of S2Cells spanned by the index and choose a level such
        // that the entire index can be covered with just a few cells.  These are
        // the "top-level" cells.  There are two cases:
        //
        //  - If the index spans more than one face, then there is one top-level cell
        // per spanned face, just big enough to cover the index cells on that face.
        //
        //  - If the index spans only one face, then we find the smallest cell "C"
        // that covers the index cells on that face (just like the case above).
        // Then for each of the 4 children of "C", if the child contains any index
        // cells then we create a top-level cell that is big enough to just fit
        // those index cells (i.e., shrinking the child as much as possible to fit
        // its contents).  This essentially replicates what would happen if we
        // started with "C" as the top-level cell, since "C" would immediately be
        // split, except that we take the time to prune the children further since
        // this will save work on every subsequent query.

        // Don't need to reserve index_cells_ since it is an InlinedVector.
        indexCovering.ensureCapacity(6)

        // TODO(ericv): Use a single iterator (iter_) below and save position
        // information using pair<S2CellId, const S2ShapeIndexCell*> type.
        val next = index.cellIterator(InitialPosition.BEGIN)
        val last = index.cellIterator(InitialPosition.END)
        last.prev();
        if (next.id() != last.id()) {
            // The index has at least two cells.  Choose a level such that the entire
            // index can be spanned with at most 6 cells (if the index spans multiple
            // faces) or 4 cells (it the index spans a single face).
            val level = next.id().getCommonAncestorLevel(last.id()) + 1

            // Visit each potential top-level cell except the last (handled below).
            val last_id = last.id().parent(level)
            var id = next.id().parent(level)
            while (id != last_id) {
                // Skip any top-level cells that don't contain any index cells.
                if (id.rangeMax() < next.id()) continue

                // Find the range of index cells contained by this top-level cell and
                // then shrink the cell if necessary so that it just covers them.
                val cell_first = next.clone()
                next.seek(id.rangeMax().next())
                val cell_last = next.clone()
                cell_last.prev()
                addInitialRange(cell_first, cell_last);
                id = id.next()
            }
        }
        addInitialRange(next, last);
    }

    // Add an entry to index_covering_ and index_cells_ that covers the given
    // inclusive range of cells.
    //
    // REQUIRES: "first" and "last" have a common ancestor.
    private fun addInitialRange(first: S2ShapeIndexCellIterator, last: S2ShapeIndexCellIterator) {
        if (first.id() == last.id()) {
            // The range consists of a single index cell.
            indexCovering.add(first.id())
            index_cells.add(first.cell())
        } else {
            // Add the lowest common ancestor of the given range.
            val level = first.id().getCommonAncestorLevel(last.id())
            check(level >= 0)
            indexCovering.add(first.id().parent(level))
            index_cells.add(null)
        }
    }

    private fun maybeAddResult(shape: S2Shape, edge_id: Int) {
        if (avoidDuplicates && !testedEdges.add(ShapeEdgeId(shape.id, edge_id))) {
            return
        }
        val edge = shape.edge(edge_id)
        val distance = distanceFactory.distance(distanceLimit)
        if (target.updateMinDistance(edge.v0, edge.v1, distance)) {
            addResult(Result(distance, shape.id, edge_id))
        }
    }

    private fun addResult(result: Result<T>) {
        if (options.getMaxResults() == 1) {
            // Optimization for the common case where only the closest edge is wanted.
            resultSingleton = result;
            distanceLimit = result.distance - options.maxError
        } else if (options.getMaxResults() == Options.kMaxMaxResults) {
            resultVector.add(result);  // Sort/unique at end.
        } else {
            // Add this edge to result_set_.  Note that even if we already have enough
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

    private fun processEdges(entry: QueueEntry<T>) {
        val indexCell = entry.indexCell!!
        for (s in 0 until indexCell.numClipped) {
            val clipped = indexCell.clipped(s)
            val shape = index.shape(clipped.shapeId)!!
            for (j in 0 until clipped.numEdges) {
                maybeAddResult(shape, clipped.edge(j))
            }
        }
    }

    // Enqueue the given cell id.
    // REQUIRES: iter_ is positioned at a cell contained by "id".
    private fun processOrEnqueue(id: S2CellId) {
        check(id.contains(iter.id()))
        if (iter.id() == id) {
            processOrEnqueue(id, iter.cell())
        } else {
            processOrEnqueue(id, null)
        }
    }

    // Add the given cell id to the queue.  "index_cell" is the corresponding
    // S2ShapeIndexCell, or nullptr if "id" is not an index cell.
    //
    // This version is called directly only by InitQueue().
    private fun processOrEnqueue(id: S2CellId, index_cell: S2ShapeIndexCell?) {
        if (index_cell != null) {
            // If this index cell has only a few edges, then it is faster to check
            // them directly rather than computing the minimum distance to the S2Cell
            // and inserting it into the queue.
            val kMinEdgesToEnqueue = 10
            val num_edges = countEdges(index_cell)
            if (num_edges == 0) return
            if (num_edges < kMinEdgesToEnqueue) {
                // Set "distance" to zero to avoid the expense of computing it.
                processEdges(QueueEntry(distanceFactory.zero(), id, index_cell))
                return
            }
        }
        // Otherwise compute the minimum distance to any point in the cell and add
        // it to the priority queue.
        val cell = S2Cell(id)
        var distance = distanceFactory.distance(distanceLimit)
        if (!target.updateMinDistance(cell, distance)) return
        if (useConservativeCellDistance) {
            // Ensure that "distance" is a lower bound on the true distance to the cell.
            distance -= options.maxError;  // operator-=() not defined.
        }
        queue.push(QueueEntry(distance, id, index_cell))
    }

}

