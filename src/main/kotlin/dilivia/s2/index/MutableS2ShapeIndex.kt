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

import dilivia.s2.Assertions
import dilivia.s2.MutableR2Rect
import dilivia.s2.R1Interval
import dilivia.s2.R2Rect
import dilivia.s2.S2CellId
import dilivia.s2.S2CellMetrics
import dilivia.s2.S2EdgeClipping
import dilivia.s2.S2PaddedCell
import dilivia.s2.coords.S2Coords
import dilivia.s2.math.R2Point
import dilivia.s2.region.S2CellUnion
import dilivia.s2.shape.Edge
import dilivia.s2.shape.InteriorTracker
import dilivia.s2.shape.S2ClippedShape
import dilivia.s2.shape.S2Shape
import mu.KotlinLogging
import java.util.*
import java.util.concurrent.atomic.AtomicReference
import java.util.concurrent.locks.ReentrantLock
import kotlin.math.abs
import kotlin.math.max
import kotlin.math.min
import kotlin.math.pow

typealias CellMap = TreeMap<S2CellId, S2ShapeIndexCell>

// MutableS2ShapeIndex is a class for in-memory indexing of polygonal geometry.
// The objects in the index are known as "shapes", and may consist of points,
// polylines, and/or polygons, possibly overlapping.  The index makes it very
// fast to answer queries such as finding nearby shapes, measuring distances,
// testing for intersection and containment, etc.
//
// MutableS2ShapeIndex allows not only building an index, but also updating it
// incrementally by adding or removing shapes (hence its name).  It is one of
// several implementations of the S2ShapeIndex interface.  MutableS2ShapeIndex
// is designed to be compact; usually it is smaller than the underlying
// geometry being indexed.  It is capable of indexing up to hundreds of
// millions of edges.  The index is also fast to construct.
//
// There are a number of built-in classes that work with S2ShapeIndex objects.
// Generally these classes accept any collection of geometry that can be
// represented by an S2ShapeIndex, i.e. any combination of points, polylines,
// and polygons.  Such classes include:
//
// - S2ContainsPointQuery: returns the shape(s) that contain a given point.
//
// - S2ClosestEdgeQuery: returns the closest edge(s) to a given point, edge,
//                       S2CellId, or S2ShapeIndex.
//
// - S2CrossingEdgeQuery: returns the edge(s) that cross a given edge.
//
// - S2BooleanOperation: computes boolean operations such as union,
//                       and boolean predicates such as containment.
//
// - S2ShapeIndexRegion: computes approximations for a collection of geometry.
//
// - S2ShapeIndexBufferedRegion: computes approximations that have been
//                               expanded by a given radius.
//
// Here is an example showing how to build an index for a set of polygons, and
// then then determine which polygon(s) contain each of a set of query points:
//
//   void TestContainment(const vector<S2Point>& points,
//                        const vector<S2Polygon*>& polygons) {
//     MutableS2ShapeIndex index
//     for (auto polygon : polygons) {
//       index.Add(absl::make_unique<S2Polygon::Shape>(polygon))
//     }
//     auto query = MakeS2ContainsPointQuery(&index)
//     for (const auto& point : points) {
//       for (S2Shape* shape : query.GetContainingShapes(point)) {
//         S2Polygon* polygon = polygons[shape->id()]
//         ... do something with (point, polygon) ...
//       }
//     }
//   }
//
// This example uses S2Polygon::Shape, which is one example of an S2Shape
// object.  S2Polyline and S2Loop also have nested Shape classes, and there are
// additional S2Shape types defined in *_shape.h.
//
// Internally, MutableS2ShapeIndex is essentially a map from S2CellIds to the
// set of shapes that intersect each S2CellId.  It is adaptively refined to
// ensure that no cell contains more than a small number of edges.
//
// For efficiency, updates are batched together and applied lazily on the
// first subsequent query.  Locking is used to ensure that MutableS2ShapeIndex
// has the same thread-safety properties as "vector": const methods are
// thread-safe, while non-const methods are not thread-safe.  This means that
// if one thread updates the index, you must ensure that no other thread is
// reading or updating the index at the same time.
//
// TODO(ericv): MutableS2ShapeIndex has an Encode() method that allows the
// index to be serialized.  An encoded S2ShapeIndex can be decoded either into
// its original form (MutableS2ShapeIndex) or into an EncodedS2ShapeIndex.
// The key property of EncodedS2ShapeIndex is that it can be constructed
// instantaneously, since the index is kept in its original encoded form.
// Data is decoded only when an operation needs it.  For example, to determine
// which shapes(s) contain a given query point only requires decoding the data
// in the S2ShapeIndexCell that contains that point.
// @constructor
// Create a MutableS2ShapeIndex with the given options.
// Creates a MutableS2ShapeIndex that uses the default option settings.
// Option values may be changed by calling Init().
// @property otions
// The options supplied for this index.
class MutableS2ShapeIndex(val options: Options = Options()) : S2ShapeIndex() {

    // The representation of an edge that has been queued for removal.
    data class RemovedShape(val shapeId: Int, val hasInterior: Boolean, val contains_tracker_origin: Boolean, val edges: List<Edge>)

    enum class IndexStatus {
        STALE,     // There are pending updates.
        UPDATING,  // Updates are currently being applied.
        FRESH,     // There are no pending updates.
    }

    // Options that affect construction of the MutableS2ShapeIndex.
    // @property maxEdgesPerCell
    // The maximum number of edges per cell.  If a cell has more than this
    // many edges that are not considered "long" relative to the cell size,
    // then it is subdivided.  (Whether an edge is considered "long" is
    // controlled by --s2shape_index_cell_size_to_long_edge_ratio flag.)
    //
    // Values between 10 and 50 represent a reasonable balance between memory
    // usage, construction time, and query time.  Small values make queries
    // faster, while large values make construction faster and use less memory.
    // Values higher than 50 do not save significant additional memory, and
    // query times can increase substantially, especially for algorithms that
    // visit all pairs of potentially intersecting edges (such as polygon
    // validation), since this is quadratic in the number of edges per cell.
    //
    // Note that the *average* number of edges per cell is generally slightly
    // less than half of the maximum value defined here.
    //
    // Defaults to value given by --s2shape_index_default_max_edges_per_cell.
    data class Options(val maxEdgesPerCell: Int = 10)

    // FaceEdge and ClippedEdge store temporary edge data while the index is being
    // updated.  FaceEdge represents an edge that has been projected onto a given
    // face, while ClippedEdge represents the portion of that edge that has been
    // clipped to a given S2Cell.
    //
    // While it would be possible to combine all the edge information into one
    // structure, there are two good reasons for separating it:
    //
    //  - Memory usage.  Separating the two classes means that we only need to
    //    store one copy of the per-face data no matter how many times an edge is
    //    subdivided, and it also lets us delay computing bounding boxes until
    //    they are needed for processing each face (when the dataset spans
    //    multiple faces).
    //
    //  - Performance.  UpdateEdges is significantly faster on large polygons when
    //    the data is separated, because it often only needs to access the data in
    //    ClippedEdge and this data is cached more successfully.

    data class FaceEdge(
            val shapeId: Int = -1,                  // The shape that this edge belongs to
            var edgeId: Int = 0,                    // Edge id within that shape
            var maxLevel: Int = 0,                  // Not desirable to subdivide this edge beyond this level
            val hasInterior: Boolean = false,       // Belongs to a shape of dimension 2.
            var a: R2Point = R2Point(),
            var b: R2Point = R2Point(),             // The edge endpoints, clipped to a given face
            var edge: Edge = Edge() // The edge endpoints
    )

    data class ClippedEdge(
            val faceEdge: FaceEdge = FaceEdge(),     // The original unclipped edge
            val bound: R2Rect = R2Rect()           // Bounding box for the clipped portion
    )

    // A BatchDescriptor represents a set of pending updates that will be applied
    // at the same time.  The batch consists of all updates with shape ids between
    // the current value of "ShapeIndex::pending_additions_begin_" (inclusive) and
    // "additions_end" (exclusive).  The first batch to be processed also
    // implicitly includes all shapes being removed.  "num_edges" is the total
    // number of edges that will be added or removed in this batch.
    data class BatchDescriptor(
            var additionsEnd: Int,
            val numEdges: Int
    )

    // UpdateState holds temporary data related to thread synchronization.  It
    // is only allocated while updates are being applied.
    class UpdateState(
            // This mutex is used as a condition variable.  It is locked by the
            // updating thread for the entire duration of the update; other threads
            // lock it in order to wait until the update is finished.
            val waitMutex: ReentrantLock = ReentrantLock(),

            // The number of threads currently waiting on "wait_mutex_".  The
            // UpdateState can only be freed when this number reaches zero.
            //
            // Reads and writes to this field are guarded by "lock_".
            var numWaiting: Int = 0) {

        fun destroy() {
            Assertions.assertEQ(0, numWaiting)
        }
    }

    private var updateState = UpdateState()

    // The shapes in the index, accessed by their shape id.  Removed shapes are
    // replaced by nullptr pointers.
    private val shapes: MutableList<S2Shape?> = mutableListOf()

    // A map from S2CellId to the set of clipped shapes that intersect that
    // cell.  The cell ids cover a set of non-overlapping regions on the
    // sphere.  Note that this field is updated lazily (see below).  Const
    // methods *must* call MaybeApplyUpdates() before accessing this field.
    // (The easiest way to achieve this is simply to use an Iterator.)
    private val cellMap = CellMap()

    // The id of the first shape that has been queued for addition but not
    // processed yet.
    private var pendingAdditionsBegin = 0

    // The set of shapes that have been queued for removal but not processed
    // yet.  Note that we need to copy the edge data since the caller is free to
    // destroy the shape once Release() has been called.  This field is present
    // only when there are removed shapes to process (to save memory).
    private val pendingRemovals = mutableListOf<RemovedShape>()

    // Reads and writes to this field are guarded by "lock_".
    private val indexStatus = AtomicReference(IndexStatus.FRESH)

    // Additions and removals are queued and processed on the first subsequent
    // query.  There are several reasons to do this:
    //
    //  - It is significantly more efficient to process updates in batches.
    //  - Often the index will never be queried, in which case we can save both
    //    the time and memory required to build it.  Examples:
    //     + S2Loops that are created simply to pass to an S2Polygon.  (We don't
    //       need the S2Loop index, because S2Polygon builds its own index.)
    //     + Applications that load a database of geometry and then query only
    //       a small fraction of it.
    //     + Applications that only read and write geometry (Decode/Encode).
    //
    // The main drawback is that we need to go to some extra work to ensure that
    // "const" methods are still thread-safe.  Note that the goal is *not* to
    // make this class thread-safe in general, but simply to hide the fact that
    // we defer some of the indexing work until query time.
    //
    // The textbook approach to this problem would be to use a mutex and a
    // condition variable.  Unfortunately pthread mutexes are huge (40 bytes).
    // Instead we use spinlock (which is only 4 bytes) to guard a few small
    // fields representing the current update status, and only create additional
    // state while the update is actually occurring.
    private val lock: ReentrantLock = ReentrantLock()

    // The number of distinct shape ids that have been assigned.  This equals
    // the number of shapes in the index provided that no shapes have ever been
    // removed.  (Shape ids are not reused.)
    override fun numShapeIds(): Int = shapes.size

    // Returns a pointer to the shape with the given id, or nullptr if the shape
    // has been removed from the index.
    override fun shape(id: Int): S2Shape? = shapes[id]

    // @constructor
    // Constructs an iterator positioned as specified.  By default iterators
    // are unpositioned, since this avoids an extra seek in this situation
    // where one of the seek methods (such as Locate) is immediately called.
    //
    // If you want to position the iterator at the beginning, e.g. in order to
    // loop through the entire index, do this instead:
    //
    //   for (MutableS2ShapeIndex::Iterator it(&index, S2ShapeIndex::BEGIN)
    //        !it.done(); it.Next()) { ... }
    inner class Iterator private constructor() : S2ShapeIndexIteratorBase() {

        private lateinit var end: S2CellId
        private lateinit var keySet: NavigableSet<S2CellId>
        private var currentCellId: S2CellId = S2CellId.sentinel()

        constructor(pos: InitialPosition = InitialPosition.UNPOSITIONED): this() {
            init(pos)
        }

        override fun clone(): Iterator {
            val clone = Iterator()
            clone.id = this.id()
            clone.cell.set(this.cell())
            clone.end = this.end
            clone.keySet = this.keySet
            clone.currentCellId = this.currentCellId
            return clone
        }

        // Initializes an iterator for the given MutableS2ShapeIndex.  This method
        // may also be called in order to restore an iterator to a valid state
        // after the underlying index has been updated (although it is usually
        // easier just to declare a new iterator whenever required, since iterator
        // construction is cheap).
        fun init(pos: InitialPosition) {
            logger.trace { "Iterator.init(pos = $pos)" }
            maybeApplyUpdates()
            initStale(pos)
        }

        // Initialize an iterator for the given MutableS2ShapeIndex without
        // applying any pending updates.  This can be used to observe the actual
        // current state of the index without modifying it in any way.
        fun initStale(pos: InitialPosition = InitialPosition.UNPOSITIONED) {
            logger.trace { "InitStale(pos = $pos): cellMap size = ${cellMap.size}" }
            end = if (cellMap.isEmpty()) S2CellId.sentinel() else cellMap.lastKey()
            keySet = cellMap.navigableKeySet()
            currentCellId = when {
                cellMap.isEmpty() -> S2CellId.sentinel()
                pos == InitialPosition.BEGIN -> keySet.first()
                else -> keySet.last()
            }

            logger.trace { "InitStale(pos = $pos): current = $currentCellId, end = $end" }
            refresh()
        }

        override fun begin() {
            // Make sure that the index has not been modified since Init() was called.
            Assertions.assert { isFresh() }
            currentCellId = keySet.first()
            refresh()
        }

        override fun finish() {
            currentCellId = S2CellId.sentinel()
            refresh()
        }

        override fun next() {
            Assertions.assert { !done() }
            currentCellId = keySet.higher(currentCellId) ?: S2CellId.sentinel()
            refresh()
        }

        override fun prev(): Boolean {
            if (keySet.isNotEmpty() && currentCellId == keySet.first()) return false
            currentCellId = keySet.lower(currentCellId) ?: keySet.firstOrNull() ?: S2CellId.sentinel()
            refresh()
            return true
        }

        override fun seek(target: S2CellId) {
            currentCellId = keySet.ceiling(target) ?: keySet.firstOrNull() ?: S2CellId.sentinel()
            refresh()
        }

        private fun refresh() {
            if (currentCellId == S2CellId.sentinel()) {
                setFinished()
            } else {
                setState(currentCellId, cellMap.getValue(currentCellId))
            }
        }

    }

    override fun iterator(pos: InitialPosition): Iterator = Iterator(pos)

    override fun begin(): ListIterator<S2Shape?> = shapes.listIterator(0)

    override fun end(): ListIterator<S2Shape?> = shapes.listIterator(shapes.lastIndex)

    // Takes ownership of the given shape and adds it to the index.  Also
    // assigns a unique id to the shape (shape->id()) and returns that id.
    // Shape ids are assigned sequentially starting from 0 in the order shapes
    // are added.  Invalidates all iterators and their associated data.
    fun add(shape: S2Shape): Int {
        // Additions are processed lazily by ApplyUpdates().
        shape.id = shapes.size
        shapes.add(shape)
        indexStatus.set(IndexStatus.STALE)
        return shape.id
    }

    // Removes the given shape from the index and return ownership to the caller.
    // Invalidates all iterators and their associated data.
    fun release(shapeId: Int): S2Shape {
        // This class updates itself lazily, because it is much more efficient to
        // process additions and removals in batches.  However this means that when
        // a shape is removed, we need to make a copy of all its edges, since the
        // client is free to delete "shape" once this call is finished.
        Assertions.assert { shapes[shapeId] != null }
        val shape = shapes.removeAt(shapeId)!!
        if (shapeId >= pendingAdditionsBegin) {
            // We are removing a shape that has not yet been added to the index,
            // so there is nothing else to do.
        } else {
            // We build the new RemovedShape in place, since it includes a potentially
            // large vector of edges that might be expensive to copy.
            val numEdges = shape.numEdges
            val removedEdges = mutableListOf<Edge>()
            for (e in 0 until numEdges) {
                removedEdges.add(shape.edge(e))
            }
            pendingRemovals.add(RemovedShape(
                    shapeId = shape.id,
                    hasInterior = shape.dimension == 2,
                    contains_tracker_origin = S2Shape.containsBruteForce(shape, kInteriorTrackerOrigin()),
                    edges = removedEdges
            ))
        }
        indexStatus.set(IndexStatus.STALE)
        return shape
    }

    // Resets the index to its original state and returns ownership of all
    // shapes to the caller.  This method is much more efficient than removing
    // all shapes one at a time.
    fun releaseAll(): List<S2Shape> {
        cellMap.clear()
        pendingAdditionsBegin = 0
        pendingRemovals.clear()
        //Assertions.assert { updateState == null }
        indexStatus.set(IndexStatus.FRESH)
        val removeShapes = shapes.toList().filterNotNull()
        shapes.clear()
        return removeShapes
    }

    // Resets the index to its original state and deletes all shapes.  Any
    // options specified via Init() are preserved.
    fun clear() {
        releaseAll()
    }

    // Calls to Add() and Release() are normally queued and processed on the
    // first subsequent query (in a thread-safe way).  This has many advantages,
    // the most important of which is that sometimes there *is* no subsequent
    // query, which lets us avoid building the index completely.
    //
    // This method forces any pending updates to be applied immediately.
    // Calling this method is rarely a good idea.  (One valid reason is to
    // exclude the cost of building the index from benchmark results.)
    fun forceBuild() {
        // No locks required because this is not a const method.  It is the client's
        // responsibility to ensure correct thread synchronization.
        if (indexStatus.get() != IndexStatus.FRESH) {
            applyUpdatesInternal()
            indexStatus.set(IndexStatus.FRESH)
        }
    }

    // Returns true if there are no pending updates that need to be applied.
    // This can be useful to avoid building the index unnecessarily, or for
    // choosing between two different algorithms depending on whether the index
    // is available.
    //
    // The returned index status may be slightly out of date if the index was
    // built in a different thread.  This is fine for the intended use (as an
    // efficiency hint), but it should not be used by internal methods  (see
    // MaybeApplyUpdates).
    fun isFresh(): Boolean = indexStatus.get() == IndexStatus.FRESH

    // Ensure that any pending updates have been applied.  This method must be
    // called before accessing the cell_map_ field, even if the index_status_
    // appears to be FRESH, because a memory barrier is required in order to
    // ensure that all the index updates are visible if the updates were done in
    // another thread.
    private fun maybeApplyUpdates() {
        // To avoid acquiring and releasing the spinlock on every query, we use
        // atomic operations when testing whether the status is FRESH and when
        // updating the status to be FRESH.  This guarantees that any thread that
        // sees a status of FRESH will also see the corresponding index updates.
        if (indexStatus.get() != IndexStatus.FRESH) {
            applyUpdatesThreadSafe()
        }
    }

    // Return true if this is the first update to the index.
    private fun isFirstUpdate(): Boolean {
        // Note that it is not sufficient to check whether cell_map_ is empty, since
        // entries are added during the update process.
        return pendingAdditionsBegin == 0
    }

    // Given that the given shape is being updated, return true if it is being
    // removed (as opposed to being added).
    private fun isShapeBeingRemoved(shape_id: Int): Boolean {
        // All shape ids being removed are less than all shape ids being added.
        return shape_id < pendingAdditionsBegin
    }

    // Apply any pending updates in a thread-safe way.
    private fun applyUpdatesThreadSafe() {
        logger.trace { "Iterator.applyUpdatesThreadSafe()" }
        lock.lock()
        when (indexStatus.get()) {
            IndexStatus.FRESH -> lock.unlock()
            IndexStatus.UPDATING -> {
                // Wait until the updating thread is finished.  We do this by attempting
                // to lock a mutex that is held by the updating thread.  When this mutex
                // is unlocked the index_status_ is guaranteed to be FRESH.
                ++updateState.numWaiting
                lock.unlock()
                updateState.waitMutex.lock()
                lock.lock()
                --updateState.numWaiting
                unlockAndSignal()  // Notify other waiting threads.
            }
            else -> {
                Assertions.assertEQ(IndexStatus.STALE, indexStatus.get())
                indexStatus.set(IndexStatus.UPDATING)
                // Allocate the extra state needed for thread synchronization.  We keep
                // the spinlock held while doing this, because (1) memory allocation is
                // fast, so the chance of a context switch while holding the lock is low
                // (2) by far the most common situation is that there is no contention,
                // and this saves an extra lock and unlock step; (3) even in the rare case
                // where there is contention, the main side effect is that some other
                // thread will burn a few CPU cycles rather than sleeping.
                updateState = UpdateState()
                // lock_.Lock wait_mutex *before* calling Unlock() to ensure that all other
                // threads will block on it.
                updateState.waitMutex.lock()
                // Release the spinlock before doing any real work.
                lock.unlock()
                applyUpdatesInternal()
                lock.lock()
                // index_status_ can be updated to FRESH only while locked *and* using
                // an atomic store operation, so that MaybeApplyUpdates() can check
                // whether the index is FRESH without acquiring the spinlock.
                indexStatus.set(IndexStatus.FRESH)
                unlockAndSignal()  // Notify any waiting threads.
            }
        }
    }

    // This method updates the index by applying all pending additions and
    // removals.  It does *not* update index_status_ (see ApplyUpdatesThreadSafe).
    private fun applyUpdatesInternal() {
        logger.trace { "--> Iterator.applyUpdatesInternal()" }
        // Check whether we have so many edges to process that we should process
        // them in multiple batches to save memory.  Building the index can use up
        // to 20x as much memory (per edge) as the final index size.
        val batches = getUpdateBatches()
        logger.trace { "${batches.size} batches to process." }
        for ((i, batch) in batches.withIndex()) {
            val allEdges: Array<MutableList<FaceEdge>> = Array(6) { mutableListOf() }
            logger.trace { "Batch $i : shape_limit=${batch.additionsEnd}, edges=${batch.numEdges}" }

            reserveSpace(batch, allEdges)
            val tracker = InteriorTracker()
            if (pendingRemovals.isNotEmpty()) {
                logger.trace { "Execute ${pendingRemovals.size} pending removals." }
                // The first batch implicitly includes all shapes being removed.
                for (pendingRemoval in pendingRemovals) {
                    removeShape(pendingRemoval, allEdges, tracker)
                }
                pendingRemovals.clear()
            }
            for (id in pendingAdditionsBegin until batch.additionsEnd) {
                addShape(id, allEdges, tracker)
            }
            for (face in 0..5) {
                updateFaceEdges(face, allEdges[face], tracker)
                // Save memory by clearing vectors after we are done with them.
                //vector<FaceEdge>().swap(all_edges[face])
            }
            pendingAdditionsBegin = batch.additionsEnd
        }
        // It is the caller's responsibility to update index_status_.
        logger.trace { "<-- Iterator.applyUpdatesInternal()" }
    }

    // Count the number of edges being updated, and break them into several
    // batches if necessary to reduce the amount of memory needed.  (See the
    // documentation for FLAGS_s2shape_index_tmp_memory_budget_mb.)
    private fun getUpdateBatches(): List<BatchDescriptor> {
        val batches = mutableListOf<BatchDescriptor>()
        // Count the edges being removed and added.
        var numEdgesRemoved = 0
        if (pendingRemovals.isNotEmpty()) {
            for (pending_removal in pendingRemovals) {
                numEdgesRemoved += pending_removal.edges.size
            }
        }
        var numEdgesAdded = 0
        for (id in pendingAdditionsBegin until shapes.size) {
            val shape = shape(id) ?: continue
            numEdgesAdded += shape.numEdges
        }
        var numEdges = numEdgesRemoved + numEdgesAdded

        if (numEdges * kTmpBytesPerEdge <= kTmpMemoryBudgetBytes) {
            // We can update all edges at once without exceeding kTmpMemoryBudgetBytes.
            batches.add(BatchDescriptor(shapes.size, numEdges))
            return batches
        }
        // Otherwise, break the updates into up to several batches, where the size
        // of each batch is chosen so that all batches use approximately the same
        // high-water memory.  GetBatchSizes() returns the recommended number of
        // edges in each batch.
        val batchSizes = getBatchSizes(
                numEdges,
                kMaxUpdateBatches,
                kFinalBytesPerEdge.toDouble(),
                kTmpBytesPerEdge.toDouble(),
                kTmpMemoryBudgetBytes.toDouble()
        )

        // We always process removed edges in a single batch, since (1) they already
        // take up a lot of memory because we have copied all their edges, and (2)
        // AbsorbIndexCell() uses (shapes_[id] == nullptr) to detect when a shape is
        // being removed, so in order to split the removals into batches we would
        // need a different approach (e.g., temporarily add fake entries to shapes_
        // and restore them back to nullptr as shapes are actually removed).
        numEdges = 0
        if (pendingRemovals.isNotEmpty()) {
            numEdges += numEdgesRemoved
            if (numEdges >= batchSizes[0]) {
                batches.add(BatchDescriptor(pendingAdditionsBegin, numEdges))
                numEdges = 0
            }
        }
        // Keep adding shapes to each batch until the recommended number of edges
        // for that batch is reached, then move on to the next batch.
        for (id in pendingAdditionsBegin until shapes.size) {
            val shape = shape(id) ?: continue
            numEdges += shape.numEdges
            if (numEdges >= batchSizes[batches.size]) {
                batches.add(BatchDescriptor(id + 1, numEdges))
                numEdges = 0
            }
        }
        // Some shapes have no edges.  If a shape with no edges is the last shape to
        // be added or removed, then the final batch may not include it, so we fix
        // that problem here.
        batches.last().additionsEnd = shapes.size
        Assertions.assertLE(batches.size, kMaxUpdateBatches)
        return batches
    }

    // Given "num_items" items, each of which uses "tmp_bytes_per_item" while it
    // is being updated but only "final_bytes_per_item" in the end, divide the
    // items into batches that have approximately the same *total* memory usage
    // consisting of the temporary memory needed for the items in the current
    // batch plus the final size of all the items that have already been
    // processed.  Use the fewest number of batches (but never more than
    // "max_batches") such that the total memory usage does not exceed the
    // combined final size of all the items plus "tmp_memory_budget_bytes".
    /* static */
    private fun getBatchSizes(numItems: Int, maxBatches: Int, finalBytesPerItem: Double, tmpBytesPerItem: Double, tmpMemoryBudgetBytes: Double): List<Int> {
        val batchSizes = mutableListOf<Int>()
        var remainingItems = numItems
        // This code tries to fit all the data into the same memory space
        // ("total_budget_bytes") at every iteration.  The data consists of some
        // number of processed items (at "final_bytes_per_item" each), plus some
        // number being updated (at "tmp_bytes_per_item" each).  The space occupied
        // by the items being updated is the "free space".  At each iteration, the
        // free space is multiplied by (1 - final_bytes_per_item/tmp_bytes_per_item)
        // as the items are converted into their final form.
        val finalBytes = remainingItems * finalBytesPerItem
        val finalBytesRatio = finalBytesPerItem / tmpBytesPerItem
        val freeSpaceMultiplier = 1 - finalBytesRatio

        // The total memory budget is the greater of the final size plus the allowed
        // temporary memory, or the minimum amount of memory required to limit the
        // number of batches to "max_batches".
        val totalBudgetBytes = max(finalBytes + tmpMemoryBudgetBytes, finalBytes / (1 - freeSpaceMultiplier.pow(maxBatches.toDouble())))

        // "max_batch_items" is the number of items in the current batch.
        var maxBatchItems = totalBudgetBytes / tmpBytesPerItem
        var i = 0
        while (i + 1 < maxBatches && remainingItems > 0) {
            val batchItems = min(remainingItems, (maxBatchItems + 1).toInt())
            batchSizes.add(batchItems)
            remainingItems -= batchItems
            maxBatchItems *= freeSpaceMultiplier
            ++i
        }
        Assertions.assertLE(batchSizes.size, maxBatches)
        return batchSizes
    }

    fun reserveSpace(batch: BatchDescriptor, allEdges: Array<MutableList<FaceEdge>>) {}

    // Clip all edges of the given shape to the six cube faces, add the clipped
    // edges to "all_edges", and start tracking its interior if necessary.
    private fun addShape(id: Int, all_edges: Array<MutableList<FaceEdge>>, tracker: InteriorTracker) {
        val shape = shape(id) ?: return  // This shape has already been removed.
        // Construct a template for the edges to be added.
        val edge = FaceEdge(
                shapeId = id,
                hasInterior = shape.dimension == 2
        )
        if (edge.hasInterior) {
            tracker.addShape(id, S2Shape.containsBruteForce(shape, tracker.focus()))
        }
        val numEdges = shape.numEdges
        for (e in 0 until numEdges) {
            addFaceEdge(edge.copy(edgeId = e, edge = shape.edge(e), maxLevel = getEdgeMaxLevel(edge.edge)), all_edges)
        }
    }

    private fun removeShape(removed: RemovedShape, allEdges: Array<MutableList<FaceEdge>>, tracker: InteriorTracker) {
        val edge = FaceEdge(
                shapeId = removed.shapeId,
                edgeId = -1,  // Not used or needed for removed edges.
                hasInterior = removed.hasInterior,
        )
        if (edge.hasInterior) {
            tracker.addShape(edge.shapeId, removed.contains_tracker_origin)
        }
        for (removed_edge in removed.edges) {
            addFaceEdge(edge.copy(edge = removed_edge, maxLevel = getEdgeMaxLevel(edge.edge)), allEdges)
        }
    }

    private fun addFaceEdge(edge: FaceEdge, allEdges: Array<MutableList<FaceEdge>>) {
        // Fast path: both endpoints are on the same face, and are far enough from
        // the edge of the face that don't intersect any (padded) adjacent face.
        val aFace = S2Coords.getFace(edge.edge.v0)
        if (aFace == S2Coords.getFace(edge.edge.v1)) {
            edge.a = S2Coords.validFaceXYZtoUV(aFace, edge.edge.v0)
            edge.b = S2Coords.validFaceXYZtoUV(aFace, edge.edge.v1)
            val kMaxUV = 1 - kCellPadding
            if (abs(edge.a[0]) <= kMaxUV && abs(edge.a[1]) <= kMaxUV && abs(edge.b[0]) <= kMaxUV && abs(edge.b[1]) <= kMaxUV) {
                allEdges[aFace].add(edge)
                return
            }
        }
        // Otherwise we simply clip the edge to all six faces.
        for (face in 0..5) {
            val abClippedUV = S2EdgeClipping.clipToPaddedFace(edge.edge.v0, edge.edge.v1, face, kCellPadding)
            if (abClippedUV != null) {
                edge.a = abClippedUV.first
                edge.b = abClippedUV.second
                allEdges[face].add(edge)
            }
        }
    }

    // Return the first level at which the edge will *not* contribute towards
    // the decision to subdivide.
    private fun getEdgeMaxLevel(edge: Edge): Int {
        // Compute the maximum cell size for which this edge is considered "long".
        // The calculation does not need to be perfectly accurate, so we use Norm()
        // rather than Angle() for speed.
        val cellSize = ((edge.v0 - edge.v1).norm() * s2shape_index_cell_size_to_long_edge_ratio)
        // Now return the first level encountered during subdivision where the
        // average cell size is at most "cell_size".
        return S2CellMetrics.kAvgEdge.getLevelForMaxValue(cellSize)
    }

    // Given a face and a vector of edges that intersect that face, add or remove
    // all the edges from the index.  (An edge is added if shapes_[id] is not
    // nullptr, and removed otherwise.)
    private fun updateFaceEdges(face: Int, face_edges: List<FaceEdge>, tracker: InteriorTracker) {
        val numEdges = face_edges.size
        if (numEdges == 0 && tracker.shapeIds().isEmpty()) return

        // Create the initial ClippedEdge for each FaceEdge.  Additional clipped
        // edges are created when edges are split between child cells.  We create
        // two arrays, one containing the edge data and another containing pointers
        // to those edges, so that during the recursion we only need to copy
        // pointers in order to propagate an edge to the correct child.
        val clippedEdgeStorage = mutableListOf<ClippedEdge>()
        val clippedEdges = mutableListOf<ClippedEdge>()
        val bound = MutableR2Rect()
        for (e in 0 until numEdges) {
            val clipped = ClippedEdge(
                    faceEdge = face_edges[e],
                    bound = R2Rect.fromPointPair(face_edges[e].a, face_edges[e].b)
            )
            //clipped.face_edge = &face_edges[e]
            //clipped.bound = R2Rect::FromPointPair(face_edges[e].a, face_edges[e].b)
            clippedEdgeStorage.add(clipped)
            clippedEdges.add(clippedEdgeStorage.last())
            bound.addRect(clipped.bound)
        }
        // Construct the initial face cell containing all the edges, and then update
        // all the edges in the index recursively.
        //EdgeAllocator alloc
        val faceId = S2CellId.fromFace(face)
        var pcell = S2PaddedCell(faceId, kCellPadding)

        // "disjoint_from_index" means that the current cell being processed (and
        // all its descendants) are not already present in the index.
        val disjointFromIndex = isFirstUpdate()

        if (numEdges > 0) {
            val shrunkId = shrinkToFit(pcell, bound)
            if (shrunkId != pcell.id) {
                // All the edges are contained by some descendant of the face cell.  We
                // can save a lot of work by starting directly with that cell, but if we
                // are in the interior of at least one shape then we need to create
                // index entries for the cells we are skipping over.
                skipCellRange(faceId.rangeMin(), shrunkId.rangeMin(), tracker, disjointFromIndex)
                pcell = S2PaddedCell(shrunkId, kCellPadding)
                updateEdges(pcell, clippedEdges, tracker, disjointFromIndex)
                skipCellRange(shrunkId.rangeMax().next(), faceId.rangeMax().next(), tracker, disjointFromIndex)
                return
            }
        }
        // Otherwise (no edges, or no shrinking is possible), subdivide normally.
        updateEdges(pcell, clippedEdges, tracker, disjointFromIndex)
    }


    private fun shrinkToFit(pcell: S2PaddedCell, bound: R2Rect): S2CellId {
        var shrunkId = pcell.shrinkToFit(bound)
        if (!isFirstUpdate() && shrunkId != pcell.id) {
            // Don't shrink any smaller than the existing index cells, since we need
            // to combine the new edges with those cells.
            // Use InitStale() to avoid applying updated recursively.
            val iter = iterator()
            iter.initStale()
            val r = iter.locate(shrunkId)
            if (r == CellRelation.INDEXED) {
                shrunkId = iter.id()
            }
        }
        return shrunkId
    }

    // Skip over the cells in the given range, creating index cells if we are
    // currently in the interior of at least one shape.
    private fun skipCellRange(begin: S2CellId, end: S2CellId, tracker: InteriorTracker/*, EdgeAllocator* alloc*/, disjoint_from_index: Boolean) {
        // If we aren't in the interior of a shape, then skipping over cells is easy.
        if (tracker.shapeIds().isEmpty()) return

        // Otherwise generate the list of cell ids that we need to visit, and create
        // an index entry for each one.
        for (skipped_id in S2CellUnion.fromBeginEnd(begin, end)) {
            updateEdges(S2PaddedCell(skipped_id, kCellPadding), mutableListOf(), tracker, disjoint_from_index)
        }
    }

    // Given a cell and a set of ClippedEdges whose bounding boxes intersect that
    // cell, add or remove all the edges from the index.  Temporary space for
    // edges that need to be subdivided is allocated from the given EdgeAllocator.
    // "disjoint_from_index" is an optimization hint indicating that cell_map_
    // does not contain any entries that overlap the given cell.
    private fun updateEdges(pcell: S2PaddedCell, edges: MutableList<ClippedEdge>, tracker: InteriorTracker, disjointFromIndex: Boolean) {
        var disjointFromIndex = disjointFromIndex
        // Cases where an index cell is not needed should be detected before this.
        Assertions.assert { (edges.isNotEmpty() || tracker.shapeIds().isNotEmpty()) }

        // This function is recursive with a maximum recursion depth of 30
        // (S2CellId::kMaxLevel).  Note that using an explicit stack does not seem
        // to be any faster based on profiling.

        // Incremental updates are handled as follows.  All edges being added or
        // removed are combined together in "edges", and all shapes with interiors
        // are tracked using "tracker".  We subdivide recursively as usual until we
        // encounter an existing index cell.  At this point we "absorb" the index
        // cell as follows:
        //
        //   - Edges and shapes that are being removed are deleted from "edges" and
        //     "tracker".
        //   - All remaining edges and shapes from the index cell are added to
        //     "edges" and "tracker".
        //   - Continue subdividing recursively, creating new index cells as needed.
        //   - When the recursion gets back to the cell that was absorbed, we
        //     restore "edges" and "tracker" to their previous state.
        //
        // Note that the only reason that we include removed shapes in the recursive
        // subdivision process is so that we can find all of the index cells that
        // contain those shapes efficiently, without maintaining an explicit list of
        // index cells for each shape (which would be expensive in terms of memory).
        var indexCellAbsorbed = false
        if (!disjointFromIndex) {
            // There may be existing index cells contained inside "pcell".  If we
            // encounter such a cell, we need to combine the edges being updated with
            // the existing cell contents by "absorbing" the cell.
            // Use InitStale() to avoid applying updated recursively.
            val iter = iterator()
            iter.initStale()
            when (val r = iter.locate(pcell.id)) {
                CellRelation.DISJOINT -> disjointFromIndex = true
                CellRelation.INDEXED -> {
                    // Absorb the index cell by transferring its contents to "edges" and
                    // deleting it.  We also start tracking the interior of any new shapes.
                    absorbIndexCell(pcell, iter, edges, tracker)
                    indexCellAbsorbed = true
                    disjointFromIndex = true
                }
                else -> Assertions.assertEQ(CellRelation.SUBDIVIDED, r)
            }
        }

        // If there are existing index cells below us, then we need to keep
        // subdividing so that we can merge with those cells.  Otherwise,
        // MakeIndexCell checks if the number of edges is small enough, and creates
        // an index cell if possible (returning true when it does so).
        if (!disjointFromIndex || !makeIndexCell(pcell, edges, tracker)) {
            // Reserve space for the edges that will be passed to each child.  This is
            // important since otherwise the running time is dominated by the time
            // required to grow the vectors.  The amount of memory involved is
            // relatively small, so we simply reserve the maximum space for every child.
            val childEdges = Array(2) { Array(2) { mutableListOf<ClippedEdge>() } }  // [i][j]
            val numEdges = edges.size
            //for (i in 0..1) {
            //    for (j in 0..1) {
            //    child_edges[i][j].reserve(num_edges)
            //}
            //}

            // Remember the current size of the EdgeAllocator so that we can free any
            // edges that are allocated during edge splitting.
            //size_t alloc_size = alloc->size()

            // Compute the middle of the padded cell, defined as the rectangle in
            // (u,v)-space that belongs to all four (padded) children.  By comparing
            // against the four boundaries of "middle" we can determine which children
            // each edge needs to be propagated to.
            val middle = pcell.middle()

            // Build up a vector edges to be passed to each child cell.  The (i,j)
            // directions are left (i=0), right (i=1), lower (j=0), and upper (j=1).
            // Note that the vast majority of edges are propagated to a single child.
            // This case is very fast, consisting of between 2 and 4 floating-point
            // comparisons and copying one pointer.  (ClipVAxis is inline.)
            for (e in 0 until numEdges) {
                val edge = edges[e]
                when {
                    edge.bound[0].hi <= middle[0].lo -> {
                        // Edge is entirely contained in the two left children.
                        clipVAxis(edge, middle[1], childEdges[0])
                    }
                    edge.bound[0].lo >= middle[0].hi -> {
                        // Edge is entirely contained in the two right children.
                        clipVAxis(edge, middle[1], childEdges[1])
                    }
                    edge.bound[1].hi <= middle[1].lo -> {
                        // Edge is entirely contained in the two lower children.
                        childEdges[0][0].add(clipUBound(edge, 1, middle[0].hi))
                        childEdges[1][0].add(clipUBound(edge, 0, middle[0].lo))
                    }
                    edge.bound[1].lo >= middle[1].hi -> {
                        // Edge is entirely contained in the two upper children.
                        childEdges[0][1].add(clipUBound(edge, 1, middle[0].hi))
                        childEdges[1][1].add(clipUBound(edge, 0, middle[0].lo))
                    }
                    else -> {
                        // The edge bound spans all four children.  The edge itself intersects
                        // either three or four (padded) children.
                        val left = clipUBound(edge, 1, middle[0].hi)
                        clipVAxis(left, middle[1], childEdges[0])
                        val right = clipUBound(edge, 0, middle[0].lo)
                        clipVAxis(right, middle[1], childEdges[1])
                    }
                }
            }
            // Free any memory reserved for children that turned out to be empty.  This
            // step is cheap and reduces peak memory usage by about 10% when building
            // large indexes (> 10M edges).
            //for (int i = 0; i < 2; ++i) {
            //    for (int j = 0; j < 2; ++j) {
            //    if (child_edges[i][j].empty()) {
            //        vector<const ClippedEdge*>().swap(child_edges[i][j])
            //    }
            //}
            //}

            // Now recursively update the edges in each child.  We call the children in
            // increasing order of S2CellId so that when the index is first constructed,
            // all insertions into cell_map_ are at the end (which is much faster).
            for (pos in 0..3) {
                val (i, j) = pcell.getChildIJ(pos)
                if (childEdges[i][j].isNotEmpty() || tracker.shapeIds().isNotEmpty()) {
                    updateEdges(S2PaddedCell(pcell, i, j), childEdges[i][j], tracker, disjointFromIndex)
                }
            }
            // Free any temporary edges that were allocated during clipping.
            //alloc->Reset(alloc_size)
        }
        if (indexCellAbsorbed) {
            // Restore the state for any edges being removed that we are tracking.
            tracker.restoreStateBefore(pendingAdditionsBegin)
        }
    }

    // Given an edge and an interval "middle" along the v-axis, clip the edge
    // against the boundaries of "middle" and add the edge to the corresponding
    // children.
    /* static */
    private fun clipVAxis(edge: ClippedEdge, middle: R1Interval, child_edges: Array<MutableList<ClippedEdge>>) {
        when {
            edge.bound[1].hi <= middle.lo -> child_edges[0].add(edge) // Edge is entirely contained in the lower child.
            edge.bound[1].lo >= middle.hi -> child_edges[1].add(edge)  // Edge is entirely contained in the upper child.
            else -> {
                // The edge bound spans both children.
                child_edges[0].add(clipVBound(edge, 1, middle.hi))
                child_edges[1].add(clipVBound(edge, 0, middle.lo))
            }
        }
    }

    // Given an edge, clip the given endpoint (lo=0, hi=1) of the u-axis so that
    // it does not extend past the given value.
    /* static */
    private fun clipUBound(edge: ClippedEdge, u_end: Int, u: Double): ClippedEdge {
        // First check whether the edge actually requires any clipping.  (Sometimes
        // this method is called when clipping is not necessary, e.g. when one edge
        // endpoint is in the overlap area between two padded child cells.)
        if (u_end == 0) {
            if (edge.bound[0].lo >= u) return edge
        } else {
            if (edge.bound[0].hi <= u) return edge
        }
        // We interpolate the new v-value from the endpoints of the original edge.
        // This has two advantages: (1) we don't need to store the clipped endpoints
        // at all, just their bounding box; and (2) it avoids the accumulation of
        // roundoff errors due to repeated interpolations.  The result needs to be
        // clamped to ensure that it is in the appropriate range.
        val e = edge.faceEdge
        val v = edge.bound[1].project(S2EdgeClipping.interpolateDouble(u, e.a[0], e.b[0], e.a[1], e.b[1]))

        // Determine which endpoint of the v-axis bound to update.  If the edge
        // slope is positive we update the same endpoint, otherwise we update the
        // opposite endpoint.
        val vEnd = u_end xor (if ((e.a[0] > e.b[0]) != (e.a[1] > e.b[1])) 1 else 0)
        return updateBound(edge, u_end, u, vEnd, v)
    }

    // Given an edge, clip the given endpoint (lo=0, hi=1) of the v-axis so that
    // it does not extend past the given value.
    /* static */
    private fun clipVBound(edge: ClippedEdge, v_end: Int, v: Double): ClippedEdge {
        // See comments in ClipUBound.
        if (v_end == 0) {
            if (edge.bound[1].lo >= v) return edge
        } else {
            if (edge.bound[1].hi <= v) return edge
        }
        val e = edge.faceEdge
        val u = edge.bound[0].project(S2EdgeClipping.interpolateDouble(v, e.a[1], e.b[1], e.a[0], e.b[0]))
        val uEnd = v_end xor (if ((e.a[0] > e.b[0]) != (e.a[1] > e.b[1])) 1 else 0)
        return updateBound(edge, uEnd, u, v_end, v)
    }

    // Given an edge and two bound endpoints that need to be updated, allocate and
    // return a new edge with the updated bound.
    /* static */
    private fun updateBound(edge: ClippedEdge, u_end: Int, u: Double, v_end: Int, v: Double): ClippedEdge {
        logger.trace { "--> updateBound(edge = $edge, uEnd = $u_end, u = $u, vEnd = $v_end, v = $v)" }
        val bound = edge.bound.toMutable()
        bound[0][u_end] = u
        bound[1][v_end] = v
        bound[0][1 - u_end] = bound[0][1 - u_end]
        bound[1][1 - v_end] = bound[1][1 - v_end]
        logger.trace { "Bound = $bound" }
        val clipped = ClippedEdge(faceEdge = edge.faceEdge, bound = bound)
        logger.trace { "<-- updateBound() => $clipped" }
        Assertions.assert { !clipped.bound.isEmpty }
        Assertions.assert { edge.bound.contains(clipped.bound) }
        return clipped
    }

    // Absorb an index cell by transferring its contents to "edges" and/or
    // "tracker", and then delete this cell from the index.  If "edges" includes
    // any edges that are being removed, this method also updates their
    // InteriorTracker state to correspond to the exit vertex of this cell, and
    // saves the InteriorTracker state by calling SaveAndClearStateBefore().  It
    // is the caller's responsibility to restore this state by calling
    // RestoreStateBefore() when processing of this cell is finished.
    private fun absorbIndexCell(pcell: S2PaddedCell, iter: Iterator, edges: MutableList<ClippedEdge>, tracker: InteriorTracker) {
        Assertions.assertEQ(pcell.id, iter.id())

        // When we absorb a cell, we erase all the edges that are being removed.
        // However when we are finished with this cell, we want to restore the state
        // of those edges (since that is how we find all the index cells that need
        // to be updated).  The edges themselves are restored automatically when
        // UpdateEdges returns from its recursive call, but the InteriorTracker
        // state needs to be restored explicitly.
        //
        // Here we first update the InteriorTracker state for removed edges to
        // correspond to the exit vertex of this cell, and then save the
        // InteriorTracker state.  This state will be restored by UpdateEdges when
        // it is finished processing the contents of this cell.
        if (tracker.isActive() && edges.isNotEmpty() && isShapeBeingRemoved(edges[0].faceEdge.shapeId)) {
            // We probably need to update the InteriorTracker.  ("Probably" because
            // it's possible that all shapes being removed do not have interiors.)
            if (!tracker.atCellId(pcell.id)) {
                tracker.moveTo(pcell.getEntryVertex())
            }
            tracker.drawTo(pcell.getExitVertex())
            tracker.setNextCellId(pcell.id.next())
            for (edge in edges) {
                val faceEdge = edge.faceEdge
                if (!isShapeBeingRemoved(faceEdge.shapeId)) {
                    break  // All shapes being removed come first.
                }
                if (faceEdge.hasInterior) {
                    tracker.testEdge(faceEdge.shapeId, faceEdge.edge)
                }
            }
        }
        // Save the state of the edges being removed, so that it can be restored
        // when we are finished processing this cell and its children.  We don't
        // need to save the state of the edges being added because they aren't being
        // removed from "edges" and will therefore be updated normally as we visit
        // this cell and its children.
        tracker.saveAndClearStateBefore(pendingAdditionsBegin)

        // Create a FaceEdge for each edge in this cell that isn't being removed.
        val faceEdges = mutableListOf<FaceEdge>()
        var trackerMoved = false
        val cell = iter.cell() ?: S2ShapeIndexCell()
        for (s in 0 until cell.numClipped()) {
            val clipped = cell.clipped(s)
            val shapeId = clipped.shapeId
            val shape = shape(shapeId) ?: continue // This shape is being removed.
            val numEdges = clipped.numEdges()

            // If this shape has an interior, start tracking whether we are inside the
            // shape.  UpdateEdges() wants to know whether the entry vertex of this
            // cell is inside the shape, but we only know whether the center of the
            // cell is inside the shape, so we need to test all the edges against the
            // line segment from the cell center to the entry vertex.
            val edge = FaceEdge(
                    shapeId = shape.id,
                    hasInterior = (shape.dimension == 2)
            )
            if (edge.hasInterior) {
                tracker.addShape(shapeId, clipped.containsCenter)
                // There might not be any edges in this entire cell (i.e., it might be
                // in the interior of all shapes), so we delay updating the tracker
                // until we see the first edge.
                if (!trackerMoved && numEdges > 0) {
                    tracker.moveTo(pcell.getCenter())
                    tracker.drawTo(pcell.getEntryVertex())
                    tracker.setNextCellId(pcell.id)
                    trackerMoved = true
                }
            }
            for (i in 0 until numEdges) {
                val e = clipped.edge(i)
                edge.edgeId = e
                edge.edge = shape.edge(e)
                edge.maxLevel = getEdgeMaxLevel(edge.edge)
                if (edge.hasInterior) tracker.testEdge(shapeId, edge.edge)
                val clippedToFace = S2EdgeClipping.clipToPaddedFace(edge.edge.v0, edge.edge.v1, pcell.id.face(), kCellPadding)
                if (clippedToFace == null) {
                    logger.error { "Invariant failure in MutableS2ShapeIndex" }
                }
                faceEdges.add(edge)
            }
        }
        // Now create a ClippedEdge for each FaceEdge, and put them in "new_edges".
        val newEdges = mutableListOf<ClippedEdge>()
        for (face_edge in faceEdges) {
            val clipped = ClippedEdge(
                    faceEdge = face_edge,
                    bound = S2EdgeClipping.getClippedEdgeBound(face_edge.a, face_edge.b, pcell.bound)
            )
            newEdges.add(clipped)
        }
        // Discard any edges from "edges" that are being removed, and append the
        // remainder to "new_edges".  (This keeps the edges sorted by shape id.)
        for (i in edges.indices) {
            val clipped = edges[i]
            if (!isShapeBeingRemoved(clipped.faceEdge.shapeId)) {
                newEdges.addAll(edges.subList(i, edges.size))
                break
            }
        }
        // Update the edge list and delete this cell from the index.
        edges.clear()
        edges.addAll(newEdges)
        cellMap.remove(pcell.id)
        //delete &cell
    }

    // Attempt to build an index cell containing the given edges, and return true
    // if successful.  (Otherwise the edges should be subdivided further.)
    private fun makeIndexCell(pcell: S2PaddedCell, edges: List<ClippedEdge>, tracker: InteriorTracker): Boolean {
        logger.trace { "--> makeIndexCell(pcell = $pcell, edges = $edges)" }
        if (edges.isEmpty() && tracker.shapeIds().isEmpty()) {
            // No index cell is needed.  (In most cases this situation is detected
            // before we get to this point, but this can happen when all shapes in a
            // cell are removed.)
            logger.trace { "<-- makeIndexCell(pcell = $pcell, edges = $edges) = true : isEdgesEmpty = ${edges.isEmpty()}, isTrackerShapeIdsEmpty = ${tracker.shapeIds().isEmpty()}" }
            return true
        }

        // Count the number of edges that have not reached their maximum level yet.
        // Return false if there are too many such edges.
        var count = 0
        for (edge in edges) {
            count += if (pcell.level < edge.faceEdge.maxLevel) 1 else 0
            if (count > options.maxEdgesPerCell) {
                logger.trace { "<-- makeIndexCell(pcell = $pcell, edges = $edges) = false: edges count = $count > ${options.maxEdgesPerCell}" }
                return false
            }
        }
        logger.trace { "Edges that have not reached their maximum level = $count" }

        // Possible optimization: Continue subdividing as long as exactly one child
        // of "pcell" intersects the given edges.  This can be done by finding the
        // bounding box of all the edges and calling ShrinkToFit():
        //
        // S2CellId cellid = pcell.ShrinkToFit(GetRectBound(edges))
        //
        // Currently this is not beneficial; it slows down construction by 4-25%
        // (mainly computing the union of the bounding rectangles) and also slows
        // down queries (since more recursive clipping is required to get down to
        // the level of a spatial index cell).  But it may be worth trying again
        // once "contains_center" is computed and all algorithms are modified to
        // take advantage of it.

        // We update the InteriorTracker as follows.  For every S2Cell in the index
        // we construct two edges: one edge from entry vertex of the cell to its
        // center, and one from the cell center to its exit vertex.  Here "entry"
        // and "exit" refer the S2CellId ordering, i.e. the order in which points
        // are encountered along the S2 space-filling curve.  The exit vertex then
        // becomes the entry vertex for the next cell in the index, unless there are
        // one or more empty intervening cells, in which case the InteriorTracker
        // state is unchanged because the intervening cells have no edges.

        // Shift the InteriorTracker focus point to the center of the current cell.
        if (tracker.isActive() && edges.isNotEmpty()) {
            if (!tracker.atCellId(pcell.id)) {
                tracker.moveTo(pcell.getEntryVertex())
            }
            tracker.drawTo(pcell.getCenter())
            testAllEdges(edges, tracker)
        }
        // Allocate and fill a new index cell.  To get the total number of shapes we
        // need to merge the shapes associated with the intersecting edges together
        // with the shapes that happen to contain the cell center.
        val cshapeIds = tracker.shapeIds()
        val numShapes = countShapes(edges, cshapeIds)
        val cell = S2ShapeIndexCell()
        //val base = cell.add_shapes(num_shapes)

        // To fill the index cell we merge the two sources of shapes: "edge shapes"
        // (those that have at least one edge that intersects this cell), and
        // "containing shapes" (those that contain the cell center).  We keep track
        // of the index of the next intersecting edge and the next containing shape
        // as we go along.  Both sets of shape ids are already sorted.
        var enext = 0
        val cShapeIterator = cshapeIds.iterator()
        var cnext = if (cShapeIterator.hasNext()) cShapeIterator.next() else -1
        for (i in 0 until numShapes) {
            //S2ClippedShape* clipped = base + i;
            var eshape_id = numShapeIds()
            var cshape_id = eshape_id;  // Sentinels
            if (enext != edges.size) {
                eshape_id = edges[enext].faceEdge.shapeId
            }
            if (cnext != -1) {
                cshape_id = cnext
            }
            val ebegin = enext
            if (cshape_id < eshape_id) {
                // The entire cell is in the shape interior.
                cell.addClipped(S2ClippedShape(cshape_id, emptyList(), true))
                cnext = if (cShapeIterator.hasNext()) cShapeIterator.next() else -1
              //  ++cnext;
            } else {
                // Count the number of edges for this shape and allocate space for them.
                while (enext < edges.size && edges[enext].faceEdge.shapeId == eshape_id) {
                    ++enext;
                }

                val clippedShapeEdges = mutableListOf<Int>()
                for (e in ebegin until enext) {
                    clippedShapeEdges.add(edges[e].faceEdge.edgeId)
                }
                var containsCenter = false
                if (cshape_id == eshape_id) {
                    containsCenter = true
                    cnext = if (cShapeIterator.hasNext()) cShapeIterator.next() else -1
                }
                cell.addClipped(S2ClippedShape(eshape_id, clippedShapeEdges, containsCenter))
            }
        }

        // UpdateEdges() visits cells in increasing order of S2CellId, so during
        // initial construction of the index all insertions happen at the end.  It
        // is much faster to give an insertion hint in this case.  Otherwise the
        // hint doesn't do much harm.  With more effort we could provide a hint even
        // during incremental updates, but this is probably not worth the effort.
        logger.trace { "Add cell ${pcell.id} => $cell" }
        cellMap[pcell.id] = cell

        // Shift the InteriorTracker focus point to the exit vertex of this cell.
        if (tracker.isActive() && edges.isNotEmpty()) {
            tracker.drawTo(pcell.getExitVertex())
            testAllEdges(edges, tracker)
            tracker.setNextCellId(pcell.id.next())
        }

        logger.trace { "<-- makeIndexCell(pcell = $pcell, edges = $edges) = true" }
        return true
    }

    // Call tracker->TestEdge() on all edges from shapes that have interiors.
/* static */
    fun testAllEdges(edges: List<ClippedEdge>, tracker: InteriorTracker) {
        for (edge in edges) {
            val faceEdge = edge.faceEdge
            if (faceEdge.hasInterior) {
                tracker.testEdge(faceEdge.shapeId, faceEdge.edge)
            }
        }
    }

    // Return the number of distinct shapes that are either associated with the
// given edges, or that are currently stored in the InteriorTracker.
/* static */
    fun countShapes(edges: List<ClippedEdge>, cshape_ids: List<Int>): Int {
        logger.trace { "--> countShapes(edges = $edges, shapeIds = $cshape_ids)" }

        var count = 0
        var last_shape_id = -1
        val shapeIdsIterator = cshape_ids.iterator()
        //ShapeIdSet::const_iterator cnext = cshape_ids.begin();  // Next shape
        for (edge in edges) {
            logger.trace { "Process edge: $edge" }
            if (edge.faceEdge.shapeId != last_shape_id) {
                ++count
                logger.trace { "Edge is own by a new shape: edge shape = ${edge.faceEdge.shapeId}, last shape id = $last_shape_id => count = $count" }
                last_shape_id = edge.faceEdge.shapeId
                // Skip over any containing shapes up to and including this one,
                // updating "count" appropriately.
                logger.trace { "Skip over any containing shapes up to and including this one." }
                while (shapeIdsIterator.hasNext()) {
                    val cnext = shapeIdsIterator.next()
                    logger.trace { "next shape id = $cnext" }
                    if (cnext > last_shape_id) {
                        logger.trace { "It is greater, break." }
                        break
                    }
                    if (cnext < last_shape_id) {
                        ++count
                        logger.trace { "It is smaller, count shape => count = $count" }
                    }
                }
            }
        }
        // Count any remaining containing shapes.
        count += shapeIdsIterator.asSequence().count()
        logger.trace { "<-- countShapes = $count" }
        return count
    }


/*
    protected:
    std::unique_ptr<IteratorBase> NewIterator(InitialPosition pos) const override

    // Internal methods are documented with their definitions.
    static int CountShapes(const std::vector<const ClippedEdge*>& edges,const ShapeIdSet& cshape_ids)
    bool MakeIndexCell(const S2PaddedCell& pcell,const std::vector<const ClippedEdge*>& edges,InteriorTracker* tracker)
    static void TestAllEdges(const std::vector<const ClippedEdge*>& edges,InteriorTracker* tracker)
*/


    // Releases lock_ and wakes up any waiting threads by releasing wait_mutex.
// If this was the last waiting thread, also deletes update_state_.
// REQUIRES: lock_ is held.
// REQUIRES: wait_mutex is held.
    private fun unlockAndSignal() {
        Assertions.assertEQ(IndexStatus.FRESH, indexStatus.get())
        //val num_waiting = updateState.numWaiting
        lock.unlock()
        // Allow another waiting thread to proceed.  Note that no new threads can
        // start waiting because the index_status_ is now FRESH, and the caller is
        // required to prevent any new mutations from occurring while these const
        // methods are running.
        //
        // We need to unlock wait_mutex before destroying it even if there are no
        // waiting threads.
        updateState.waitMutex.unlock()
        //if (num_waiting == 0) {
        //updateState.reset()
        //}
    }

    /*
        UNLOCK_FUNCTION(lock_)
        UNLOCK_FUNCTION(update_state_->wait_mutex)
    */
    companion object {

        private val logger = KotlinLogging.logger(MutableS2ShapeIndex::class.java.name)

        // FLAGS_s2shape_index_cell_size_to_long_edge_ratio
        //
        // The cell size relative to the length of an edge at which it is first
        // considered to be "long".  Long edges do not contribute toward the decision
        // to subdivide a cell further.  For example, a value of 2.0 means that the
        // cell must be at least twice the size of the edge in order for that edge to
        // be counted.  There are two reasons for not counting long edges: (1) such
        // edges typically need to be propagated to several children, which increases
        // time and memory costs without much benefit, and (2) in pathological cases,
        // many long edges close together could force subdivision to continue all the
        // way to the leaf cell level.
        var s2shape_index_cell_size_to_long_edge_ratio = 1.0

        // The amount by which cells are "padded" to compensate for numerical errors
        // when clipping line segments to cell boundaries.
        // The total error when clipping an edge comes from two sources:
        // (1) Clipping the original spherical edge to a cube face (the "face edge").
        //     The maximum error in this step is S2::kFaceClipErrorUVCoord.
        // (2) Clipping the face edge to the u- or v-coordinate of a cell boundary.
        //     The maximum error in this step is S2::kEdgeClipErrorUVCoord.
        // Finally, since we encounter the same errors when clipping query edges, we
        // double the total error so that we only need to pad edges during indexing
        // and not at query time.
        val kCellPadding: Double = 2 * (S2EdgeClipping.kFaceClipErrorUVCoord + S2EdgeClipping.kEdgeClipErrorUVCoord)

        // Defines the initial focus point of MutableS2ShapeIndex::InteriorTracker
        // (the start of the S2CellId space-filling curve).
        //
        // TODO(ericv): Move InteriorTracker here to avoid the need for this method.
        fun kInteriorTrackerOrigin() = S2Coords.faceUVtoXYZ(0, -1.0, -1.0).normalize()

        // TODO(fmeurisse) check memory estimation with kotlin objects.
        // The following memory estimates are based on heap profiling.
        //
        // The final size of a MutableS2ShapeIndex depends mainly on how finely the
        // index is subdivided, as controlled by Options::max_edges_per_cell() and
        // --s2shape_index_default_max_edges_per_cell. For realistic values of
        // max_edges_per_cell() and shapes with moderate numbers of edges, it is
        // difficult to get much below 8 bytes per edge.  [The minimum possible size
        // is 4 bytes per edge (to store a 32-bit edge id in an S2ClippedShape) plus
        // 24 bytes per shape (for the S2ClippedShape itself plus a pointer in the
        // shapes_ vector.]
        //
        // The temporary memory consists mainly of the FaceEdge and ClippedEdge
        // structures plus a ClippedEdge pointer for every level of recursive
        // subdivision.  For very large indexes this can be 200 bytes per edge.
        const val kFinalBytesPerEdge = 8
        const val kTmpBytesPerEdge = 200
        const val kTmpMemoryBudgetBytes = 100 shl 20// static_cast<size_t>(FLAGS_s2shape_index_tmp_memory_budget_mb) << 20

        // We arbitrarily limit the number of batches just as a safety measure.
        // With the current default memory budget of 100 MB, this limit is not
        // reached even when building an index of 350 million edges.
        const val kMaxUpdateBatches = 100

    }
}
