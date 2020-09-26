package dilivia.s2.shape

import dilivia.s2.*
import dilivia.s2.math.R2Point
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
//     MutableS2ShapeIndex index;
//     for (auto polygon : polygons) {
//       index.Add(absl::make_unique<S2Polygon::Shape>(polygon));
//     }
//     auto query = MakeS2ContainsPointQuery(&index);
//     for (const auto& point : points) {
//       for (S2Shape* shape : query.GetContainingShapes(point)) {
//         S2Polygon* polygon = polygons[shape->id()];
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
    data class RemovedShape(val shapeId: Int, val hasInterior: Boolean, val contains_tracker_origin: Boolean, val edges: List<S2Shape.Edge>)

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
            val shapeId: Int,                       // The shape that this edge belongs to
            val edgeId: Int = 0,                    // Edge id within that shape
            val maxLevel: Int = 0,                  // Not desirable to subdivide this edge beyond this level
            val hasInterior: Boolean = false,       // Belongs to a shape of dimension 2.
            var a: R2Point = R2Point(),
            var b: R2Point = R2Point(),             // The edge endpoints, clipped to a given face
            val edge: S2Shape.Edge = S2Shape.Edge() // The edge endpoints
    )

    data class ClippedEdge(
        val faceEdge: FaceEdge,     // The original unclipped edge
        val bound: R2Rect           // Bounding box for the clipped portion
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
    private var pendingAdditionsBegin = 0;

    // The set of shapes that have been queued for removal but not processed
    // yet.  Note that we need to copy the edge data since the caller is free to
    // destroy the shape once Release() has been called.  This field is present
    // only when there are removed shapes to process (to save memory).
    private val pendingRemovals = mutableListOf<RemovedShape>()

    // Reads and writes to this field are guarded by "lock_".
    private val indexStatus = AtomicReference<IndexStatus>(IndexStatus.FRESH)

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
    //   for (MutableS2ShapeIndex::Iterator it(&index, S2ShapeIndex::BEGIN);
    //        !it.done(); it.Next()) { ... }
    inner class Iterator(pos: InitialPosition = InitialPosition.UNPOSITIONED) : S2ShapeIndexIteratorBase() {

        private lateinit var end: S2CellId
        private lateinit var keySet: NavigableSet<S2CellId>
        private var currentCellId: S2CellId = S2CellId.sentinel()

        init {
            init(pos)
        }

        // Initializes an iterator for the given MutableS2ShapeIndex.  This method
        // may also be called in order to restore an iterator to a valid state
        // after the underlying index has been updated (although it is usually
        // easier just to declare a new iterator whenever required, since iterator
        // construction is cheap).
        fun init(pos: InitialPosition) {
            maybeApplyUpdates()
            initStale(pos)
        }

        // Initialize an iterator for the given MutableS2ShapeIndex without
        // applying any pending updates.  This can be used to observe the actual
        // current state of the index without modifying it in any way.
        fun initStale(pos: InitialPosition) {
            end = cellMap.lastKey()
            cellMap.navigableKeySet()
            if (pos == InitialPosition.BEGIN) {
                currentCellId = keySet.first()
            } else {
                currentCellId = keySet.last();
            }
            refresh()
        }

        override fun begin() {
            // Make sure that the index has not been modified since Init() was called.
            Assertions.assert { isFresh() }
            currentCellId = keySet.first()
            refresh();
        }

        override fun finish() {
            currentCellId = keySet.last()
            refresh()
        }

        override fun next() {
            Assertions.assert { !done() }
            currentCellId = keySet.higher(currentCellId) ?: keySet.last()
            refresh();
        }

        override fun prev(): Boolean {
            if (currentCellId == keySet.first()) return false
            currentCellId = keySet.lower(currentCellId) ?: keySet.first()
            refresh()
            return true
        }

        override fun seek(target: S2CellId) {
            currentCellId = keySet.ceiling(target) ?: keySet.first()
            refresh();
        }

        private fun refresh() {
            if (currentCellId == end) {
                setFinished();
            } else {
                setState(currentCellId, cellMap.getValue(currentCellId))
            }
        }

    }

    override fun iterator(pos: InitialPosition): Iterator = Iterator(pos)

    // Takes ownership of the given shape and adds it to the index.  Also
    // assigns a unique id to the shape (shape->id()) and returns that id.
    // Shape ids are assigned sequentially starting from 0 in the order shapes
    // are added.  Invalidates all iterators and their associated data.
    fun add(shape: S2Shape): Int {
        // Additions are processed lazily by ApplyUpdates().
        shape.id = shapes.size
        shapes.add(shape);
        indexStatus.set(IndexStatus.STALE)
        return shape.id
    }

    // Removes the given shape from the index and return ownership to the caller.
    // Invalidates all iterators and their associated data.
    fun release(shapeId: Int): S2Shape = TODO()

    // Resets the index to its original state and returns ownership of all
    // shapes to the caller.  This method is much more efficient than removing
    // all shapes one at a time.
    fun releaseAll(): List<S2Shape> = TODO()

    // Resets the index to its original state and deletes all shapes.  Any
    // options specified via Init() are preserved.
    fun clear(): Unit = TODO()

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
    private fun is_first_update(): Boolean {
        // Note that it is not sufficient to check whether cell_map_ is empty, since
        // entries are added during the update process.
        return pendingAdditionsBegin == 0;
    }

    // Given that the given shape is being updated, return true if it is being
    // removed (as opposed to being added).
    private fun is_shape_being_removed(shape_id: Int): Boolean {
        // All shape ids being removed are less than all shape ids being added.
        return shape_id < pendingAdditionsBegin
    }

    // Apply any pending updates in a thread-safe way.
    private fun applyUpdatesThreadSafe(): Unit {
        lock.lock()
        if (indexStatus.get() == IndexStatus.FRESH) {
            lock.unlock()
        } else if (indexStatus.get() == IndexStatus.UPDATING) {
            // Wait until the updating thread is finished.  We do this by attempting
            // to lock a mutex that is held by the updating thread.  When this mutex
            // is unlocked the index_status_ is guaranteed to be FRESH.
            ++updateState.numWaiting
            lock.unlock()
            updateState.waitMutex.lock()
            lock.lock()
            --updateState.numWaiting
            unlockAndSignal();  // Notify other waiting threads.
        } else {
            Assertions.assertEQ(IndexStatus.STALE, indexStatus)
            indexStatus.set(IndexStatus.UPDATING)
            // Allocate the extra state needed for thread synchronization.  We keep
            // the spinlock held while doing this, because (1) memory allocation is
            // fast, so the chance of a context switch while holding the lock is low;
            // (2) by far the most common situation is that there is no contention,
            // and this saves an extra lock and unlock step; (3) even in the rare case
            // where there is contention, the main side effect is that some other
            // thread will burn a few CPU cycles rather than sleeping.
            updateState = UpdateState()
            // lock_.Lock wait_mutex *before* calling Unlock() to ensure that all other
            // threads will block on it.
            updateState.waitMutex.lock()
            // Release the spinlock before doing any real work.
            lock.unlock();
            applyUpdatesInternal();
            lock.lock();
            // index_status_ can be updated to FRESH only while locked *and* using
            // an atomic store operation, so that MaybeApplyUpdates() can check
            // whether the index is FRESH without acquiring the spinlock.
            indexStatus.set(IndexStatus.FRESH)
            unlockAndSignal();  // Notify any waiting threads.
        }
    }

    // This method updates the index by applying all pending additions and
    // removals.  It does *not* update index_status_ (see ApplyUpdatesThreadSafe).
    private fun applyUpdatesInternal() {
        // Check whether we have so many edges to process that we should process
        // them in multiple batches to save memory.  Building the index can use up
        // to 20x as much memory (per edge) as the final index size.
        val batches = getUpdateBatches()
        var i = 0
        for (batch in batches) {
            i++
            val all_edges: Array<MutableList<FaceEdge>> = Array(6) { mutableListOf<FaceEdge>() }
            logger.info { "Batch $i : shape_limit=${batch.additionsEnd}, edges=${batch.numEdges}" }

            reserveSpace(batch, all_edges)
            val tracker = InteriorTracker()
            if (pendingRemovals.isNotEmpty()) {
                // The first batch implicitly includes all shapes being removed.
                for (pendingRemoval in pendingRemovals) {
                    removeShape(pendingRemoval, all_edges, tracker)
                }
                pendingRemovals.clear()
            }
            for (id in pendingAdditionsBegin until batch.additionsEnd) {
                addShape(id, all_edges, tracker)
            }
            for (face in 0..5) {
                updateFaceEdges(face, all_edges[face], tracker)
                // Save memory by clearing vectors after we are done with them.
                //vector<FaceEdge>().swap(all_edges[face]);
            }
            pendingAdditionsBegin = batch.additionsEnd
        }
        // It is the caller's responsibility to update index_status_.
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
        var numEdges = numEdgesRemoved + numEdgesAdded;

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
        numEdges = 0;
        if (pendingRemovals.isNotEmpty()) {
            numEdges += numEdgesRemoved;
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
            numEdges = 0;
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
    fun getBatchSizes(numItems: Int, maxBatches: Int, finalBytesPerItem: Double, tmpBytesPerItem: Double, tmpMemoryBudgetBytes: Double) : List<Int> {
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
        val freeSpaceMultiplier = 1 - finalBytesRatio;

        // The total memory budget is the greater of the final size plus the allowed
        // temporary memory, or the minimum amount of memory required to limit the
        // number of batches to "max_batches".
        val totalBudgetBytes = max(finalBytes + tmpMemoryBudgetBytes, finalBytes / (1 - freeSpaceMultiplier.pow(maxBatches.toDouble())));

        // "max_batch_items" is the number of items in the current batch.
        var maxBatchItems = totalBudgetBytes / tmpBytesPerItem
        var i = 0
        while (i + 1 < maxBatches && remainingItems > 0) {
            val batchItems = min(remainingItems, (maxBatchItems + 1).toInt())
            batchSizes.add(batchItems)
            remainingItems -= batchItems;
            maxBatchItems *= freeSpaceMultiplier;
            ++i
        }
        Assertions.assertLE(batchSizes.size, maxBatches)
        return batchSizes
    }

    fun reserveSpace(batch: BatchDescriptor, allEdges: Array<MutableList<FaceEdge>>) {}

    // Clip all edges of the given shape to the six cube faces, add the clipped
    // edges to "all_edges", and start tracking its interior if necessary.
    fun addShape(id: Int, all_edges: Array<MutableList<FaceEdge>>, tracker:InteriorTracker) {
        val shape = shape(id) ?: return  // This shape has already been removed.
        // Construct a template for the edges to be added.
        val edge = FaceEdge(
                shapeId = id,
                hasInterior = shape.dimension == 2
        )
        if (edge.hasInterior) {
            tracker.addShape(id, S2ShapeUtil.containsBruteForce(shape, tracker.focus()))
        }
        val num_edges = shape.numEdges
        for (e in 0 until num_edges) {
            addFaceEdge(edge.copy(edgeId = e, edge = shape.edge(e), maxLevel = getEdgeMaxLevel(edge.edge)), all_edges)
        }
    }


    fun removeShape(removed: RemovedShape, allEdges: Array<MutableList<FaceEdge>>, tracker: InteriorTracker) {
        val edge = FaceEdge(
                shapeId = removed.shapeId,
                edgeId =  -1,  // Not used or needed for removed edges.
                hasInterior = removed.hasInterior,
        )
        if (edge.hasInterior) {
            tracker.addShape(edge.shapeId, removed.contains_tracker_origin);
        }
        for (removed_edge in removed.edges) {
            addFaceEdge(edge.copy(edge = removed_edge, maxLevel = getEdgeMaxLevel(edge.edge)), allEdges)
        }
    }

    fun addFaceEdge(edge: FaceEdge, allEdges: Array<MutableList<FaceEdge>>) {
        // Fast path: both endpoints are on the same face, and are far enough from
        // the edge of the face that don't intersect any (padded) adjacent face.
        val a_face = S2Coords.getFace(edge.edge.v0)
        if (a_face == S2Coords.getFace(edge.edge.v1)) {
            edge.a = S2Coords.validFaceXYZtoUV(a_face, edge.edge.v0)
            edge.b = S2Coords.validFaceXYZtoUV(a_face, edge.edge.v1)
            val kMaxUV = 1 - kCellPadding;
            if (abs(edge.a[0]) <= kMaxUV && abs(edge.a[1]) <= kMaxUV && abs(edge.b[0]) <= kMaxUV && abs(edge.b[1]) <= kMaxUV) {
            allEdges[a_face].add(edge);
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
    fun getEdgeMaxLevel(edge: S2Shape.Edge): Int {
        // Compute the maximum cell size for which this edge is considered "long".
        // The calculation does not need to be perfectly accurate, so we use Norm()
        // rather than Angle() for speed.
        val cell_size = ((edge.v0 - edge.v1).norm() * s2shape_index_cell_size_to_long_edge_ratio)
        // Now return the first level encountered during subdivision where the
        // average cell size is at most "cell_size".
        return S2CellMetrics.kAvgEdge.getLevelForMaxValue(cell_size);
    }

    fun updateFaceEdges(face: Int, face_edges: List<FaceEdge>, tracker: InteriorTracker): Unit = TODO()

/*
    protected:
    std::unique_ptr<IteratorBase> NewIterator(InitialPosition pos) const override;

    private:
    friend class EncodedS2ShapeIndex;
    friend class Iterator;
    friend class MutableS2ShapeIndexTest;
    friend class S2Stats;

    struct BatchDescriptor;
    struct ClippedEdge;
    class EdgeAllocator;
    struct FaceEdge;
    class InteriorTracker;
    struct RemovedShape;


    // Internal methods are documented with their definitions.
    void RemoveShape(const RemovedShape& removed,std::vector<FaceEdge> all_edges[6],InteriorTracker* tracker) const;
    void AddFaceEdge(FaceEdge* edge, std::vector<FaceEdge> all_edges[6]) const;
    void UpdateFaceEdges(int face, const std::vector<FaceEdge>& face_edges,InteriorTracker* tracker);
    S2CellId ShrinkToFit(const S2PaddedCell& pcell, const R2Rect& bound) const;
    void SkipCellRange(S2CellId begin, S2CellId end, InteriorTracker* tracker,EdgeAllocator* alloc, bool disjoint_from_index);
    void UpdateEdges(const S2PaddedCell& pcell,std::vector<const ClippedEdge*>* edges,InteriorTracker* tracker, EdgeAllocator* alloc,bool disjoint_from_index);
    void AbsorbIndexCell(const S2PaddedCell& pcell,const Iterator& iter,std::vector<const ClippedEdge*>* edges,InteriorTracker* tracker,EdgeAllocator* alloc);
    int GetEdgeMaxLevel(const S2Shape::Edge& edge) const;
    static int CountShapes(const std::vector<const ClippedEdge*>& edges,const ShapeIdSet& cshape_ids);
    bool MakeIndexCell(const S2PaddedCell& pcell,const std::vector<const ClippedEdge*>& edges,InteriorTracker* tracker);
    static void TestAllEdges(const std::vector<const ClippedEdge*>& edges,InteriorTracker* tracker);
    inline static const ClippedEdge* UpdateBound(const ClippedEdge* edge,int u_end, double u,int v_end, double v,EdgeAllocator* alloc);
    static const ClippedEdge* ClipUBound(const ClippedEdge* edge,int u_end, double u,EdgeAllocator* alloc);
    static const ClippedEdge* ClipVBound(const ClippedEdge* edge,int v_end, double v,EdgeAllocator* alloc);
    static void ClipVAxis(const ClippedEdge* edge, const R1Interval& middle, std::vector<const ClippedEdge*> child_edges[2], EdgeAllocator* alloc);
*/
    // Documented in the .cc file.
    fun unlockAndSignal(): Unit = TODO()
/*
    UNLOCK_FUNCTION(lock_)
    UNLOCK_FUNCTION(update_state_->wait_mutex);
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
        val kCellPadding: Double = 2 * (dilivia.s2.S2EdgeClipping.kFaceClipErrorUVCoord + dilivia.s2.S2EdgeClipping.kEdgeClipErrorUVCoord)

        // Defines the initial focus point of MutableS2ShapeIndex::InteriorTracker
        // (the start of the S2CellId space-filling curve).
        //
        // TODO(ericv): Move InteriorTracker here to avoid the need for this method.
        fun kInteriorTrackerOrigin() = S2Coords.faceUVtoXYZ(0, -1.0, -1.0).normalize();

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
        val kFinalBytesPerEdge = 8
        val kTmpBytesPerEdge = 200
        val kTmpMemoryBudgetBytes = 100 shl 20// static_cast<size_t>(FLAGS_s2shape_index_tmp_memory_budget_mb) << 20;

        // We arbitrarily limit the number of batches just as a safety measure.
        // With the current default memory budget of 100 MB, this limit is not
        // reached even when building an index of 350 million edges.
        val kMaxUpdateBatches = 100;

    }
}