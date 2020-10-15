package dilivia.s2.index

import com.google.common.collect.ComparisonChain
import dilivia.s2.S2CellId
import dilivia.s2.region.S2CellUnion
import dilivia.s2.collections.sortAndRemoveDuplicates
import dilivia.s2.collections.upperBound
import mu.KotlinLogging

/**
 * Labels are 32-bit non-negative integers.  To support other label types, you can use ValueLexicon to map label values
 * to integers:
 *
 * val my_label_lexicon = ValueLexicon<MyLabel>()
 * index.add(cell_id, my_label_lexicon.add(label))
 */
typealias Label = Int

/**
// S2CellIndex stores a collection of (cell_id, label) pairs.  The S2CellIds
// may be overlapping or contain duplicate values.  For example, an
// S2CellIndex could store a collection of S2CellUnions, where each
// S2CellUnion has its own label.
//
// Labels are 32-bit non-negative integers, and are typically used to map the
// results of queries back to client data structures.  Labels other than
// integers can be supported by using a ValueLexicon, which maintains a set of
// distinct labels and maps them to sequentially numbered integers.  For
// example, the following code uses strings as labels:
//
//   ValueLexicon<string> my_label_lexicon;
//   string label_str = ...;
//   cell_index.Add(cell_id, my_label_lexicon.Add(label_str));
//   ...
//   int32 label = ...;
//   string label_str = my_label_lexicon.value(label);
//
// To build an S2CellIndex, call Add() for each (cell_id, label) pair, and
// then call the Build() method.  For example:
//
//   vector<S2CellId> contents = ...;
//   for (int i = 0; i < contents.size(); ++i) {
//     index.Add(contents[i], i /*label*/);
//   }
//   index.Build();
//
// There is also a convenience method that adds an S2CellUnion:
//
//     index.Add(cell_union, label);
//
// Note that the index is not dynamic; the contents of the index cannot be
// changed once it has been built.
//
// There are several options for retrieving data from the index.  The simplest
// is to use a built-in method such as GetIntersectingLabels (which returns
// the labels of all cells that intersect a given target S2CellUnion):
//
//   vector<Label> labels = index.GetIntersectingLabels(target_union);
//
// Alternatively, you can use an external class such as S2ClosestCellQuery,
// which computes the cell(s) that are closest to a given target geometry.
// For example, here is how to find all cells that are closer than
// "distance_limit" to a given target point:
//
//   S2ClosestCellQuery query(&index);
//   query.mutable_options()->set_max_distance(distance_limit);
//   S2ClosestCellQuery::PointTarget target(target_point);
//   for (const auto& result : query.FindClosestCells(&target)) {
//     // result.distance() is the distance to the target.
//     // result.cell_id() is the indexed S2CellId.
//     // result.label() is the integer label associated with the S2CellId.
//     DoSomething(target_point, result);
//   }
//
// Finally, you can access the index contents directly.  Internally, the index
// consists of a set of non-overlapping leaf cell ranges that subdivide the
// sphere and such that each range intersects a particular set of (cell_id,
// label) pairs.  Data is accessed using the following iterator types:
//
//   RangeIterator:
//    - used to seek and iterate over the non-overlapping leaf cell ranges.
//   NonEmptyRangeIterator:
//    - like RangeIterator, but skips ranges whose contents are empty.
//   ContentsIterator:
//    - iterates over the (cell_id, label) pairs that intersect a given range.
//   CellIterator:
//    - iterates over the entire set of (cell_id, label) pairs.
//
// Note that these are low-level, efficient types intended mainly for
// implementing new query classes.  Most clients should use either the
// built-in methods such as VisitIntersectingCells and GetIntersectingLabels,
// or a helper such as S2ClosestCellQuery or S2Closest*Query::CellUnionTarget.
 */
class S2CellIndex() {

    // A tree of (cell_id, label) pairs such that if X is an ancestor of Y, then
    // X.cell_id contains Y.cell_id.  The contents of a given range of leaf
    // cells can be represented by pointing to a node of this tree.
    private val cellTree = mutableListOf<CellNode>()

    // The last element of range_nodes_ is a sentinel value, which is necessary
    // in order to represent the range covered by the previous element.
    private val rangeNodes = ArrayList<RangeNode>()

    // Returns the number of (cell_id, label) pairs in the index.
    fun numCells(): Int = cellTree.size

    // Adds the given (cell_id, label) pair to the index.  Note that the index
    // is not valid until Build() is called.
    //
    // The S2CellIds in the index may overlap (including duplicate values).
    // Duplicate (cell_id, label) pairs are also allowed, although be aware that
    // S2ClosestCellQuery will eliminate such duplicates anyway.
    //
    // REQUIRES: cell_id.is_valid()
    fun add(cellId: S2CellId, label: Label) {
        check(cellId.isValid())
        check(label >= 0)
        cellTree.add(CellNode(cellId, label, -1))
    }

    // Convenience function that adds a collection of cells with the same label.
    fun add(cellIds: S2CellUnion, label: Label) {
        cellIds.forEach { cellId -> add(cellId, label) }
    }

    // Constructs the index.  This method may only be called once.  No iterators
    // may be used until the index is built.
    fun build() {
        val startTime = System.currentTimeMillis()
        logger.debug { "Start building index" }

        val deltas = ArrayList<Delta>(2 * cellTree.size + 2)
        // Create two deltas for each (cell_id, label) pair: one to add the pair to
        // the stack (at the start of its leaf cell range), and one to remove it from
        // the stack (at the end of its leaf cell range).
        for (node in cellTree) {
            deltas.add(Delta(node.cellId.rangeMin(), node.cellId, node.label))
            deltas.add(Delta(node.cellId.rangeMax().next(), S2CellId.sentinel(), -1))
        }
        // We also create two special deltas to ensure that a RangeNode is emitted at
        // the beginning and end of the S2CellId range.
        deltas.add(Delta(S2CellId.begin(S2CellId.kMaxLevel), S2CellId.none(), -1))
        deltas.add(Delta(S2CellId.end(S2CellId.kMaxLevel), S2CellId.none(), -1))
        deltas.sort()
        logger.trace { "Deltas: ${deltas.size} elements\n-----------------------------\n${deltas.joinToString("\n")}" }

        // Now walk through the deltas to build the leaf cell ranges and cell tree
        // (which is essentially a permanent form of the "stack" described above).
        cellTree.clear()
        rangeNodes.ensureCapacity(deltas.size)
        var contents = -1
        var i = 0
        while (i < deltas.size) {
            val start_id = deltas[i].startId
            // Process all the deltas associated with the current start_id.
            while (i < deltas.size && deltas[i].startId == start_id) {
                if (deltas[i].label >= 0) {
                    cellTree.add(CellNode(cellId = deltas[i].cellId, label = deltas[i].label, parent = contents))
                    contents = cellTree.size - 1
                } else if (deltas[i].cellId == S2CellId.sentinel()) {
                    contents = cellTree[contents].parent
                }
                ++i
            }
            rangeNodes.add(RangeNode(startId = start_id, contents = contents))
        }

        logger.trace { "Cell tree: ${cellTree.size} elements\n-----------------------------\n${cellTree.joinToString("\n")}" }
        logger.trace { "Range nodes: ${rangeNodes.size} elements\n-----------------------------\n${rangeNodes.joinToString("\n")}" }
        logger.debug { "Index build in ${System.currentTimeMillis() - startTime} ms" }
    }

    // Clears the index so that it can be re-used.
    fun clear() {
        cellTree.clear()
        rangeNodes.clear()
    }

    // Visits all (cell_id, label) pairs in the given index that intersect the
    // given S2CellUnion "target", terminating early if the given CellVisitor
    // function returns false (in which case VisitIntersectingCells returns false
    // as well).  Each (cell_id, label) pair in the index is visited at most
    // once.  (If the index contains duplicates, then each copy is visited.)
    fun visitIntersectingCells(target: S2CellUnion, visitor: CellVisitor): Boolean {
        logger.trace { "--> visitIntersectingCells(target = $target)" }
        if (target.isEmpty()) return true
        val targetIterator = target.listIterator()
        var it = targetIterator.next()
        val contents = ContentsIterator(this)
        val range = RangeIterator(this)
        range.begin()
        do {
            logger.trace { """
                |
                |Current target cell: $it [ range min = ${it.rangeMin()} ; range max = ${it.rangeMax()} ]
                |=====================================
                |Current range: 
                |  - start id:  ${range.startId()}
                |  - limit id : ${range.limitId()}
                |---------------------------------------
            """.trimMargin() }


            if (range.limitId() <= it.rangeMin()) {
                logger.trace { "Range limit id <= current cell range min id => seek to cell range min = ${it.rangeMin()}" }
                range.seek(it.rangeMin())  // Only seek when necessary.
            }
            while (range.startId() <= it.rangeMax()) {
                logger.trace { "Range start id = ${range.startId()} <= target cell range max = ${it.rangeMax()} => visit range contents" }
                contents.startUnion(range)
                while (!contents.done()) {
                    logger.trace { "Visit cell ${LabelledCell(cellId = contents.cellId(), label = contents.label())}" }
                    if (!visitor.visit(contents.cellId(), contents.label())) {
                        return false;
                    }
                    contents.next()
                }
                range.next()
            }
            // Check whether the next target cell is also contained by the leaf cell
            // range that we just processed.  If so, we can skip over all such cells
            // using binary search.  This speeds up benchmarks by between 2x and 10x
            // when the average number of intersecting cells is small (< 1).
            if (!range.done() && targetIterator.hasNext()) {
                it = targetIterator.next()
                if (it.rangeMax() < range.startId()) {
                    logger.trace { "Next cell ($it) range max = ${it.rangeMax()} < range start id = ${range.startId()}" }
                    // Skip to the first target cell that extends past the previous range.
                    // Skip to the first target cell that extends past the previous range.
                    //it = std::lower_bound(it + 1, target.end(), range.start_id());
                    //if ((it - 1)->range_max() >= range.start_id())--it;
                    while (targetIterator.hasNext() && target[targetIterator.nextIndex()] < range.startId()) {
                        it = targetIterator.next()
                    }
                    logger.trace { "Current cell = $it (idx = ${targetIterator.previousIndex() + 1}) after Skip to the first target cell that extends past the previous range" }
                    if (target[targetIterator.previousIndex()].rangeMax() >= range.startId()) {
                        logger.trace { "Previous cell = ${target[targetIterator.previousIndex()]} range max >= range start id = ${range.startId()} => previous" }
                        it = targetIterator.previous()
                    }
                }
            } else break
        } while (true)
        return true
    }

    // Convenience function that returns the labels of all indexed cells that
    // intersect the given S2CellUnion "target".  The output contains each label
    // at most once, but is not sorted.
    fun getIntersectingLabels(target: S2CellUnion): List<Label> {
        val labels = mutableListOf<Label>()
        getIntersectingLabels(target, labels)
        return labels
    }

    // This version can be more efficient when it is called many times, since it
    // does not require allocating a new vector on each call.
    fun getIntersectingLabels(target: S2CellUnion, labels: MutableList<Label>) {
        labels.clear()
        visitIntersectingCells(target, object : CellVisitor {

            override fun visit(cellId: S2CellId, label: Label): Boolean {
                labels.add(label)
                return true
            }
        })
        labels.sortAndRemoveDuplicates()
    }

    // A function that is called with each (cell_id, label) pair to be visited.
    // The function may return false in order to indicate that no further
    // (cell_id, label) pairs are needed.
    interface CellVisitor {
        fun visit(cellId: S2CellId, label: Label): Boolean
    }

    // Convenience class that represents a (cell_id, label) pair.
    data class LabelledCell(
            val cellId: S2CellId = S2CellId.none(),
            val label: Label = -1
    ) : Comparable<LabelledCell> {

        override fun compareTo(other: LabelledCell): Int = ComparisonChain.start()
                .compare(cellId, other.cellId)
                .compare(label, other.label)
                .result()

    }

    // Represents a node in the (cell_id, label) tree.  Cells are organized in a
    // tree such that the ancestors of a given node contain that node.
    data class CellNode(
            val cellId: S2CellId = S2CellId.none(),
            val label: Label = kDoneContents,
            val parent: Int = -1
    )

    // An iterator that visits the entire set of indexed (cell_id, label) pairs
    // in an unspecified order.
    // Initializes a CellIterator for the given S2CellIndex, positioned at the
    // first cell (if any).
    class CellIterator(val index: S2CellIndex) {

        // NOTE(ericv): There is a potential optimization that would require this
        // class to iterate over both cell_tree_ *and* range_nodes_.
        private var cellIdx: Int = -1

        init {
            cellIdx = 0
            check(!index.rangeNodes.isEmpty()) { "Call build() first." }
        }

        // The S2CellId of the current (cell_id, label) pair.
        // REQUIRES: !done()
        fun cellId(): S2CellId {
            check(!done()) { "End of cell iterator is reached. " }
            return index.cellTree[cellIdx].cellId
        }

        // The Label of the current (cell_id, label) pair.
        // REQUIRES: !done()
        fun label(): Label {
            check(!done()) { "End of cell iterator is reached. " }
            return index.cellTree[cellIdx].label
        }

        // Returns the current (cell_id, label) pair.
        fun labelledCell(): LabelledCell {
            check(!done()) { "End of cell iterator is reached. " }
            val cellNode = index.cellTree[cellIdx]
            return LabelledCell(cellId = cellNode.cellId, label = cellNode.label)
        }

        // Returns true if all (cell_id, label) pairs have been visited.
        fun done(): Boolean = cellIdx == index.cellTree.size

        // Advances the iterator to the next (cell_id, label) pair.
        // REQUIRES: !done()
        fun next(): Unit {
            check(!done()) { "End of cell iterator is reached. " }
            ++cellIdx
        }

    }

    // A RangeNode represents a range of leaf S2CellIds.  The range starts at
    // "start_id" (a leaf cell) and ends at the "start_id" field of the next
    // RangeNode.  "contents" points to the node of cell_tree_ representing the
    // cells that overlap this range.
    private data class RangeNode(
            val startId: S2CellId,  // First leaf cell contained by this range.
            val contents: Int     // Contents of this node (an index within cell_tree_).
    ) : Comparable<Any?> {

        override fun compareTo(other: Any?): Int {
            return when (other) {
                is S2CellId -> compareToCell(other)
                is RangeNode -> compareToRange(other)
                else -> throw IllegalArgumentException("")
            }

        }

        fun compareToCell(other: S2CellId): Int = startId.compareTo(other)

        fun compareToRange(other: RangeNode): Int = ComparisonChain.start()
                .compare(startId, other.startId)
                .compare(contents, other.contents)
                .result()
    }

    // An iterator that seeks and iterates over a set of non-overlapping leaf
    // cell ranges that cover the entire sphere.  The indexed (s2cell_id, label)
    // pairs that intersect the current leaf cell range can be visited using
    // ContentsIterator (see below).
    // Initializes a RangeIterator for the given S2CellIndex.  The iterator is
    // initially *unpositioned*; you must call a positioning method such as
    // Begin() or Seek() before accessing its contents.
    open class RangeIterator(val index: S2CellIndex) {

        private var rangeNodeIdx: Int

        init {
            check(index.rangeNodes.isNotEmpty()) { "Call build() first." }
            rangeNodeIdx = kUninitialized
        }

        constructor(iterator: RangeIterator): this(iterator.index) {
            this.rangeNodeIdx = iterator.rangeNodeIdx
        }

        // The start of the current range of leaf S2CellIds.
        //
        // If done() is true, returns S2CellId::End(S2CellId::kMaxLevel).  This
        // property means that most loops do not need to test done() explicitly.
        fun startId(): S2CellId = index.rangeNodes[rangeNodeIdx].startId

        // The (non-inclusive) end of the current range of leaf S2CellIds.
        // REQUIRES: !done()
        fun limitId(): S2CellId {
            check(!done())
            return index.rangeNodes[rangeNodeIdx + 1].startId
        }

        // Returns true if the iterator is positioned beyond the last valid range.
        fun done(): Boolean {
            check(rangeNodeIdx != kUninitialized) { "Call begin() or seek() first." }

            // Note that the last element of range_nodes_ is a sentinel value.
            val done = rangeNodeIdx >= index.rangeNodes.lastIndex
            return done
        }

        // Positions the iterator at the first range of leaf cells (if any).
        open fun begin() {
            rangeNodeIdx = 0
        }

        // Positions the iterator so that done() is true.
        fun finish(): Unit {
            rangeNodeIdx = index.rangeNodes.lastIndex
        }

        // Advances the iterator to the next range of leaf cells.
        // REQUIRES: !done()
        open fun next() {
            check(!done())
            ++rangeNodeIdx
        }

        // If the iterator is already positioned at the beginning, returns false.
        // Otherwise positions the iterator at the previous entry and returns true.
        open fun prev(): Boolean {
            if (rangeNodeIdx == 0) return false
            --rangeNodeIdx
            return true
        }

        // Positions the iterator at the first range with start_id() >= target.
        // (Such an entry always exists as long as "target" is a valid leaf cell.
        // Note that it is valid to access start_id() even when done() is true.)
        //
        // REQUIRES: target.is_leaf()
        open fun seek(target: S2CellId) {
            check(target.isLeaf())
            rangeNodeIdx = index.rangeNodes.upperBound(value = target) - 1// std::upper_bound(range_nodes_->begin(), range_nodes_->end(), target)-1;
        }

        // Returns true if no (s2cell_id, label) pairs intersect this range.
        // Also returns true if done() is true.
        fun isEmpty(): Boolean = index.rangeNodes[rangeNodeIdx].contents == kDoneContents

        // If advancing the iterator "n" times would leave it positioned on a
        // valid range, does so and returns true.  Otherwise leaves the iterator
        // unmodified and returns false.
        fun advance(n: Int): Boolean {
            // Note that the last element of range_nodes_ is a sentinel value.
            if (n >= index.rangeNodes.lastIndex - 1 - rangeNodeIdx) return false
            rangeNodeIdx += n
            return true
        }

        fun contents(): Int = index.rangeNodes[rangeNodeIdx].contents

        companion object {
            // A special value used to indicate that the RangeIterator has not yet
            // been initialized by calling Begin() or Seek().
            // Note that since the last element of range_nodes_ is a sentinel value,
            // it_ will never legitimately be positioned at range_nodes_->end().
            const val kUninitialized = Int.MAX_VALUE
        }

    }

    // Like RangeIterator, but only visits leaf cell ranges that overlap at
    // least one (cell_id, label) pair.
    // Initializes a NonEmptyRangeIterator for the given S2CellIndex.
    // The iterator is initially *unpositioned*; you must call a positioning
    // method such as Begin() or Seek() before accessing its contents.
    class NonEmptyRangeIterator(index: S2CellIndex) : RangeIterator(index) {

        // Positions the iterator at the first non-empty range of leaf cells.
        override fun begin() {
            super.begin()
            while (isEmpty() && !done()) next()
        }

        // Advances the iterator to the next non-empty range of leaf cells.
        // REQUIRES: !done()
        override fun next() {
            do {
                super.next()
            } while (isEmpty() && !done())
        }

        // If the iterator is already positioned at the beginning, returns false.
        // Otherwise positions the iterator at the previous entry and returns true.
        override fun prev(): Boolean {
            while (super.prev()) {
                if (!isEmpty()) return true
            }
            // Return the iterator to its original position.
            if (isEmpty() && !done()) next()
            return false
        }

        // Positions the iterator at the first non-empty range with
        // start_id() >= target.
        //
        // REQUIRES: target.is_leaf()
        override fun seek(target: S2CellId) {
            super.seek(target)
            while (isEmpty() && !done()) super.next()
        }

    }

    // An iterator that visits the (cell_id, label) pairs that cover a set of
    // leaf cell ranges (see RangeIterator).  Note that when multiple leaf cell
    // ranges are visited, this class only guarantees that each result will be
    // reported at least once, i.e. duplicate values may be suppressed.  If you
    // want duplicate values to be reported again, be sure to call Clear() first.
    //
    // [In particular, the implementation guarantees that when multiple leaf
    // cell ranges are visited in monotonically increasing order, then each
    // (cell_id, label) pair is reported exactly once.]
    class ContentsIterator {

        private var index: S2CellIndex? = null

        // The value of it.start_id() from the previous call to StartUnion().
        // This is used to check whether these values are monotonically
        // increasing.
        private var prevStartId: S2CellId = S2CellId.none()

        // The maximum index within the cell_tree_ vector visited during the
        // previous call to StartUnion().  This is used to eliminate duplicate
        // values when StartUnion() is called multiple times.
        private var nodeCutoff: Int = -1

        // The maximum index within the cell_tree_ vector visited during the
        // current call to StartUnion().  This is used to update node_cutoff_.
        private var nextNodeCutoff: Int = -1

        // A copy of the current node in the cell tree.
        private var node: CellNode = CellNode()

        // Default constructor; must be followed by a call to Init().
        constructor() {

        }

        // Convenience constructor that calls Init().
        constructor(index: S2CellIndex) {
            init(index)
        }

        // Initializes the iterator.  Should be followed by a call to UnionWith()
        // to visit the contents of each desired leaf cell range.
        fun init(index: S2CellIndex): Unit {
            this.index = index
            clear()
        }

        // Clears all state with respect to which range(s) have been visited.
        fun clear() {
            prevStartId = S2CellId.none()
            nodeCutoff = -1
            nextNodeCutoff = -1
            setDone()
        }

        // Positions the ContentsIterator at the first (cell_id, label) pair that
        // covers the given leaf cell range.  Note that when multiple leaf cell
        // ranges are visited using the same ContentsIterator, duplicate values
        // may be suppressed.  If you don't want this behavior, call Clear() first.
        fun startUnion(range: RangeIterator) {
            check(index != null)
            if (range.startId() < prevStartId) {
                nodeCutoff = -1;  // Can't automatically eliminate duplicates.
            }
            prevStartId = range.startId()

            // TODO(ericv): Since RangeNode only uses 12 of its 16 bytes, we could add a
            // "label" field without using any extra space.  Then we could store a leaf
            // node of cell_tree_ directly in each RangeNode, where the cell_id is
            // implicitly defined as the one that covers the current leaf cell range.
            // This would save quite a bit of space; e.g. if the given cells are
            // non-overlapping, then cell_tree_ would be empty (since every node is a
            // leaf node and could therefore be stored directly in a RangeNode).  It
            // would also be faster because cell_tree_ would rarely be accessed.
            val contents = range.contents()
            if (contents <= nodeCutoff) {
                setDone()
            } else {
                node = index!!.cellTree[contents].copy()
            }

            // When visiting ancestors, we can stop as soon as the node index is smaller
            // than any previously visited node index.  Because indexes are assigned
            // using a preorder traversal, such nodes are guaranteed to have already
            // been reported.
            nextNodeCutoff = contents

        }

        // The S2CellId of the current (cell_id, label) pair.
        // REQUIRES: !done()
        fun cellId(): S2CellId {
            check(!done())
            return node.cellId
        }

        // The Label of the current (cell_id, label) pair.
        // REQUIRES: !done()
        fun label(): Label {
            check(!done())
            return node.label
        }

        // Returns the current (cell_id, label) pair.
        // REQUIRES: !done()
        fun labelledCell(): LabelledCell {
            check(!done())
            return LabelledCell(node.cellId, node.label)
        }

        // Returns true if all (cell_id, label) pairs have been visited.
        fun done(): Boolean = node.label == kDoneContents

        // Advances the iterator to the next (cell_id, label) pair covered by the
        // current leaf cell range.
        // REQUIRES: !done()
        fun next() {
            check(index != null)
            check(!done())
            if (node.parent <= nodeCutoff) {
                // We have already processed this node and its ancestors.
                nodeCutoff = nextNodeCutoff
                setDone()
            } else {
                node = index!!.cellTree[node.parent].copy()
            }

        }

        // node_.label == kDoneContents indicates that done() is true.
        private fun setDone() {
            node = node.copy(label = kDoneContents)
        }

    }

    // To build the cell tree and leaf cell ranges, we maintain a stack of
    // (cell_id, label) pairs that contain the current leaf cell.  This class
    // represents an instruction to push or pop a (cell_id, label) pair.
    //
    // If label >= 0, the (cell_id, label) pair is pushed on the stack.
    // If cell_id == S2CellId::Sentinel(), a pair is popped from the stack.
    // Otherwise the stack is unchanged but a RangeNode is still emitted.
    private data class Delta(val startId: S2CellId, val cellId: S2CellId, val label: Label) : Comparable<Delta> {

        // Deltas are sorted first by start_id, then in reverse order by cell_id,
        // and then by label.  This is necessary to ensure that (1) larger cells
        // are pushed on the stack before smaller cells, and (2) cells are popped
        // off the stack before any new cells are added.

        override fun compareTo(other: Delta): Int {
            return ComparisonChain.start()
                    .compare(startId, other.startId)
                    .compare(other.cellId, cellId)
                    .compare(label, other.label)
                    .result()
        }

    }

    companion object {

        private val logger = KotlinLogging.logger(S2CellIndex::class.java.name)

        // A special label indicating that ContentsIterator::done() is true.
        const val kDoneContents: Label = -1;

    }


}
