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

import dilivia.s2.region.S2CellUnion

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
class S2CellIndex<T> where T : Comparable<T> {

    // Convenience class that represents a (cell_id, label) pair.
    data class LabelledCell<T>(val cellId: S2CellId = S2CellId.none(), val label: T) : Comparable<LabelledCell<T>> where T : Comparable<T> {

        /**
         * Compares this object with the specified object for order. Returns zero if this object is equal
         * to the specified [other] object, a negative number if it's less than [other], or a positive number
         * if it's greater than [other].
         */
        override fun compareTo(other: LabelledCell<T>): Int {
            val cellIdComparison = cellId.compareTo(other.cellId)
            return if (cellIdComparison != 0) cellIdComparison else label.compareTo(other.label)
        }

    }

    // Represents a node in the (cell_id, label) tree.  Cells are organized in a
    // tree such that the ancestors of a given node contain that node.
    data class CellNode<T>(val cellId: S2CellId = S2CellId.none(), val label: T, val parent: Int = -1)

    // A RangeNode represents a range of leaf S2CellIds.  The range starts at
    // "start_id" (a leaf cell) and ends at the "start_id" field of the next
    // RangeNode.  "contents" points to the node of cell_tree_ representing the
    // cells that overlap this range.
    // @property startId First leaf cell contained by this range.
    // @property contents Contents of this node (an index within cell_tree).
    data class RangeNode(val startId: S2CellId, val contents: Int) {

        // TODO(fmeurisse)
        // Comparison operator needed for std::upper_bound().
//        friend bool operator <(S2CellId x, const RangeNode& y)
//        {
//            return x < y.start_id;
//        }
    }

    // A function that is called with each (cell_id, label) pair to be visited.
    // The function may return false in order to indicate that no further
    // (cell_id, label) pairs are needed.
    @FunctionalInterface
    interface CellVisitor<T> {

        fun apply(cellId: S2CellId, label: T): Boolean

    }


    // An iterator that visits the entire set of indexed (cell_id, label) pairs
    // in an unspecified order.
    // NOTE(ericv): There is a potential optimization that would require this
    // class to iterate over both cell_tree_ *and* range_nodes_.
    inner class CellIterator(private val iterator: Iterator<CellNode<T>>) {

        private lateinit var currentCellNode: CellNode<T>

        // The S2CellId of the current (cell_id, label) pair.
        fun currentCellId(): S2CellId {
            check(this::currentCellNode.isInitialized)
            return currentCellNode.cellId
        }

        // The Label of the current (cell_id, label) pair.
        fun currentLabel(): T {
            check(this::currentCellNode.isInitialized)
            return currentCellNode.label
        }

        // Returns the current (cell_id, label) pair.
        fun currentLabelledCell(): LabelledCell<T> {
            check(this::currentCellNode.isInitialized)
            return LabelledCell(cellId = currentCellNode.cellId, label = currentCellNode.label)
        }

        // Returns true if all (cell_id, label) pairs have been visited.
        fun hasNext(): Boolean = iterator.hasNext()

        // Advances the iterator to the next (cell_id, label) pair.
        fun next(): LabelledCell<T> {
            currentCellNode = iterator.next()
            return currentLabelledCell()
        }

    }

    // An iterator that seeks and iterates over a set of non-overlapping leaf
    // cell ranges that cover the entire sphere.  The indexed (s2cell_id, label)
    // pairs that intersect the current leaf cell range can be visited using
    // ContentsIterator (see below).
    open inner class RangeIterator() {

        protected var currentRangeNodeIdx: Int = -1
        protected lateinit var currentRangeNode: RangeNode

        // The start of the current range of leaf S2CellIds.
        //
        // If done() is true, returns S2CellId::End(S2CellId::kMaxLevel).  This
        // property means that most loops do not need to test done() explicitly.
        fun startId(): S2CellId {
            check(this::currentRangeNode.isInitialized)
            return currentRangeNode.startId
        }

        // The (non-inclusive) end of the current range of leaf S2CellIds.
        // REQUIRES: !done()
        fun limitId(): S2CellId {
            check(this.currentRangeNodeIdx != -1)
            check(this.currentRangeNodeIdx <= range_nodes.lastIndex - 1)
            return range_nodes[currentRangeNodeIdx + 1].startId
        }

        // Returns true if the iterator is positioned beyond the last valid range.
        fun done(): Boolean {
            check(currentRangeNodeIdx != -1) { "Call begin() or seek() first." }
            return currentRangeNodeIdx >= range_nodes.lastIndex
        }

        // Positions the iterator at the first range of leaf cells (if any).
        open fun begin(): Unit {
            currentRangeNodeIdx = 0
            currentRangeNode = range_nodes[currentRangeNodeIdx]
        }

        // Positions the iterator so that done() is true.
        fun finish(): Unit {
            currentRangeNodeIdx = range_nodes.lastIndex
            currentRangeNode = range_nodes[currentRangeNodeIdx]
        }

        // Advances the iterator to the next range of leaf cells.
        // REQUIRES: !done()
        open fun next(): Unit {
            check(!done())
            currentRangeNodeIdx++
            currentRangeNode = range_nodes[currentRangeNodeIdx]
        }

        // If the iterator is already positioned at the beginning, returns false.
        // Otherwise positions the iterator at the previous entry and returns true.
        open fun previous(): Boolean {
            if (currentRangeNodeIdx != 0) {
                currentRangeNodeIdx--
                currentRangeNode = range_nodes[currentRangeNodeIdx]
                return true
            }
            return false
        }

        // Positions the iterator at the first range with start_id() >= target.
        // (Such an entry always exists as long as "target" is a valid leaf cell.
        // Note that it is valid to access start_id() even when done() is true.)
        //
        // REQUIRES: target.is_leaf()
        open fun seek(target: S2CellId): Unit = TODO()

        // Returns true if no (s2cell_id, label) pairs intersect this range.
        // Also returns true if done() is true.
        fun isEmpty(): Boolean {
            return currentRangeNode.contents == -1
        }

        // If advancing the iterator "n" times would leave it positioned on a
        // valid range, does so and returns true.  Otherwise leaves the iterator
        // unmodified and returns false.
        fun advance(n: Int): Boolean {
            check(currentRangeNodeIdx != -1) { "Call begin() or seek() first." }
            // Note that the last element of range_nodes_ is a sentinel value.
            if (n >= range_nodes.lastIndex - currentRangeNodeIdx) return false;
            currentRangeNodeIdx += n
            currentRangeNode = range_nodes[currentRangeNodeIdx]
            return true
        }

    }

    // Like RangeIterator, but only visits leaf cell ranges that overlap at
    // least one (cell_id, label) pair.
    inner class NonEmptyRangeIterator() : RangeIterator() {

        // Positions the iterator at the first non-empty range of leaf cells.
        override fun begin() {
            super.begin()
            while (isEmpty() && !done()) super.next()
        }

        // Advances the iterator to the next non-empty range of leaf cells.
        // REQUIRES: !done()
        override fun next() {
            do {
                super.next();
            } while (isEmpty() && !done())
        }

        // If the iterator is already positioned at the beginning, returns false.
        // Otherwise positions the iterator at the previous entry and returns true.
        override fun previous(): Boolean {
            while (super.previous()) {
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
        override fun seek(target: S2CellId): Unit {
            super.seek(target);
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
    inner class ContentsIterator {

        // The value of it.start_id() from the previous call to StartUnion().
        // This is used to check whether these values are monotonically
        // increasing.
        private var prev_start_id: S2CellId = S2CellId.none()

        // The maximum index within the cell_tree_ vector visited during the
        // previous call to StartUnion().  This is used to eliminate duplicate
        // values when StartUnion() is called multiple times.
        private var node_cutoff: Int = -1

        // The maximum index within the cell_tree_ vector visited during the
        // current call to StartUnion().  This is used to update node_cutoff_.
        private var next_node_cutoff: Int = -1

        // A copy of the current node in the cell tree.
        protected var currentNodeIdx: Int = -1
        private lateinit var currentNode: CellNode<T>

        // Clears all state with respect to which range(s) have been visited.
        fun clear() {
            prev_start_id = S2CellId.none()
            node_cutoff = -1;
            next_node_cutoff = -1;
        }

        // Positions the ContentsIterator at the first (cell_id, label) pair that
        // covers the given leaf cell range.  Note that when multiple leaf cell
        // ranges are visited using the same ContentsIterator, duplicate values
        // may be suppressed.  If you don't want this behavior, call Clear() first.
        fun startUnion(range: RangeIterator): Unit = TODO()

        // The S2CellId of the current (cell_id, label) pair.
        // REQUIRES: !done()
        fun currentCellId(): S2CellId {
            check(this::currentNode.isInitialized)
            check(!done())
            return currentNode.cellId
        }

        // The Label of the current (cell_id, label) pair.
        // REQUIRES: !done()
        fun currentLabel(): T {
            check(this::currentNode.isInitialized)
            check(!done())
            return currentNode.label
        }

        // Returns the current (cell_id, label) pair.
        // REQUIRES: !done()
        fun currentLabelledCell(): LabelledCell<T> {
            check(this::currentNode.isInitialized)
            check(!done())
            return LabelledCell(cellId = currentNode.cellId, label = currentNode.label)
        }

        // Returns true if all (cell_id, label) pairs have been visited.
        fun done(): Boolean {
            return currentNodeIdx >= cell_tree.lastIndex
        }

        // Advances the iterator to the next (cell_id, label) pair covered by the
        // current leaf cell range.
        // REQUIRES: !done()
        fun next() {
            check(!done());
            if (currentNode.parent <= node_cutoff) {
                // We have already processed this node and its ancestors.
                node_cutoff = next_node_cutoff;
                setDone()
            } else {
                currentNode = cell_tree[currentNode.parent]
            }
        }

        // node_.label == kDoneContents indicates that done() is true.
        fun setDone() {
            currentNodeIdx = cell_tree.lastIndex
        }

    }

    // A tree of (cell_id, label) pairs such that if X is an ancestor of Y, then
    // X.cell_id contains Y.cell_id.  The contents of a given range of leaf
    // cells can be represented by pointing to a node of this tree.
    private val cell_tree = mutableListOf<CellNode<T>>();

    // The last element of range_nodes_ is a sentinel value, which is necessary
    // in order to represent the range covered by the previous element.
    private val range_nodes = mutableListOf<RangeNode>()

    // Returns the number of (cell_id, label) pairs in the index.
    fun numCells(): Int = cell_tree.size

    // Adds the given (cell_id, label) pair to the index.  Note that the index
    // is not valid until Build() is called.
    //
    // The S2CellIds in the index may overlap (including duplicate values).
    // Duplicate (cell_id, label) pairs are also allowed, although be aware that
    // S2ClosestCellQuery will eliminate such duplicates anyway.
    //
    // REQUIRES: cell_id.is_valid()
    fun add(cell_id: S2CellId, label: T) {
        check(cell_id.isValid())
        cell_tree.add(CellNode(cell_id, label, -1))
    }

    // Convenience function that adds a collection of cells with the same label.
    fun add(cellIds: S2CellUnion, label: T): Unit = TODO()

    // Constructs the index.  This method may only be called once.  No iterators
    // may be used until the index is built.
    fun build() {
        TODO()
    }

    // Clears the index so that it can be re-used.
    fun clear() {
        cell_tree.clear()
        range_nodes.clear()
    }

    // Visits all (cell_id, label) pairs in the given index that intersect the
    // given S2CellUnion "target", terminating early if the given CellVisitor
    // function returns false (in which case VisitIntersectingCells returns false
    // as well).  Each (cell_id, label) pair in the index is visited at most
    // once.  (If the index contains duplicates, then each copy is visited.)
    fun visitIntersectingCells(target: S2CellUnion, visitor: CellVisitor<T>): Boolean {
        if (target.isEmpty()) return true
        val targetCellIds = target.cellIds()
        val contents = ContentsIterator()
        val range = RangeIterator()
        range.begin();
        var currentIdx: Int = 0
        var current: S2CellId = targetCellIds[currentIdx]
        do {
            if (range.limitId() <= current.rangeMin()) {
                range.seek(current.rangeMin())  // Only seek when necessary.
            }
            while (range.startId() <= current.rangeMax()) {
                contents.startUnion(range)
                while (!contents.done()) {
                if (!visitor.apply(contents.currentCellId(), contents.currentLabel())) {
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
            if (++currentIdx < targetCellIds.size) {
                current = targetCellIds[currentIdx]
                if (current.rangeMax() < range.startId()) {
                    // Skip to the first target cell that extends past the previous range.
                    currentIdx = targetCellIds.lowerBound(currentIdx + 1, targetCellIds.size, range.startId())
                    if (targetCellIds[currentIdx - 1].rangeMax() >= range.startId()) {
                        --currentIdx
                        current = targetCellIds[currentIdx]
                    }
                }
            }
        } while (currentIdx < targetCellIds.size)
        return true
    }

    // Convenience function that returns the labels of all indexed cells that
// intersect the given S2CellUnion "target".  The output contains each label
// at most once, but is not sorted.
    fun getIntersectingLabels(target: S2CellUnion): List<T> = TODO()

    // This version can be more efficient when it is called many times, since it
// does not require allocating a new vector on each call.
    fun getIntersectingLabels(target: S2CellUnion, labels: MutableList<T>): Unit = TODO()


}