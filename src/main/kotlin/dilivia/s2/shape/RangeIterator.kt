package dilivia.s2.shape

import dilivia.s2.S2CellId

// RangeIterator is a wrapper over S2ShapeIndex::Iterator with extra methods
// that are useful for merging the contents of two or more S2ShapeIndexes.
class RangeIterator(index: S2ShapeIndex) {

    private val iterator: S2ShapeIndexIteratorBase = index.iterator(S2ShapeIndex.InitialPosition.BEGIN)
    // The min and max leaf cell ids covered by the current cell.  If done() is
    // true, these methods return a value larger than any valid cell id.
    private var rangeMin: S2CellId = S2CellId.sentinel()
    private var rangeMax: S2CellId = S2CellId.sentinel()
    
    init {
        refresh()
    }
    
    // The current S2CellId and cell contents.
    fun id(): S2CellId = iterator.id()
    fun cell(): S2ShapeIndexCell? = iterator.cell()

    // The min and max leaf cell ids covered by the current cell.  If done() is
    // true, these methods return a value larger than any valid cell id.
    fun rangeMin(): S2CellId = rangeMin
    fun rangeMax(): S2CellId = rangeMax

    // Various other convenience methods for the current cell.
    fun clipped(i: Int): S2ClippedShape = cell()!!.clipped(i)
    fun numEdges(i: Int): Int = clipped(i).numEdges()
    fun containsCenter(i: Int): Boolean = clipped(i).containsCenter

    fun next() {
        iterator.next()
        refresh()
    }
    
    fun done(): Boolean = iterator.done()

    // Position the iterator at the first cell that overlaps or follows
    // "target", i.e. such that range_max() >= target.range_min().
    fun seekTo(target: RangeIterator) {
        iterator.seek(target.rangeMin)
        // If the current cell does not overlap "target", it is possible that the
        // previous cell is the one we are looking for.  This can only happen when
        // the previous cell contains "target" but has a smaller S2CellId.
        if (iterator.done() || iterator.id().rangeMin() > target.rangeMax) {
            if (iterator.prev() && iterator.id().rangeMax() < target.id()) iterator.next()
        }
        refresh()
    }

    // Position the iterator at the first cell that follows "target", i.e. the
    // first cell such that range_min() > target.range_max().
    fun seekBeyond(target: RangeIterator) {
        iterator.seek(target.rangeMax.next())
        if (!iterator.done() && iterator.id().rangeMin() <= target.rangeMax) {
            iterator.next()
        }
        refresh()
    }

    // Updates internal state after the iterator has been repositioned.
    private fun refresh() {
        rangeMin = id().rangeMin()
        rangeMax = id().rangeMax()
    }
}