package dilivia.s2.index.point

import dilivia.s2.S2CellId
import dilivia.s2.S2Point
import mu.KotlinLogging


class S2PointIndexIterator<T : Comparable<T>>(index: S2PointIndex<T>) {

    private lateinit var index: S2PointIndex<T>
    private var currentCellId: S2CellId? = null
    private var currentPointData: PointData<T>? = null
    private var currentOccurence: Int = 0

    init {
        init(index)
    }

    // Initializes an iterator for the given S2PointIndex.  If the index is
    // non-empty, the iterator is positioned at the first cell.
    //
    // This method may be called multiple times, e.g. to make an iterator
    // valid again after the index is modified.
    fun init(index: S2PointIndex<T>) {
        this.index = index
        begin()
    }

    // The S2CellId for the current index entry.
    // REQUIRES: !done()
    fun id(): S2CellId = currentCellId!!

    // The point associated with the current index entry.
    // REQUIRES: !done()
    fun point(): S2Point = currentPointData!!.point

    // The client-supplied data associated with the current index entry.
    // REQUIRES: !done()
    fun data(): T = currentPointData!!.data

    // The (S2Point, data) pair associated with the current index entry.
    fun pointData(): PointData<T> = currentPointData!!

    // Returns true if the iterator is positioned past the last index entry.
    fun done(): Boolean = currentPointData == null

    // Positions the iterator at the first index entry (if any).
    fun begin() {
        currentPointData = null
        currentCellId = if (index.map.isNotEmpty()) index.map.firstKey() else null
        while (currentPointData == null && currentCellId != null) {
            currentPointData = index.map[currentCellId]?.firstOrNull()
            if (currentPointData == null) {
                currentCellId = index.map.higherKey(currentCellId)
            }
        }
        if (currentPointData != null) currentOccurence = 1

        logger.trace { """
                |Iterator.begin()
                |--------------------------
                | Current cell id: $currentCellId
                | Current point data: $currentPointData
                | Current occurence: $currentOccurence
            """.trimMargin() }
    }

    // Positions the iterator so that done() is true.
    fun finish() {
        currentCellId = null
        currentPointData = null
        currentOccurence = 0

        logger.trace { """
                |Iterator.finish()
                |--------------------------
                | Current cell id: $currentCellId
                | Current point data: $currentPointData
                | Current occurence: $currentOccurence
            """.trimMargin() }
    }

    // Advances the iterator to the next index entry.
    // REQUIRES: !done()
    fun next() {
        var nextPointData: PointData<T>? = null
        var cellId = currentCellId
        var nextOccurence = currentOccurence + 1
        while (nextPointData == null && cellId != null) {
            val pointMultiset = index.map[cellId]
            if (pointMultiset != null) {
                if (nextOccurence > pointMultiset.count(currentPointData)) {
                    nextPointData = pointMultiset.elementSet()?.higher(currentPointData)
                    nextOccurence = 1
                } else {
                    nextPointData = currentPointData
                }
            }
            if (nextPointData == null) {
                cellId = index.map.higherKey(cellId)
                nextPointData = cellId?.let { index.map[it]?.firstOrNull() }
                nextOccurence = if (nextPointData == null) 0 else 1
            }
        }
        currentCellId = cellId
        currentPointData = nextPointData
        currentOccurence = nextOccurence

        logger.trace { """
                |Iterator.next()
                |--------------------------
                | Current cell id: $currentCellId
                | Current point data: $currentPointData
                | Current occurence: $currentOccurence
            """.trimMargin() }

    }

    // If the iterator is already positioned at the beginning, returns false.
    // Otherwise positions the iterator at the previous entry and returns true.
    fun prev(): Boolean {
        var previousPointData: PointData<T>? = null
        var cellId = currentCellId
        var previousOccurence = currentOccurence - 1
        while (previousPointData == null && cellId != null) {
            var pointMultiset = index.map[cellId]
            if (pointMultiset != null) {
                if (previousOccurence <= 0) {
                    previousPointData = index.map[cellId]?.elementSet()?.lower(currentPointData)
                    previousOccurence = pointMultiset.count(previousPointData)
                } else {
                    previousPointData = currentPointData
                }
            }
            if (previousPointData == null) {
                cellId = index.map.lowerKey(cellId)
                pointMultiset = cellId?.let { index.map[it] }
                previousPointData = pointMultiset?.lastOrNull()
                previousOccurence = pointMultiset?.size ?: 0
            }
        }
        val result = if (previousPointData != null) {
            currentCellId = cellId
            currentPointData = previousPointData
            currentOccurence = previousOccurence
            true
        } else false

        logger.trace { """
                |Iterator.prev()
                |--------------------------
                | Current cell id: $currentCellId
                | Current point data: $currentPointData
                | Current occurence: $currentOccurence
                | Has moved: $result
            """.trimMargin() }

        return result
    }

    // Positions the iterator at the first entry with id() >= target, or at the
    // end of the index if no such entry exists.
    fun seek(target: S2CellId) {
        currentCellId = index.map.ceilingKey(target)
        currentPointData = null
        currentOccurence = 0
        while (currentPointData == null && currentCellId != null) {
            val pointList = index.map.getValue(currentCellId!!)
            if (pointList.isNotEmpty()) {
                currentPointData = pointList.first()
                currentOccurence = 1
            }
            else {
                currentCellId = index.map.higherKey(currentCellId)
            }
        }


        logger.trace { """
                |Iterator.seek($target)
                |--------------------------
                | Current cell id: $currentCellId
                | Current point data: $currentPointData
                | Current occurence: $currentOccurence
            """.trimMargin() }
    }

    companion object {
        private val logger = KotlinLogging.logger(S2PointIndexIterator::class.java.name)
    }

}
