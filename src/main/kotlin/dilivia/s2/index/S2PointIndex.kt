package dilivia.s2.index

import com.google.common.collect.SortedMultiset
import com.google.common.collect.TreeMultiset
import dilivia.s2.S2CellId
import dilivia.s2.S2Point
import mu.KotlinLogging
import java.util.*

// S2PointIndex maintains an index of points sorted by leaf S2CellId.  Each
// point can optionally store auxiliary data such as an integer or pointer.
// This can be used to map results back to client data structures.
//
// The class supports adding or removing points dynamically, and provides a
// seekable iterator interface for navigating the index.
//
// You can use this class in conjuction with S2ClosestPointQuery to find the
// closest index points to a given query point.  For example:
//
// void Test(const vector<S2Point>& index_points,
//           const vector<S2Point>& target_points) {
//   // The template argument allows auxiliary data to be attached to each
//   // point (in this case, the array index).
//   S2PointIndex<int> index;
//   for (int i = 0; i < index_points.size(); ++i) {
//     index.Add(index_points[i], i);
//   }
//   S2ClosestPointQuery<int> query(&index);
//   query.mutable_options()->set_max_results(5);
//   for (const S2Point& target_point : target_points) {
//     S2ClosestPointQueryPointTarget target(target_point);
//     for (const auto& result : query.FindClosestPoints(&target)) {
//       // The Result class contains the following methods:
//       //   distance() is the distance to the target.
//       //   point() is the indexed point.
//       //   data() is the auxiliary data.
//       DoSomething(target_point, result);
//     }
//   }
// }
//
// The Data argument defaults to an empty class, which uses no additional
// space beyond the S2Point itself.  In this case the Data argument is
// required.  For example:
//
//   S2PointIndex<> index;
//   index.Add(point);
//
// Points can be added or removed from the index at any time by calling Add()
// or Remove().  However when the index is modified, you must call Init() on
// each iterator before using it again (or simply create a new iterator).
//
//   index.Add(new_point, 123456);
//   it.Init(&index);
//   it.Seek(target.range_min());
//
// You can also access the index directly using the iterator interface.  For
// example, here is how to iterate through all the points in a given S2CellId
// "target_id":
//
//   S2PointIndex<int>::Iterator it(&index);
//   it.Seek(target_id.range_min());
//   for (; !it.done() && it.id() <= target_id.range_max(); it.Next()) {
//     DoSomething(it.id(), it.point(), it.data());
//   }
//
// TODO(ericv): Consider adding an S2PointIndexRegion class, which could be
// used to efficiently compute coverings of a collection of S2Points.
//
// REQUIRES: "Data" has default and copy constructors.
// REQUIRES: "Data" has operator== and operator<.
interface Data<T : Data<T>> : Comparable<T>

// PointData is essentially std::pair with named fields.  It stores an
// S2Point and its associated data, taking advantage of the "empty base
// optimization" to ensure that no extra space is used when Data is empty.
data class PointData<T : Comparable<T>>(val point: S2Point, val data: T) : Data<PointData<T>> {

    override fun compareTo(other: PointData<T>): Int {
        val pointComparison = point.compareTo(other.point)
        return if (pointComparison == 0) data.compareTo(other.data) else pointComparison
    }

}

class S2PointIndex<T : Comparable<T>>() {

    private val map: TreeMap<S2CellId, SortedMultiset<PointData<T>>> = TreeMap()

    // Returns the number of points in the index.
    fun numPoints(): Int = map.map { entry -> entry.value.size }.sum()

    // Adds the given point to the index.  Invalidates all iterators.
    fun add(point: S2Point, data: T) = add(PointData(point, data))
    fun add(point_data: PointData<T>) {
        val id = S2CellId.fromPoint(point_data.point)
        map.getOrPut(id, { TreeMultiset.create() }).add(point_data)
    }

    // Removes the given point from the index.  Both the "point" and "data"
    // fields must match the point to be removed.  Returns false if the given
    // point was not present.  Invalidates all iterators.
    fun remove(point: S2Point, data: T): Boolean = remove(PointData(point, data))
    fun remove(point_data: PointData<T>): Boolean {
        val id = S2CellId.fromPoint(point_data.point)
        val dataSet = map[id] ?: return false
        val removed = dataSet.remove(point_data)
        if (removed && dataSet.isEmpty()) {
            map.remove(id)
        }
        return removed
    }

    // Resets the index to its original empty state.  Invalidates all iterators.
    fun clear(): Unit = map.clear()

    fun iterator(): Iterator = Iterator()

    inner class Iterator() {

        private var currentCellId: S2CellId? = null
        private var currentPointData: PointData<T>? = null
        private var currentOccurence: Int = 0

        /*
              private:
              typename Map::const_iterator iter_, end_;

         */
        init {
            init()
        }

        // Initializes an iterator for the given S2PointIndex.  If the index is
        // non-empty, the iterator is positioned at the first cell.
        //
        // This method may be called multiple times, e.g. to make an iterator
        // valid again after the index is modified.
        fun init(): Unit {
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
            currentCellId = if (map.isNotEmpty()) map.firstKey() else null
            while (currentPointData == null && currentCellId != null) {
                currentPointData = map[currentCellId]?.firstOrNull()
                if (currentPointData == null) {
                    currentCellId = map.higherKey(currentCellId)
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
                val pointMultiset = map[cellId]
                if (pointMultiset != null) {
                    if (nextOccurence > pointMultiset.count(currentPointData)) {
                        nextPointData = pointMultiset.elementSet()?.higher(currentPointData)
                        nextOccurence = 1
                    } else {
                        nextPointData = currentPointData
                    }
                }
                if (nextPointData == null) {
                    cellId = map.higherKey(cellId)
                    nextPointData = cellId?.let { map[it]?.firstOrNull() }
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
                var pointMultiset = map[cellId]
                if (pointMultiset != null) {
                    if (previousOccurence <= 0) {
                        previousPointData = map[cellId]?.elementSet()?.lower(currentPointData)
                        previousOccurence = pointMultiset.count(previousPointData)
                    } else {
                        previousPointData = currentPointData
                    }
                }
                if (previousPointData == null) {
                    cellId = map.lowerKey(cellId)
                    pointMultiset = cellId?.let { map[it] }
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
        fun seek(target: S2CellId): Unit {
            currentCellId = map.ceilingKey(target)
            currentPointData = null
            currentOccurence = 0
            while (currentPointData == null && currentCellId != null) {
                val pointList = map.getValue(currentCellId!!)
                if (pointList.isNotEmpty()) {
                    currentPointData = pointList.first()
                    currentOccurence = 1
                }
                else {
                    currentCellId = map.higherKey(currentCellId)
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

    }

    companion object {
        private val logger = KotlinLogging.logger(S2PointIndex::class.java.name)
    }

}
