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
package dilivia.s2.index.point

import com.google.common.collect.SortedMultiset
import com.google.common.collect.TreeMultiset
import dilivia.s2.S2CellId
import dilivia.s2.S2Point
import mu.KotlinLogging
import java.util.*

/**
 * S2PointIndex maintains an index of points sorted by leaf S2CellId.
 *
 * Each point can optionally store auxiliary data such as an integer or object. This can be used to map results back
 * to client data structures.
 *
 * The class supports adding or removing points dynamically, and provides a seekable iterator interface for navigating
 * the index.
 *
 *  You can use this class in conjuction with S2ClosestPointQuery to find the closest index points to a given query
 *  point.  For example:
 *
 *  <pre>
 * fun test(indexPoints: List<S2Point>, target_points: List<S2Point) {
 *   // The template argument allows auxiliary data to be attached to each
 *   // point (in this case, the array index).
 *   val index = S2PointIndex<Int>()
 *   indexPoints.forEachIndexed { i, p -> index.add(p, i) }
 *
 *   val query = S2ClosestPointQuery<Int>(index)
 *   query.mutable_options()->set_max_results(5);
 *   for (const S2Point& target_point : target_points) {
 *     S2ClosestPointQueryPointTarget target(target_point);
 *     for (const auto& result : query.FindClosestPoints(&target)) {
 *       // The Result class contains the following methods:
 *       //   distance() is the distance to the target.
 *       //   point() is the indexed point.
 *       //   data() is the auxiliary data.
 *       DoSomething(target_point, result);
 *     }
 *   }
 * }
 * </pre>
 *
 *
 * Points can be added or removed from the index at any time by calling add() or remove(). However when the index is
 * modified, you must call init() on each iterator before using it again (or simply create a new iterator).
 *
 * <pre>
 *   index.add(newPoint, 123456)
 *   it.init(index)
 *   it.seek(target.rangeMin())
 * </pre>
 *
 * You can also access the index directly using the iterator interface. For example, here is how to iterate through
 * all the points in a given S2CellId "targetId":
 *
 * <pre>
 *   val it = S2PointIndexIterator<Int>(index)
 *   it.seek(targetId.rangeMin());
 *   while (!it.done() && it.id() <= targetId.rangeMax()) {
 *     doSomething(it.id(), it.point(), it.data())
 *     it.next()
 *   }
 * </pre>
 *
 * TODO(ericv): Consider adding an S2PointIndexRegion class, which could be used to efficiently compute coverings of a collection of S2Points.
 */
class S2PointIndex<T : Comparable<T>> {

    internal val map: TreeMap<S2CellId, SortedMultiset<PointData<T>>> = TreeMap()

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

    fun iterator(): S2PointIndexIterator<T> = S2PointIndexIterator(this)

    companion object {
        private val logger = KotlinLogging.logger(S2PointIndex::class.java.name)
    }

}
