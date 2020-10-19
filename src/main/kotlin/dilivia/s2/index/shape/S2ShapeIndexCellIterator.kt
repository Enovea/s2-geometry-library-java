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
package dilivia.s2.index.shape

import dilivia.s2.Assertions
import dilivia.s2.S2CellId
import dilivia.s2.S2Point
import dilivia.s2.index.shape.InitialPosition.UNPOSITIONED

/**
 * A random access iterator that provides low-level access to the cells of the index.
 * Cells are sorted in increasing order of S2CellId.
 *
 * This shape index iterator delegates the iteration action to a concrete iterator base instance associated with the
 * shape index type we are scanning.
 *
 * @constructor Default constructor; must be followed by a call to Init().
 *
 * @author Google S2Geometry Project
 * @author Fabien Meurisse (fabien.meurisse@enovea.net)
 * @since 1.0
 */
class S2ShapeIndexCellIterator(): Cloneable {

    /** The concrete iterator. */
    private lateinit var iter: IteratorBase

    /**
     * Constructs an iterator positioned as specified. By default iterators are unpositioned, since this avoids an
     * extra seek in this situation where one of the seek methods (such as locate) is immediately called.
     *
     * If you want to position the iterator at the beginning, e.g. in order to loop through the entire index,
     * do this:
     *
     * <pre>
     *     val iter = S2ShapeIndexCellIterator(index, BEGIN)
     *     while(!iter.done() {
     *          ...
     *          iter.next()
     *     }
     * </pre>
     *
     * @param index The shape index.
     * @param pos The initial position of the iterator (Default = UNPOSITIONED).
     */
    constructor(index: S2ShapeIndex, pos: InitialPosition = UNPOSITIONED): this() {
        this.iter = index.newIterator(pos)
    }

    /**
     * Makes a copy of the given source iterator.
     *
     * @param cellIter The iterator to copy.
     */
    constructor(cellIter: S2ShapeIndexCellIterator): this() {
        this.iter = cellIter.iter.clone()
    }

    /**
     * Initializes an iterator for the given S2ShapeIndex.  This method may also be called in order to restore an
     * iterator to a valid state after the underlying index has been updated (although it is usually easier just to
     * declare a new iterator whenever required, since iterator construction is cheap).
     *
     * @param index The shape index.
     * @param pos The initial position of the iterator (Default = UNPOSITIONED).
     */
    fun init(index: S2ShapeIndex, pos: InitialPosition = UNPOSITIONED) {
        this.iter = index.newIterator(pos)
    }

    /**
     * Gets the S2CellId of the current index cell. If done() is true, returns a value larger than any valid S2CellId
     * (S2CellId.sentinel()).
     *
     * @return The current sphere cell id.
     */
    fun id(): S2CellId = iter.id()

    /**
     * Gets the center point of the sphere cell.
     *
     * REQUIRES: !done()
     * @return The center point.
     */
    fun center(): S2Point {
        Assertions.assert { !done() }
        return id().toPoint()
    }

    /**
     * Gets the contents of the current index cell.
     *
     * REQUIRES: !done()
     * @return the contents of the current cell.
     */
    fun cell(): S2ShapeIndexCell {
        Assertions.assert { !done() }
        return iter.cell()
    }

    /**
     * Indicates if the iterator is positioned past the last index cell.
     *
     * @return true if the iterator is at the end and false otherwise.
     */
    fun done(): Boolean { return iter.done() }

    /** Positions the iterator at the first index cell (if any). */
    fun begin() { iter.begin() }

    /** Positions the iterator past the last index cell. */
    fun finish() { iter.finish() }

    /**
     * Positions the iterator at the next index cell.
     * REQUIRES: !done()
     */
    fun next() { iter.next() }

    /**
     * If the iterator is already positioned at the beginning, returns false.
     * Otherwise positions the iterator at the previous entry and returns true.
     *
     * @return true if the iterator have been moved.
     */
    fun prev(): Boolean = iter.prev()

    /**
     * Positions the iterator at the first cell with id() >= target, or at the end of the index if no such cell exists.
     *
     * @param target The target cell.
     */
    fun seek(target: S2CellId) { iter.seek(target) }

    /**
     * Positions the iterator at the cell containing "target". If no such cell exists, returns false and leaves the
     * iterator positioned arbitrarily. The returned index cell is guaranteed to contain all edges that might intersect
     * the line segment between "target" and the cell center.
     *
     * @param target The target point
     * @return true if the iterator has been moved to target point, false otherwise.
     */
    fun locate(target: S2Point): Boolean = iter.locate(target)

    /**
     * Moves the iterator to a target cell.
     *
     * Let T be the target S2CellId.
     * - If T is contained by some index cell I (including equality), this method positions the iterator at I and
     *   returns INDEXED.
     * - Otherwise if T contains one or more (smaller) index cells, it positions the iterator at the first such cell I
     *   and returns SUBDIVIDED.
     * - Otherwise it returns DISJOINT and leaves the iterator positioned arbitrarily.
     *
     * @param target The target cell.
     * @return The cell relation between the target and the underlying index.
     */
    fun locate(target: S2CellId): CellRelation = iter.locate(target)

    public override fun clone(): S2ShapeIndexCellIterator {
        val cellIter = S2ShapeIndexCellIterator()
        cellIter.iter = iter.clone()
        return cellIter
    }
}
