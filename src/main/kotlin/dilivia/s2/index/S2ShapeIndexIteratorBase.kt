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
import dilivia.s2.S2CellId
import dilivia.s2.S2Point
import dilivia.s2.index.S2ShapeIndex
import dilivia.s2.index.S2ShapeIndexCell
import java.util.concurrent.atomic.AtomicReference

// Each subtype of S2ShapeIndex should define an Iterator type derived
// from the following base class.
abstract class S2ShapeIndexIteratorBase(): Cloneable {

    protected var id: S2CellId = S2CellId.sentinel()
    protected val cell: AtomicReference<S2ShapeIndexCell?> = AtomicReference()

    // Returns the S2CellId of the current index cell.  If done() is true,
    // returns a value larger than any valid S2CellId (S2CellId::Sentinel()).
    fun id(): S2CellId = id

    // Returns the center point of the cell.
    // REQUIRES: !done()
    fun center(): S2Point {
        Assertions.assert { !done() }
        return id().toPoint()
    }

    // Returns a reference to the contents of the current index cell.
    // REQUIRES: !done()
    fun cell(): S2ShapeIndexCell? = cell.get()

    // Returns true if the iterator is positioned past the last index cell.
    fun done(): Boolean = id == S2CellId.sentinel()

    // Positions the iterator at the first index cell (if any).
    abstract fun begin(): Unit

    // Positions the iterator past the last index cell.
    abstract fun finish(): Unit

    // Positions the iterator at the next index cell.
    // REQUIRES: !done()
    abstract fun next(): Unit

    // If the iterator is already positioned at the beginning, returns false.
    // Otherwise positions the iterator at the previous entry and returns true.
    abstract fun prev(): Boolean

    // Positions the iterator at the first cell with id() >= target, or at the
    // end of the index if no such cell exists.
    abstract fun seek(target: S2CellId): Unit

    // Positions the iterator at the cell containing "target".  If no such cell
    // exists, returns false and leaves the iterator positioned arbitrarily.
    // The returned index cell is guaranteed to contain all edges that might
    // intersect the line segment between "target" and the cell center.
    open fun locate(target: S2Point): Boolean {
        // Let I = cell_map_->lower_bound(T), where T is the leaf cell containing
        // "target_point".  Then if T is contained by an index cell, then the
        // containing cell is either I or I'.  We test for containment by comparing
        // the ranges of leaf cells spanned by T, I, and I'.

        val targetCellId = S2CellId.fromPoint(target)
        seek(targetCellId)
        if (!done() && id().rangeMin() <= targetCellId) return true
        if (prev() && id().rangeMax() >= targetCellId) return true
        return false
    }

    // Let T be the target S2CellId.  If T is contained by some index cell I
    // (including equality), this method positions the iterator at I and
    // returns INDEXED.  Otherwise if T contains one or more (smaller) index
    // cells, it positions the iterator at the first such cell I and returns
    // SUBDIVIDED.  Otherwise it returns DISJOINT and leaves the iterator
    // positioned arbitrarily.
    open fun locate(target: S2CellId): S2ShapeIndex.CellRelation {
        // Let T be the target, let I = cell_map_->lower_bound(T.range_min()), and
        // let I' be the predecessor of I.  If T contains any index cells, then T
        // contains I.  Similarly, if T is contained by an index cell, then the
        // containing cell is either I or I'.  We test for containment by comparing
        // the ranges of leaf cells spanned by T, I, and I'.
        seek(target.rangeMin())
        if (!done()) {
            if (id() >= target && id().rangeMin() <= target) return S2ShapeIndex.CellRelation.INDEXED
            if (id() <= target.rangeMax()) return S2ShapeIndex.CellRelation.SUBDIVIDED
        }
        if (prev() && id().rangeMax() >= target) return S2ShapeIndex.CellRelation.INDEXED
        return S2ShapeIndex.CellRelation.DISJOINT;
    }

    // Sets the iterator state.  "cell" typically points to the cell contents,
    // but may also be given as "nullptr" in order to implement decoding on
    // demand.  In that situation, the first that the client attempts to
    // access the cell contents, the GetCell() method is called and "cell_" is
    // updated in a thread-safe way.
    protected fun setState(id: S2CellId, cell: S2ShapeIndexCell) {
        this.id = id
        this.cell.set(cell)
    }

    // Sets the iterator state so that done() is true.
    protected fun setFinished(): Unit {
        this.id = S2CellId.sentinel()
        this.cell.set(null)
    }

    public abstract override fun clone(): S2ShapeIndexIteratorBase

}
