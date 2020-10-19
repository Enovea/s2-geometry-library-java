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
import dilivia.s2.Assertions.assertTrue
import dilivia.s2.S2CellId
import dilivia.s2.S2Point
import dilivia.s2.index.shape.CellRelation.*
import java.util.concurrent.atomic.AtomicReference
import kotlin.reflect.full.createInstance

/**
 * Base class for ShapeIndex iterators.
 * Each subtype of S2ShapeIndex should define an Iterator type derived from the following base class.
 *
 * @author Google S2Geometry Project
 * @author Fabien Meurisse (fabien.meurisse@enovea.net)
 * @since 1.0
*/
abstract class IteratorBase(): Cloneable {

    /** The current sphere cell identifier. */
    protected var id: S2CellId = S2CellId.sentinel()
    /** The contents of the current index cell. */
    private val cell: AtomicReference<S2ShapeIndexCell?> = AtomicReference()

    /**
     * Makes a copy of the given source iterator.
     */
    constructor(iter: IteratorBase): this() {
        this.id = iter.id
        this.cell.set(iter.cell())
    }

    /**
     * Gets the S2CellId of the current index cell. If done() is true, returns a value larger than any valid S2CellId
     * (S2CellId.sentinel()).
     *
     * @return The current cell.
    */
    fun id(): S2CellId = id

    /**
     * Gets the center point of the current sphere cell.
     * REQUIRES: !done()
     *
     * @return The center of the current cell.
     */
    fun center(): S2Point {
        Assertions.assert { !done() }
        return id().toPoint()
    }

    /**
     * Gets the contents of the current index cell.
     * REQUIRES: !done()
     *
     * @return The current index cell.
     */
    open fun cell(): S2ShapeIndexCell {
        // Like other const methods, this method is thread-safe provided that it
        // does not overlap with calls to non-const methods.
        assertTrue(!done())
        var currentCell = rawCell()
        if (currentCell == null) {
            currentCell = getCell()
            this.cell.set(currentCell)
        }
        return currentCell!!
    }


    /**
     * Gets the current contents of the "cell" field, which may be null if the cell contents have not been decoded yet.
     *
     * @return The current cell.
     */
    protected fun rawCell(): S2ShapeIndexCell? = cell.get()

    /**
     * Indicates if the iterator is positioned past the last index cell.
     *
     * @return true if the iterator is on the last index cell.
     */
    fun done(): Boolean = id == S2CellId.sentinel()

    /** Positions the iterator at the first index cell (if any). */
    abstract fun begin()

    /** Positions the iterator past the last index cell. */
    abstract fun finish()

    /**
     * Positions the iterator at the next index cell.
     * REQUIRES: !done()
     */
    abstract fun next()

    /**
     * If the iterator is already positioned at the beginning, returns false. Otherwise positions the iterator at the
     * previous entry and returns true.
     */
    abstract fun prev(): Boolean

    /**
     * Positions the iterator at the first cell with id() >= target, or at the end of the index if no such cell exists.
     *
     * @param target The target sphere cell.
     */
    abstract fun seek(target: S2CellId)

    /**
     * Positions the iterator at the cell containing "target".  If no such cell exists, returns false and leaves the
     * iterator positioned arbitrarily.
     * The returned index cell is guaranteed to contain all edges that might intersect the line segment between
     * "target" and the cell center.
     *
     * @param targetPoint The target point.
     * @return true if a cell containing the "target" exists and false otherwise.
     */
    open fun locate(targetPoint: S2Point): Boolean {
        // Let I = cellMap.lowerBound(T), where T is the leaf cell containing "target_point".  Then if T is
        // contained by an index cell, then the containing cell is either I or I'.  We test for containment by comparing
        // the ranges of leaf cells spanned by T, I, and I'.

        val targetCellId = S2CellId.fromPoint(targetPoint)
        seek(targetCellId)
        if (!done() && id().rangeMin() <= targetCellId) return true
        if (prev() && id().rangeMax() >= targetCellId) return true
        return false
    }

    /**
     * Let T be the target S2CellId. If T is contained by some index cell I (including equality), this method positions
     * the iterator at I and returns INDEXED.  Otherwise if T contains one or more (smaller) index cells, it positions
     * the iterator at the first such cell I and returns SUBDIVIDED. Otherwise it returns DISJOINT and leaves the
     * iterator positioned arbitrarily.
     *
     * @param target The target sphere cell.
     */
    open fun locate(target: S2CellId): CellRelation {
        // Let T be the target, let I = cellMap.lowerBound(T.rangeMin()), and  let I' be the predecessor of I.
        // If T contains any index cells, then T contains I. Similarly, if T is contained by an index cell, then the
        // containing cell is either I or I'.  We test for containment by comparing the ranges of leaf cells spanned
        // by T, I, and I'.
        seek(target.rangeMin())
        if (!done()) {
            if (id() >= target && id().rangeMin() <= target) return INDEXED
            if (id() <= target.rangeMax()) return SUBDIVIDED
        }
        if (prev() && id().rangeMax() >= target) return INDEXED
        return DISJOINT;
    }

    /**
     * Sets the iterator state. "cell" typically points to the cell contents, but may also be given as "null" in order
     * to implement decoding on demand.  In that situation, the first that the client attempts to access the cell
     * contents, the GetCell() method is called and "cell_" is updated in a thread-safe way.
     *
     * @param id The current sphere cell.
     * @param cell The current index cell.
     */
    protected fun setState(id: S2CellId, cell: S2ShapeIndexCell?) {
        this.id = id
        this.cell.set(cell)
    }

    /**
     * This method is called to decode the contents of the current cell, if setState() was previously called with a
     * null "cell" argument.  This allows decoding on demand for subtypes that keep the cell contents in an encoded
     * state. It does not need to be implemented at all if setState() is always called with (cell != nullptr).
     *
     * REQUIRES: This method is thread-safe.
     * REQUIRES: Multiple calls to this method return the same value.
     */
    protected open fun getCell(): S2ShapeIndexCell? = null

    /** Sets the iterator state so that done() is true. */
    protected fun setFinished() {
        this.id = S2CellId.sentinel()
        this.cell.set(null)
    }

    public override fun clone(): IteratorBase {
        val clone = this::class.createInstance()
        clone.setState(id, cell())
        return clone
    }

    override fun equals(other: Any?): Boolean {
        if (this === other) return true
        if (other !is IteratorBase) return false

        if (id != other.id) return false
        if (cell.get() != other.cell.get()) return false

        return true
    }

    override fun hashCode(): Int {
        var result = id.hashCode()
        result = 31 * result + cell.hashCode()
        return result
    }


}
