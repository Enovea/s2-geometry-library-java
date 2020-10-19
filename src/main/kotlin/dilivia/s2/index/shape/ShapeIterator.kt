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

import dilivia.s2.shape.S2Shape

/**
 * Allows iterating over the indexed shapes.
 * CAVEAT: Returns null for shapes that have been removed from the index.
 *
 * @author Google S2Geometry Project
 * @author Fabien Meurisse (fabien.meurisse@enovea.net)
 * @since 1.0
 */
class ShapeIterator(val index: S2ShapeIndex) {

    /** Current shape id. */
    private var shapeId: Int = 0

    internal constructor(index: S2ShapeIndex, shapeId: Int): this(index) {
        this.shapeId = shapeId
    }

    /**
     * Gets the shape at the current position, null if the shape has been removed.
     *
     * @return The current shape.
     */
    fun shape(): S2Shape? = index.shape(shapeId)

    /**
     * Increments the iterator.
     *
     * @return A new iterator instance that points to the next shape id.
     */
    operator fun inc(): ShapeIterator = ShapeIterator(index = index, shapeId = shapeId + 1)

    fun next(): ShapeIterator { ++shapeId; return this }

    /**
     * Decrements the iterator.
     *
     * @return A new iterator instance that points to the previous shape id.
     */
    operator fun dec(): ShapeIterator = ShapeIterator(index = index, shapeId = shapeId - 1)

    fun previous(): ShapeIterator { --shapeId; return this }

    override fun equals(other: Any?): Boolean {
        if (this === other) return true
        if (other !is ShapeIterator) return false

        if (index != other.index) return false
        if (shapeId != other.shapeId) return false

        return true
    }

    override fun hashCode(): Int {
        var result = index.hashCode()
        result = 31 * result + shapeId
        return result
    }

}
