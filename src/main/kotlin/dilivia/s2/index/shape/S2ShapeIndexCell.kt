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
 * S2ShapeIndexCell stores the index contents for a particular S2CellId. It consists of a set of clipped shapes.
 *
 * @author Google S2Geometry Project
 * @author Fabien Meurisse (fabien.meurisse@enovea.net)
 * @since 1.0
 */
class S2ShapeIndexCell {

    /** Contents of the cell. */
    private val shapes = mutableListOf<S2ClippedShape>()

    /** The number of clipped shapes in this cell. */
    val numClipped: Int
        get() = shapes.size

    /**
     * Gets the clipped shape at the given index. Shapes are kept sorted in increasing order of shape id.
     *
     * REQUIRES: 0 <= i < numClipped)
     * @param i Index of the clipped shape.
     * @return The clipped shape at the given index.
     */
    fun clipped(i: Int): S2ClippedShape = shapes[i]

    /**
     * Gets the clipped shape corresponding to the given shape, or null if it does not intersect this cell.
     *
     * @param shape A shape.
     * @return The corresponding clipped shape if it exists.
     */
    fun findClipped(shape: S2Shape): S2ClippedShape? = findClipped(shape.id)


    /**
     * Gets the clipped shape corresponding to the given shape specified by its id, or null if it does not intersect
     * this cell.
     *
     * @param shapeId A shape id.
     * @return The corresponding clipped shape if it exists.
     */
    fun findClipped(shapeId: Int): S2ClippedShape? {
        // Linear search is fine because the number of shapes per cell is typically very small (most often 1), and is
        // large only for pathological inputs (e.g. very deeply nested loops).
        for (s in shapes) {
            if (s.shapeId == shapeId) return s
        }
        return null
    }

    /**
     * Convenience method that returns the total number of edges in all clipped shapes.
     *
     * @return The total number of contained edges.
     */
    fun numEdges(): Int = shapes.map { clipped -> clipped.numEdges }.sum()

    /**
     * Adds a clipped shape in this cell.
     *
     * @param shape The clipped shape to add.
     */
    fun addClipped(shape: S2ClippedShape) {
        shapes.add(shape)
    }

    /**
     * Returns a string representation of this cell.
     */
    override fun toString(): String {
        return "S2ShapeIndexCell(shapes=$shapes)"
    }

    fun toDebugString(separator: String = ", ", prefix: String = "  -"): String =
            shapes.joinToString(separator = separator) { shape ->
                "$prefix${shape.shapeId}, ${shape.containsCenter}, ${shape.edges})"
            }

}
