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

import dilivia.s2.shape.S2ClippedShape
import dilivia.s2.shape.S2Shape

// S2ShapeIndexCell stores the index contents for a particular S2CellId.
// It consists of a set of clipped shapes.
class S2ShapeIndexCell() {

    private val shapes = mutableListOf<S2ClippedShape>()

    // Returns the number of clipped shapes in this cell.
    fun numClipped(): Int {
        return shapes.size
    }

    // Returns the clipped shape at the given index.  Shapes are kept sorted in
    // increasing order of shape id.
    //
    // REQUIRES: 0 <= i < num_clipped()
    fun clipped(i: Int): S2ClippedShape {
        return shapes[i]
    }

    // Returns a pointer to the clipped shape corresponding to the given shape,
    // or nullptr if the shape does not intersect this cell.
    fun findClipped(shape: S2Shape): S2ClippedShape? = findClipped(shape.id)
    fun findClipped(shape_id: Int): S2ClippedShape? {
        // Linear search is fine because the number of shapes per cell is typically
        // very small (most often 1), and is large only for pathological inputs
        // (e.g. very deeply nested loops).
        for (s in shapes) {
            if (s.shapeId == shape_id) return s
        }
        return null
    }

    // Convenience method that returns the total number of edges in all clipped
    // shapes.
    fun num_edges(): Int {
        var n = 0
        for (i in 0 until numClipped()) n += clipped(i).numEdges()
        return n
    }

    fun addClipped(shape: S2ClippedShape) {
        shapes.add(shape)
    }

    override fun toString(): String {
        return "S2ShapeIndexCell(shapes=$shapes)"
    }


}