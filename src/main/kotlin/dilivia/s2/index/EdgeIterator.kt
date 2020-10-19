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
import dilivia.s2.index.shape.S2ShapeIndex
import dilivia.s2.shape.Edge
import dilivia.s2.shape.ShapeEdgeId

// An iterator that advances through all edges in an S2ShapeIndex.
//
// Example usage:
//
// for (EdgeIterator it(index); !it.Done(); it.Next()) {
//   auto edge = it.edge();
//   //...
// }
class EdgeIterator {

    val index: S2ShapeIndex
    private var shape_id = -1
    private var num_edges = 0
    private var edge_id = -1

    constructor(index: S2ShapeIndex) {
        this.index = index
        next()
    }

    constructor(iterator: EdgeIterator) {
        this.index = iterator.index
        this.shape_id = iterator.shape_id
        this.num_edges = iterator.num_edges
        this.edge_id = iterator.edge_id
    }

    // Returns the current shape id.
    fun shapeId() = shape_id

    // Returns the current edge id.
    fun edgeId() = edge_id

    // Returns the current (shape_id, edge_id).
    fun shape_edge_id(): ShapeEdgeId = ShapeEdgeId(shape_id, edge_id)

    // Returns the current edge.
    fun edge(): Edge {
        Assertions.assert { !done() }
        return index.shape(shape_id)!!.edge(edge_id)
    }

    // Returns true if there are no more edges in the index.
    fun done() = shapeId() >= index.numShapeIds()

    // Advances to the next edge.
    fun next(): Unit {
        while (++edge_id >= num_edges) {
            if (++shape_id >= index.numShapeIds()) break
            val shape = index.shape(shape_id)
            num_edges = shape?.numEdges ?: 0
            edge_id = -1
        }
    }

    fun debugString(): String = "(shape=$shape_id, edge=$edge_id)"

    override fun equals(other: Any?): Boolean {
        if (this === other) return true
        if (other !is EdgeIterator) return false

        if (index != other.index) return false
        if (shape_id != other.shape_id) return false
        if (edge_id != other.edge_id) return false

        return true
    }

    override fun hashCode(): Int {
        var result = index.hashCode()
        result = 31 * result + shape_id
        result = 31 * result + edge_id
        return result
    }

}

