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

import dilivia.s2.index.shape.S2ShapeIndex

object S2CountEdges {

    // Returns the total number of edges in all indexed shapes.  This method takes
    // time linear in the number of shapes.
    fun countEdges(index: S2ShapeIndex): Int = countEdgesUpTo(index, Int.MAX_VALUE)

    // Like CountEdges(), but stops once "max_edges" edges have been found (in
    // which case the current running total is returned).
    fun countEdgesUpTo(index: S2ShapeIndex, max_edges: Int): Int {
        var numEdges = 0
        index.shapeIterator().asSequence().filterNotNull().forEach { shape ->
            numEdges += shape.numEdges
            if (numEdges >= max_edges) return numEdges
        }
        return numEdges
    }


}
