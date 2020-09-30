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
package dilivia.s2.shape

import dilivia.s2.Assertions
import dilivia.s2.Assertions.assertEQ
import dilivia.s2.S2EdgeCrosser
import dilivia.s2.S2Error
import dilivia.s2.S2Point

object S2ShapeUtil {



    // Returns true if the given shape contains the given point.  Most clients
    // should not use this method, since its running time is linear in the number
    // of shape edges.  Instead clients should create an S2ShapeIndex and use
    // S2ContainsPointQuery, since this strategy is much more efficient when many
    // points need to be tested.
    //
    // Polygon boundaries are treated as being semi-open (see S2ContainsPointQuery
    // and S2VertexModel for other options).
    //
    // CAVEAT: Typically this method is only used internally.  Its running time is
    //         linear in the number of shape edges.
    fun containsBruteForce(shape: S2Shape, focus: S2Point): Boolean {
        if (shape.dimension < 2) return false

        val refPoint = shape.getReferencePoint()
        if (refPoint.point == focus) return refPoint.contained

        val crosser = S2EdgeCrosser(refPoint.point, focus);
        var inside = refPoint.contained;
        for (e in 0 until  shape.numEdges) {
            val edge = shape.edge(e)
            inside = inside xor crosser.edgeOrVertexCrossing(edge.v0, edge.v1)
        }
        return inside;
    }

}