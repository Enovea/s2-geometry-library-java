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

import dilivia.s2.Fractal
import dilivia.s2.S2Random
import dilivia.s2.math.Matrix3x3
import dilivia.s2.region.S2Cap
import dilivia.s2.region.S2Loop

// Generates a fractal loop that approximately fills the given S2Cap.
class FractalLoopShapeIndexFactory : ShapeIndexFactory {

    override fun addEdges(index_cap: S2Cap, num_edges: Int, index: MutableS2ShapeIndex) {
        val fractal = Fractal()
        fractal.setLevelForApproxMaxEdges(num_edges)
        index.add(S2Loop.Shape(loop = fractal.makeLoop(
                Matrix3x3.fromCols(S2Random.randomFrameAt(index_cap.center)),
                index_cap.radius()
        )))
    }

}
