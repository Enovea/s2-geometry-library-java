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

import dilivia.s2.S2GeometryTestCase
import dilivia.s2.S2TextParser
import dilivia.s2.index.shape.S2CountEdges

class S2CountEdgesTest : S2GeometryTestCase() {

fun testCountEdgesUpToStopsEarly() {
  val index = S2TextParser.makeIndex(
      "0:0 | 0:1 | 0:2 | 0:3 | 0:4 # 1:0, 1:1 | 1:2, 1:3 | 1:4, 1:5, 1:6 #"
  )
  // Verify the test parameters.
  assertEquals(index.numShapeIds(), 4)
  assertEquals(index.shape(0)?.numEdges, 5)
  assertEquals(index.shape(1)?.numEdges, 1)
  assertEquals(index.shape(2)?.numEdges, 1)
  assertEquals(index.shape(3)?.numEdges, 2)

  assertEquals(S2CountEdges.countEdges(index), 9)
  assertEquals(S2CountEdges.countEdgesUpTo(index, 1), 5)
  assertEquals(S2CountEdges.countEdgesUpTo(index, 5), 5)
  assertEquals(S2CountEdges.countEdgesUpTo(index, 6), 6)
  assertEquals(S2CountEdges.countEdgesUpTo(index, 8), 9)
}

}
