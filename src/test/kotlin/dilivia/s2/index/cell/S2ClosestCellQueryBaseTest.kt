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
package dilivia.s2.index.cell

import dilivia.s2.S1ChordAngle
import dilivia.s2.S2GeometryTestCase
import dilivia.s2.S2Point
import dilivia.s2.S2TextParser.makeCellId
import dilivia.s2.S2TextParser.makeCellUnion
import dilivia.s2.index.S2CellIndex
import dilivia.s2.index.S2ClosestCellQueryBase
import dilivia.s2.index.S2MaxDistance
import dilivia.s2.index.S2MaxDistanceFactory
import dilivia.s2.index.S2MaxDistancePointTarget

// This is a proof-of-concept prototype of a possible S2FurthestCellQuery
// class.  The purpose of this test is just to make sure that the code
// compiles and does something reasonable.
typealias FurthestCellQuery = S2ClosestCellQueryBase<S2MaxDistance>

class FurthestPointTarget(point: S2Point) : S2MaxDistancePointTarget(point) {

  override fun maxBruteForceIndexSize(): Int {
    return 10
  }

}

class S2ClosestCellQueryBaseTest : S2GeometryTestCase() {

fun testMaxDistance() {
  val index = S2CellIndex()
  index.add(makeCellUnion("0/123, 0/22, 0/3"), 1 /*label*/);
  index.build();
  val query = FurthestCellQuery(S2MaxDistanceFactory, index)
  val options = S2ClosestCellQueryBase.Options<S2MaxDistance>(maxResults = 1, distanceFactory = S2MaxDistanceFactory)
  val target = FurthestPointTarget(makeCellId("3/123").toPoint())
  val results = query.findClosestCells(target, options)
  assertEquals(1, results.size)
    assertEquals("0/123", results[0].cellId.toString())
    assertEquals(1, results[0].label);
    assertEquals(4.0, S1ChordAngle(results[0].distance.value).length2)
}

}  // namespace
