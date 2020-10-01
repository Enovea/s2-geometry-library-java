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
package dilivia.s2.region

import dilivia.s2.S1Angle
import dilivia.s2.S2GeometryTestCase
import dilivia.s2.S2TextParser

class S2ContainsVertexQueryTest : S2GeometryTestCase() {

    fun testUndetermined() {
        val q = S2ContainsVertexQuery(S2TextParser.makePoint("1:2"))
        q.addEdge(S2TextParser.makePoint("3:4"), 1)
        q.addEdge(S2TextParser.makePoint("3:4"), -1)
        assertEquals(0, q.containsSign())
    }

    fun testContainedWithDuplicates() {
        // The S2::Ortho reference direction points approximately due west.
        // Containment is determined by the unmatched edge immediately clockwise.
        val q = S2ContainsVertexQuery(S2TextParser.makePoint("0:0"))
        q.addEdge(S2TextParser.makePoint("3:-3"), -1)
        q.addEdge(S2TextParser.makePoint("1:-5"), 1)
        q.addEdge(S2TextParser.makePoint("2:-4"), 1)
        q.addEdge(S2TextParser.makePoint("1:-5"), -1)
        assertEquals(1, q.containsSign())
    }

    fun testNotContainedWithDuplicates() {
        // The S2::Ortho reference direction points approximately due west.
        // Containment is determined by the unmatched edge immediately clockwise.
        val q = S2ContainsVertexQuery(S2TextParser.makePoint("1:1"))
        q.addEdge(S2TextParser.makePoint("1:-5"), 1)
        q.addEdge(S2TextParser.makePoint("2:-4"), -1)
        q.addEdge(S2TextParser.makePoint("3:-3"), 1)
        q.addEdge(S2TextParser.makePoint("1:-5"), -1)
        assertEquals(-1, q.containsSign())
    }

    fun testMatchesLoopContainment() {
        // Check that the containment function defined is compatible with S2Loop
        // (which at least currently does not use this class).
        val loop = S2Loop.makeRegularLoop(S2TextParser.makePoint("89:-179"), S1Angle.degrees(10), 1000)
        for (i in 1..loop.numVertices()) {
            val q = S2ContainsVertexQuery(loop.vertex(i))
            q.addEdge(loop.vertex(i - 1), -1)
            q.addEdge(loop.vertex(i + 1), 1)
            assertEquals(q.containsSign() > 0, loop.contains(loop.vertex(i)))
        }
    }
}
  
