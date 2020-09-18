/**
 * This project is a kotlin port of the Google s2 geometry library (Copyright 2005 Google Inc. All Rights Reserved.):
 *                                 https://github.com/google/s2geometry.git
 *
 * Copyright © 2020 Dilivia (contact@dilivia.com)
 *
 * Licensed under the Apache License, Version 2.0 (the "License")
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
package dilivia.s2

import dilivia.s2.Assertions.assertGE
import dilivia.s2.Assertions.assertLT
import dilivia.s2.S2Coords.kPosToOrientation
import dilivia.s2.math.R2Point
import java.util.logging.Logger
import kotlin.math.IEEErem
import kotlin.math.abs
import kotlin.math.min
import kotlin.random.Random

/**
 */
@ExperimentalUnsignedTypes
@Strictfp
class S2CellIdTest : S2GeometryTestCase() {

    fun testDefaultConstructor() {
        val id = S2CellId()
        assertEquals(0UL, id.id)
        assertFalse(id.isValid())
    }

    fun testS2CellIdHash() {
        assertEquals(getCellId(0, 90).hashCode(), getCellId(0, 90).hashCode())
    }

    fun testFaceDefinitions() {
        assertEquals(0, getCellId(0, 0).face());
        assertEquals(1, getCellId(0, 90).face());
        assertEquals(2, getCellId(90, 0).face());
        assertEquals(3, getCellId(0, 180).face());
        assertEquals(4, getCellId(0, -90).face());
        assertEquals(5, getCellId(-90, 0).face());
    }

    fun testFromFace() {
        for (face in 0..5) {
            assertEquals(S2CellId.fromFacePosLevel(face, 0UL, 0), S2CellId.fromFace(face))
        }
    }

    fun testParentChildRelationships() {
        val id = S2CellId.fromFacePosLevel(3, 0x12345678.toULong(), S2CellId.kMaxLevel - 4);
        assertTrue(id.isValid());
        assertEquals(3, id.face());
        assertEquals(0x12345700.toULong(), id.pos());
        assertEquals(S2CellId.kMaxLevel - 4, id.level());
        assertFalse(id.isLeaf());

        assertEquals(0x12345610.toULong(), id.childBegin(id.level() + 2).pos());
        assertEquals(0x12345640.toULong(), id.childBegin().pos());
        assertEquals(0x12345400.toULong(), id.parent().pos());
        assertEquals(0x12345000.toULong(), id.parent(id.level() - 2).pos());

        // Check ordering of children relative to parents.
        assertLessThan(id.childBegin(), id)
        assertGreaterThan(id.childEnd(), id);
        assertEquals(id.childEnd(), id.childBegin().next().next().next().next());
        assertEquals(id.rangeMin(), id.childBegin(S2CellId.kMaxLevel));
        assertEquals(id.rangeMax().next(), id.childEnd(S2CellId.kMaxLevel));

        // Check that cells are represented by the position of their center
        // along the Hilbert curve.
        assertEquals(2UL * id.id, id.rangeMin().id + id.rangeMax().id);
    }

    fun testSentinelRangeMinMax() {
        assertEquals(S2CellId.sentinel(), S2CellId.sentinel().rangeMin());
        assertEquals(S2CellId.sentinel(), S2CellId.sentinel().rangeMax());
    }

    fun testCenterSiTi() {
        val id = S2CellId.fromFacePosLevel(3, 0x12345678.toULong(), S2CellId.kMaxLevel);
        // Check that the (si, ti) coordinates of the center end in a
        // 1 followed by (30 - level) 0s.

        // Leaf level, 30.

        var faceSiTi = id.getCenterSiTi()
        assertEquals(1U shl 0, faceSiTi.si and 1U);
        assertEquals(1U shl 0, faceSiTi.ti and 1U);

        // Level 29.
        faceSiTi = id.parent(S2CellId.kMaxLevel - 1).getCenterSiTi();
        assertEquals(1U shl 1, faceSiTi.si and 3U);
        assertEquals(1U shl 1, faceSiTi.ti and 3U);

        // Level 28.
        faceSiTi = id.parent(S2CellId.kMaxLevel - 2).getCenterSiTi();
        assertEquals(1U shl 2, faceSiTi.si and 7U);
        assertEquals(1U shl 2, faceSiTi.ti and 7U);

        // Level 20.
        faceSiTi = id.parent(S2CellId.kMaxLevel - 10).getCenterSiTi();
        assertEquals(1U shl 10, faceSiTi.si and ((1U shl 11) - 1U));
        assertEquals(1U shl 10, faceSiTi.ti and ((1U shl 11) - 1U));

        // Level 10.
        faceSiTi = id.parent(S2CellId.kMaxLevel - 20).getCenterSiTi();
        assertEquals(1U shl 20, faceSiTi.si and ((1U shl 21) - 1U));
        assertEquals(1U shl 20, faceSiTi.ti and ((1U shl 21) - 1U));

        // Level 0.
        faceSiTi = id.parent(0).getCenterSiTi();
        assertEquals(1U shl 30, faceSiTi.si and ((1U shl 31) - 1U));
        assertEquals(1U shl 30, faceSiTi.ti and ((1U shl 31) - 1U));
    }

    fun testWrapping() {
        // Check wrapping from beginning of Hilbert curve to end and vice versa.
        assertEquals(S2CellId.end(0).previous(), S2CellId.begin(0).prevWrap());

        assertEquals(
                S2CellId.fromFacePosLevel(5, 0UL.inv() shr S2CellId.kFaceBits, S2CellId.kMaxLevel),
                S2CellId.begin(S2CellId.kMaxLevel).prevWrap()
        )
        assertEquals(
                S2CellId.fromFacePosLevel(5, 0UL.inv() shr S2CellId.kFaceBits, S2CellId.kMaxLevel),
                S2CellId.begin(S2CellId.kMaxLevel).advanceWrap(-1)
        )

        assertEquals(S2CellId.begin(4), S2CellId.end(4).previous().nextWrap());
        assertEquals(S2CellId.begin(4), S2CellId.end(4).advance(-1).advanceWrap(1));

        assertEquals(
                S2CellId.fromFacePosLevel(0, 0UL, S2CellId.kMaxLevel),
                S2CellId.end(S2CellId.kMaxLevel).previous().nextWrap()
        );
        assertEquals(
                S2CellId.fromFacePosLevel(0, 0UL, S2CellId.kMaxLevel),
                S2CellId.end(S2CellId.kMaxLevel).advance(-1).advanceWrap(1)
        );
    }

    fun testAdvance() {
        val id = S2CellId.fromFacePosLevel(3, 0x12345678.toULong(), S2CellId.kMaxLevel - 4)
        // Check basic properties of advance().
        assertEquals(S2CellId.end(0), S2CellId.begin(0).advance(7));
        assertEquals(S2CellId.end(0), S2CellId.begin(0).advance(12));
        assertEquals(S2CellId.begin(0), S2CellId.end(0).advance(-7));
        assertEquals(S2CellId.begin(0), S2CellId.end(0).advance(-12000000));
        val num_level_5_cells = 6 shl (2 * 5)
        assertEquals(S2CellId.end(5).advance(500L - num_level_5_cells), S2CellId.begin(5).advance(500));
        assertEquals(id.next().childBegin(S2CellId.kMaxLevel), id.childBegin(S2CellId.kMaxLevel).advance(256));
        assertEquals(
                S2CellId.fromFacePosLevel(5, 0UL, S2CellId.kMaxLevel),
                S2CellId.fromFacePosLevel(1, 0UL, S2CellId.kMaxLevel).advance(4L shl (2 * S2CellId.kMaxLevel))
        );

        // Check basic properties of advance_wrap().
        assertEquals(S2CellId.fromFace(1), S2CellId.begin(0).advanceWrap(7));
        assertEquals(S2CellId.begin(0), S2CellId.begin(0).advanceWrap(12));
        assertEquals(S2CellId.fromFace(4), S2CellId.fromFace(5).advanceWrap(-7));
        assertEquals(S2CellId.begin(0), S2CellId.begin(0).advanceWrap(-12000000));
        assertEquals(S2CellId.begin(5).advanceWrap(6644), S2CellId.begin(5).advanceWrap(-11788));
        assertEquals(id.next().childBegin(S2CellId.kMaxLevel), id.childBegin(S2CellId.kMaxLevel).advanceWrap(256));
        assertEquals(
                S2CellId.fromFacePosLevel(1, 0UL, S2CellId.kMaxLevel),
                S2CellId.fromFacePosLevel(5, 0UL, S2CellId.kMaxLevel).advanceWrap(2L shl (2 * S2CellId.kMaxLevel))
        )
    }

    fun testDistanceFromBegin() {
        assertEquals(6.toULong(), S2CellId.end(0).distanceFromBegin());
        assertEquals(6UL * (1L shl (2 * S2CellId.kMaxLevel)).toULong(), S2CellId.end(S2CellId.kMaxLevel).distanceFromBegin());

        assertEquals(0UL, S2CellId.begin(0).distanceFromBegin());
        assertEquals(0UL, S2CellId.begin(S2CellId.kMaxLevel).distanceFromBegin());

        val id = S2CellId.fromFacePosLevel(3, 0x12345678.toULong(), S2CellId.kMaxLevel - 4);
        assertEquals(id, S2CellId.begin(id.level()).advance(id.distanceFromBegin().toLong()))
    }

    fun testMaximumTile() {
        // This method is tested more thoroughly in s2cell_union_test.cc.
        for (iter in 0 until 1000) {
            val id = getRandomCellId(10)

            // Check that "limit" is returned for tiles at or beyond "limit".
            assertEquals(id, id.maximumTile(id));
            assertEquals(id, id.child(0).maximumTile(id));
            assertEquals(id, id.child(1).maximumTile(id));
            assertEquals(id, id.next().maximumTile(id));
            assertEquals(id.child(0), id.maximumTile(id.child(0)));

            // Check that the tile size is increased when possible.
            assertEquals(id, id.child(0).maximumTile(id.next()));
            assertEquals(id, id.child(0).maximumTile(id.next().child(0)));
            assertEquals(id, id.child(0).maximumTile(id.next().child(1).child(0)));
            assertEquals(id, id.child(0).child(0).maximumTile(id.next()));
            assertEquals(id, id.child(0).child(0).child(0).maximumTile(id.next()));

            // Check that the tile size is decreased when necessary.
            assertEquals(id.child(0), id.maximumTile(id.child(0).next()));
            assertEquals(id.child(0), id.maximumTile(id.child(0).next().child(0)));
            assertEquals(id.child(0), id.maximumTile(id.child(0).next().child(1)));
            assertEquals(id.child(0).child(0), id.maximumTile(id.child(0).child(0).next()));
            assertEquals(id.child(0).child(0).child(0), id.maximumTile(id.child(0).child(0).child(0).next()));

            // Check that the tile size is otherwise unchanged.
            assertEquals(id, id.maximumTile(id.next()));
            assertEquals(id, id.maximumTile(id.next().child(0)));
            assertEquals(id, id.maximumTile(id.next().child(1).child(0)));
        }
    }

    fun testGetCommonAncestorLevel() {
        // Two identical cell ids.
        assertEquals(0, S2CellId.fromFace(0).getCommonAncestorLevel(S2CellId.fromFace(0)));
        assertEquals(30, S2CellId.fromFace(0).childBegin(30).getCommonAncestorLevel(S2CellId.fromFace(0).childBegin(30)));

        // One cell id is a descendant of the other.
        assertEquals(0, S2CellId.fromFace(0).childBegin(30).getCommonAncestorLevel(S2CellId.fromFace(0)));
        assertEquals(0, S2CellId.fromFace(5).getCommonAncestorLevel(S2CellId.fromFace(5).childEnd(30).previous()));

        // Two cells that have no common ancestor.
        assertEquals(-1, S2CellId.fromFace(0).getCommonAncestorLevel(S2CellId.fromFace(5)));
        assertEquals(-1, S2CellId.fromFace(2).childBegin(30).getCommonAncestorLevel(S2CellId.fromFace(3).childEnd(20)));

        // Two cells that have a common ancestor distinct from both of them.
        assertEquals(8, S2CellId.fromFace(5).childBegin(9).next().childBegin(15).getCommonAncestorLevel(
                S2CellId.fromFace(5).childBegin(9).childBegin(20)));
        assertEquals(1, S2CellId.fromFace(0).childBegin(2).childBegin(30).getCommonAncestorLevel(
                S2CellId.fromFace(0).childBegin(2).next().childBegin(5)));
    }

    fun testInverses() {
        // Check the conversion of random leaf cells to S2LatLngs and back.
        for (i in 0 until 200000) {
            val id = getRandomCellId(S2CellId.kMaxLevel)
            assertTrue(id.isLeaf())
            assertEquals(S2CellId.kMaxLevel, id.level())
            val center = id.toLatLng()
            assertEquals(id.id, S2CellId.fromLatLng(center).id);
        }
    }

    fun testTokens() {

        // Test random cell ids at all levels.
        for (i in 0 until 10000) {
            val id = randomCellId
            val token = id.toToken();
            assertLessOrEquals(token.length, 16);
            assertEquals(id, S2CellId.fromToken(token))
        }
        // Check that invalid cell ids can be encoded, and round-trip is
        // the identity operation.
        var token = S2CellId.none().toToken();
        assertEquals(S2CellId.none(), S2CellId.fromToken(token));

        // Sentinel is invalid.
        token = S2CellId.sentinel().toToken();
        assertEquals(S2CellId.fromToken(token), S2CellId.sentinel())

        // Check an invalid face.
        token = S2CellId.fromFace(7).toToken();
        assertEquals(S2CellId.fromToken(token), S2CellId.fromFace(7));

        // Check that supplying tokens with non-alphanumeric characters
        // returns S2CellId.none().
        assertEquals(S2CellId.none(), S2CellId.fromToken("876b e99"));
        assertEquals(S2CellId.none(), S2CellId.fromToken("876bee99\n"));
        assertEquals(S2CellId.none(), S2CellId.fromToken("876[ee99"));
        assertEquals(S2CellId.none(), S2CellId.fromToken(" 876bee99"));
    }

    fun testContainment() {
        // Test contains() and intersects().
        val parent_map = mutableMapOf<S2CellId, S2CellId>()
        val cells = mutableListOf<S2CellId>()
        for (face in 0..5) {
            expandCell(S2CellId.fromFace(face), cells, parent_map)
        }
        for (end_id in cells) {
            for (begin_id in cells) {
                var contained = true;
                var id = begin_id
                while (id != end_id) {
                    if (parent_map[id] == null) {
                        contained = false;
                        break;
                    }
                    id = parent_map[id]!!
                }
                assertEquals(contained, end_id.contains(begin_id));
                assertEquals(contained, begin_id in end_id.rangeMin()..end_id.rangeMax())
                assertEquals(end_id.intersects(begin_id), end_id.contains(begin_id) || begin_id.contains(end_id))
            }
        }
    }

    fun testContinuity() {
        // Make sure that sequentially increasing cell ids form a continuous
        // path over the surface of the sphere, i.e. there are no
        // discontinuous jumps from one region to another.

        val maxDist = S2CellMetrics.kMaxEdge.getValue(kMaxWalkLevel)
        val end = S2CellId.end(kMaxWalkLevel)
        var id = S2CellId.begin(kMaxWalkLevel)
        while (id != end) {
            assertLessOrEquals(id.toPointRaw().angle(id.nextWrap().toPointRaw()), maxDist);
            assertEquals(id.nextWrap(), id.advanceWrap(1));
            assertEquals(id, id.nextWrap().advanceWrap(-1));

            // Check that the ToPointRaw() returns the center of each cell
            // in (s,t) coordinates.
            val (_, u, v) = S2Coords.xyzToFaceUV(id.toPointRaw())
            val kCellSize = 1.0 / (1 shl kMaxWalkLevel)
            assertEquals(S2Coords.uvToSt(u).IEEErem(0.5 * kCellSize), 0.0, 1e-15);
            assertEquals(S2Coords.uvToSt(v).IEEErem(0.5 * kCellSize), 0.0, 1e-15);

            id = id.next()
        }
    }

    fun testCoverage() {
        // Make sure that random points on the sphere can be represented to the
        // expected level of accuracy, which in the worst case is sqrt(2/3) times
        // the maximum arc length between the points on the sphere associated with
        // adjacent values of "i" or "j".  (It is sqrt(2/3) rather than 1/2 because
        // the cells at the corners of each face are stretched -- they have 60 and
        // 120 degree angles.)
        val max_dist = 0.5 * S2CellMetrics.kMaxDiag.getValue(S2CellId.kMaxLevel);
        for (i in 0 until 1000000) {
            val p = randomPoint()
            val q = S2CellId.fromPoint(p).toPointRaw()
            assertLessOrEquals(p.angle(q), max_dist)
        }
    }

    fun testNeighbors() {
        // Check the edge neighbors of face 1.
        val out_faces = intArrayOf(5, 3, 2, 0)
        val face_nbrs = S2CellId.fromFace(1).getEdgeNeighbors()
        for (i in 0..3) {
            assertTrue(face_nbrs[i].isFace());
            assertEquals(out_faces[i], face_nbrs[i].face());
        }

        // Check the edge neighbors of the corner cells at all levels.  This case is
        // trickier because it requires projecting onto adjacent faces.
        val kMaxIJ = S2CellId.kMaxSize - 1
        for (level in 1..S2CellId.kMaxLevel) {
            val id = S2CellId.fromFaceIJ(1, 0, 0).parent(level)
            val nbrs = id.getEdgeNeighbors()
            // These neighbors were determined manually using the face and axis
            // relationships defined in s2coords.cc.
            val size_ij = S2CellId.getSizeIJ(level);
            assertEquals(S2CellId.fromFaceIJ(5, kMaxIJ, kMaxIJ).parent(level), nbrs[0]);
            assertEquals(S2CellId.fromFaceIJ(1, size_ij, 0).parent(level), nbrs[1]);
            assertEquals(S2CellId.fromFaceIJ(1, 0, size_ij).parent(level), nbrs[2]);
            assertEquals(S2CellId.fromFaceIJ(0, kMaxIJ, 0).parent(level), nbrs[3]);
        }

        // Check the vertex neighbors of the center of face 2 at level 5.
        val nbrs = mutableListOf<S2CellId>()
        S2CellId.fromPoint(S2Point(0, 0, 1)).appendVertexNeighbors(5, nbrs)
        nbrs.sort()
        for (i in 0..3) {
            assertEquals(S2CellId.fromFaceIJ(2, (1 shl 29) - if (i < 2) 1 else 0, (1 shl 29) - if (i == 0 || i == 3) 1 else 0).parent(5), nbrs[i]);
        }
        nbrs.clear();

        // Check the vertex neighbors of the corner of faces 0, 4, and 5.
        val id = S2CellId.fromFacePosLevel(0, 0UL, S2CellId.kMaxLevel);
        id.appendVertexNeighbors(0, nbrs)
        nbrs.sort()
        assertEquals(3, nbrs.size);
        assertEquals(S2CellId.fromFace(0), nbrs[0]);
        assertEquals(S2CellId.fromFace(4), nbrs[1]);
        assertEquals(S2CellId.fromFace(5), nbrs[2]);

        // Check that AppendAllNeighbors produces results that are consistent
        // with AppendVertexNeighbors for a bunch of random cells.
        for (i in 0 until 1000) {
            var id = randomCellId
            if (id.isLeaf()) id = id.parent()

            // TestAllNeighbors computes approximately 2**(2*(diff+1)) cell ids,
            // so it's not reasonable to use large values of "diff".
            val max_diff = min(5, S2CellId.kMaxLevel - id.level() - 1);
            val level = id.level() + random(max_diff + 1)
            testAllNeighbors(id, level)
        }
    }

    fun testExpandedByDistanceUV(id: S2CellId, distance: S1Angle) {
        val bound = id.getBoundUV()
        val expanded = S2CellId.expandedByDistanceUV(bound, distance);
        for (iter in 0 until 100) {
            // Choose a point on the boundary of the rectangle.
            val face = random(6);
            val center_uv = sampleBoundary(bound)
            val center = S2Coords.faceUVtoXYZ(face, center_uv.x(), center_uv.y()).normalize()

            // Now sample a point from a disc of radius (2 * distance).
            val p = randomPoint(S2Cap.fromCenterAngle(center, 2 * distance.abs()))

            // Find the closest point on the boundary to the sampled point.
            val uv = S2Coords.faceXYZtoUV(face, p) ?: continue

            val closest_uv = projectToBoundary(uv, bound)
            val closest = S2Coords.faceUVtoXYZ(face, closest_uv[0], closest_uv[1]).normalize()
            val actual_dist = S1Angle(p, closest);

            if (distance >= S1Angle.zero) {
                // "expanded" should contain all points in the original bound, and also
                // all points within "distance" of the boundary.
                if (bound.contains(uv) || actual_dist < distance) {
                    assertTrue(expanded.contains(uv));
                }
            } else {
                // "expanded" should not contain any points within "distance" of the
                // original boundary.
                if (actual_dist < -distance) {
                    assertFalse(expanded.contains(uv));
                }
            }
        }
    }

    fun testExpandedByDistanceUV() {
        val max_dist_degrees = 10.0
        for (iter in 0 until 100) {
            val id = randomCellId
            val cellSize = S2Cell(id).boundUV().size
            val min = S1Angle.radians(min(cellSize[0], cellSize[1])).degrees()
            val dist_degrees = Random.Default.nextDouble(-min, max_dist_degrees);
            testExpandedByDistanceUV(id, S1Angle.degrees(dist_degrees))
        }
    }

    fun testToString() {
        assertEquals("3/", S2CellId.fromFace(3).toString());
        assertEquals("4/000000000000000000000000000000", S2CellId.fromFace(4).rangeMin().toString());
        assertEquals("Invalid: 0000000000000000", S2CellId.none().toString());
    }

    fun testFromDebugString() {
        assertEquals(S2CellId.fromFace(3), S2CellId.fromDebugString("3/"));
        assertEquals(S2CellId.fromFace(0).child(2).child(1), S2CellId.fromDebugString("0/21"));
        assertEquals(S2CellId.fromFace(4).rangeMin(), S2CellId.fromDebugString("4/000000000000000000000000000000"));
        assertEquals(S2CellId.none(), S2CellId.fromDebugString("4/0000000000000000000000000000000"));
        assertEquals(S2CellId.none(), S2CellId.fromDebugString(""));
        assertEquals(S2CellId.none(), S2CellId.fromDebugString("7/"));
        assertEquals(S2CellId.none(), S2CellId.fromDebugString(" /"));
        assertEquals(S2CellId.none(), S2CellId.fromDebugString("3:0"));
        assertEquals(S2CellId.none(), S2CellId.fromDebugString("3/ 12"));
        assertEquals(S2CellId.none(), S2CellId.fromDebugString("3/1241"));
    }

    fun testOutputOperator() {
        val cell = S2CellId (0xbb04000000000000UL)
        assertEquals("5/31200", cell.toString())
    }

    companion object {

        private val logger = Logger.getLogger(S2CellIdTest::class.qualifiedName!!)

        const val kMaxExpandLevel = 3

        const val kMaxWalkLevel = 8

        fun getCellId(lat_degrees: Int, lng_degrees: Int): S2CellId = getCellId(lat_degrees.toDouble(), lng_degrees.toDouble())

        fun getCellId(lat_degrees: Double, lng_degrees: Double): S2CellId {
            val id = S2CellId.fromLatLng(S2LatLng.fromDegrees(lat_degrees, lng_degrees))
            logger.info { "CellId from $lat_degrees - $lng_degrees : $id = ${id.id} " }
            return id;
        }

        fun expandCell(parent: S2CellId, cells: MutableList<S2CellId>, parent_map: MutableMap<S2CellId, S2CellId>) {
            cells.add(parent)
            if (parent.level() == kMaxExpandLevel) return
            val (face, _, _, orientation) = parent.toFaceIJOrientation(true)
            assertEquals(parent.face(), face);

            var child = parent.childBegin()
            var pos = 0
            while (child != parent.childEnd()) {
                parent_map[child] = parent
                // Do some basic checks on the children.
                assertEquals(child, parent.child(pos));
                assertEquals(pos, child.childPosition());
                // Test child_position(level) on all the child's ancestors.
                var ancestor = child
                while (ancestor.level() >= 1) {
                    assertEquals(child.childPosition(ancestor.level()), ancestor.childPosition());
                    ancestor = parent_map[ancestor]!!
                }
                assertEquals(pos, child.childPosition(child.level()));
                assertEquals(parent.level() + 1, child.level());
                assertFalse(child.isLeaf());

                val (childFace, _, _, child_orientation) = child.toFaceIJOrientation(true)
                assertEquals(face, childFace);
                assertEquals(orientation!! xor kPosToOrientation[pos], child_orientation);
                expandCell(child, cells, parent_map);

                child = child.next()
                ++pos
            }
        }

        fun testAllNeighbors(id: S2CellId, level: Int) {
            assertGE(level, id.level());
            assertLT(level, S2CellId.kMaxLevel);

            // We compute AppendAllNeighbors, and then add in all the children of "id"
            // at the given level.  We then compare this against the result of finding
            // all the vertex neighbors of all the vertices of children of "id" at the
            // given level.  These should give the same result.
            val all = mutableListOf<S2CellId>()
            val expected = mutableListOf<S2CellId>()
            id.appendAllNeighbors(level, all)
            val end = id.childEnd(level + 1)
            var c = id.childBegin(level + 1)
            while (c != end) {
                all.add(c.parent())
                c.appendVertexNeighbors(level, expected)
                c = c.next()
            }
            // Sort the results and eliminate duplicates.
            assertEquals(expected.toSortedSet(), all.toSortedSet());
        }

        // Returns a random point on the boundary of the given rectangle.
        fun sampleBoundary(rect: R2Rect): R2Point {
            val uv = doubleArrayOf(0.0, 0.0)
            val d = Random.Default.nextInt(2)
            uv[d] = Random.Default.nextDouble(rect[d][0], rect[d][1])
            uv[1 - d] = if (Random.Default.nextBoolean()) rect[1 - d][0] else rect[1 - d][1]
            return R2Point(uv);
        }

        // Returns the closest point to "uv" on the boundary of "rect".
        fun projectToBoundary(uv: R2Point, rect: R2Rect): R2Point {
            val du0 = abs(uv[0] - rect[0][0]);
            val du1 = abs(uv[0] - rect[0][1]);
            val dv0 = abs(uv[1] - rect[1][0]);
            val dv1 = abs(uv[1] - rect[1][1]);
            val dmin = min(min(du0, du1), min(dv0, dv1));
            if (du0 == dmin) return R2Point(rect[0][0], rect[1].project(uv[1]));
            if (du1 == dmin) return R2Point(rect[0][1], rect[1].project(uv[1]));
            if (dv0 == dmin) return R2Point(rect[0].project(uv[0]), rect[1][0]);
            assertEquals("Bug in ProjectToBoundary", dmin, dv1)
            return R2Point(rect[0].project(uv[0]), rect[1][1]);
        }

    }

}