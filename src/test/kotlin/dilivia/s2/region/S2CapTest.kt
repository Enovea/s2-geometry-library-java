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

import dilivia.s2.*
import dilivia.s2.S2.DBL_EPSILON
import dilivia.s2.S2.M_PI_2
import dilivia.s2.S2.M_PI_4
import dilivia.s2.S2LatLng.Companion.fromDegrees
import dilivia.s2.coords.S2Coords
import kotlin.math.atan
import kotlin.math.sqrt
import kotlin.random.Random

@Strictfp
class S2CapTest : S2GeometryTestCase() {

    fun getLatLngPoint(latDegrees: Double, lngDegrees: Double): S2Point {
        return fromDegrees(latDegrees, lngDegrees).toPoint()
    }

    fun getLatLngPoint(latDegrees: Int, lngDegrees: Int): S2Point = getLatLngPoint(latDegrees.toDouble(), lngDegrees.toDouble())

    fun testBasic() {
        // Test basic properties of empty and full caps.
        val empty = S2Cap.empty
        val full = S2Cap.full
        assertTrue(empty.isValid)
        assertTrue(empty.isEmpty)
        assertTrue(empty.complement.isFull)
        assertTrue(full.isValid)
        assertTrue(full.isFull)
        assertTrue(full.complement.isEmpty)
        assertEquals(2.0, full.height)
        assertEquals(180.0, full.radius().degrees())

        // Test the S1Angle constructor using out-of-range arguments.
        assertTrue(S2Cap.fromCenterAngle(S2Point(1, 0, 0), S1Angle.radians(-20)).isEmpty)
        assertTrue(S2Cap.fromCenterAngle(S2Point(1, 0, 0), S1Angle.radians(5)).isFull)
        assertTrue(S2Cap.fromCenterAngle(S2Point(1, 0, 0), S1Angle.infinity).isFull)

        // Check that the default S2Cap is identical to Empty().
        val defaultEmpty = S2Cap()
        assertTrue(defaultEmpty.isValid)
        assertTrue(defaultEmpty.isEmpty)
        assertEquals(empty.center, defaultEmpty.center)
        assertEquals(empty.height, defaultEmpty.height)

        // Containment and intersection of empty and full caps.
        assertTrue(empty.contains(empty))
        assertTrue(full.contains(empty))
        assertTrue(full.contains(full))
        assertFalse(empty.interiorIntersects(empty))
        assertTrue(full.interiorIntersects(full))
        assertFalse(full.interiorIntersects(empty))

        // Singleton cap containing the x-axis.
        val xaxis = S2Cap.fromPoint(S2Point(1, 0, 0))
        assertTrue(xaxis.contains(S2Point(1, 0, 0)))
        assertFalse(xaxis.contains(S2Point(1.0, 1e-20, 0.0)))
        assertEquals(0.0, xaxis.radius().radians)

        // Singleton cap containing the y-axis.
        val yaxis = S2Cap.fromPoint(S2Point(0, 1, 0))
        assertFalse(yaxis.contains(xaxis.center))
        assertEquals(0.0, xaxis.height)

        // Check that the complement of a singleton cap is the full cap.
        val xcomp = xaxis.complement
        assertTrue(xcomp.isValid)
        assertTrue(xcomp.isFull)
        assertTrue(xcomp.contains(xaxis.center))

        // Check that the complement of the complement is *not* the original.
        assertTrue(xcomp.complement.isValid)
        assertTrue(xcomp.complement.isEmpty)
        assertFalse(xcomp.complement.contains(xaxis.center))

        // Check that very small caps can be represented accurately.
        // Here "kTinyRad" is small enough that unit vectors perturbed by this
        // amount along a tangent do not need to be renormalized.
        val kTinyRad = 1e-10
        val tiny = S2Cap.fromCenterAngle(S2Point(1, 2, 3).normalize(), S1Angle.radians(kTinyRad))
        val tangent = tiny.center.crossProd(S2Point(3, 2, 1)).normalize()
        assertTrue(tiny.contains(tiny.center + 0.99 * kTinyRad * tangent))
        assertFalse(tiny.contains(tiny.center + 1.01 * kTinyRad * tangent))

        // Basic tests on a hemispherical cap.
        val hemi = S2Cap.fromCenterHeight(S2Point(1, 0, 1).normalize(), 1.0)
        assertEquals(-hemi.center, hemi.complement.center)
        assertEquals(1.0, hemi.complement.height)
        assertTrue(hemi.contains(S2Point(1, 0, 0)))
        assertFalse(hemi.complement.contains(S2Point(1, 0, 0)))
        assertTrue(hemi.contains(S2Point(1.0, 0.0, -(1 - kEps)).normalize()))
        assertFalse(hemi.interiorContains(S2Point(1.0, 0.0, -(1 + kEps)).normalize()))

        // A concave cap.  Note that the error bounds for point containment tests
        // increase with the cap angle, so we need to use a larger error bound
        // here.  (It would be painful to do this everywhere, but this at least
        // gives an example of how to compute the maximum error.)
        val center = getLatLngPoint(80.0, 10.0)
        val radius = S1ChordAngle(S1Angle.degrees(150))
        val maxError = (radius.getS2PointConstructorMaxError() + radius.getS1AngleConstructorMaxError() + 3 * DBL_EPSILON)  // getLatLngPoint() error
        val concave = S2Cap(center, radius)
        val concaveMin = S2Cap(center, radius.plusError(-maxError))
        val concaveMax = S2Cap(center, radius.plusError(maxError))
        assertTrue(concaveMax.contains(getLatLngPoint(-70, 10)))
        assertFalse(concaveMin.contains(getLatLngPoint(-70, 10)))
        assertTrue(concaveMax.contains(getLatLngPoint(-50, -170)))
        assertFalse(concaveMin.contains(getLatLngPoint(-50, -170)))

        // Cap containment tests.
        assertFalse(empty.contains(xaxis))
        assertFalse(empty.interiorIntersects(xaxis))
        assertTrue(full.contains(xaxis))
        assertTrue(full.interiorIntersects(xaxis))
        assertFalse(xaxis.contains(full))
        assertFalse(xaxis.interiorIntersects(full))
        assertTrue(xaxis.contains(xaxis))
        assertFalse(xaxis.interiorIntersects(xaxis))
        assertTrue(xaxis.contains(empty))
        assertFalse(xaxis.interiorIntersects(empty))
        assertTrue(hemi.contains(tiny))
        assertTrue(hemi.contains(S2Cap.fromCenterAngle(S2Point(1, 0, 0), S1Angle.radians(M_PI_4 - kEps))))
        assertFalse(hemi.contains(S2Cap.fromCenterAngle(S2Point(1, 0, 0), S1Angle.radians(M_PI_4 + kEps))))
        assertTrue(concave.contains(hemi))
        assertTrue(concave.interiorIntersects(hemi.complement))
        assertFalse(concave.contains(S2Cap.fromCenterHeight(-concave.center, 0.1)))
    }

    fun testAddEmptyCapToNonEmptyCap() {
        val nonEmptyCap = S2Cap.fromCenterAngle(S2Point(1, 0, 0), S1Angle.degrees(10))
        val initialArea = nonEmptyCap.area
        assertEquals(initialArea, nonEmptyCap.addCap(S2Cap.empty).area)
    }

    fun testAddNonEmptyCapToEmptyCap() {
        val empty = S2Cap.empty
        val nonEmptyCap = S2Cap.fromCenterAngle(S2Point(1, 0, 0), S1Angle.degrees(10))
        assertEquals(nonEmptyCap.area, empty.addCap(nonEmptyCap).area)
    }

    fun testGetRectBound() {
        // Empty and full caps.
        assertTrue(S2Cap.empty.rectBound.isEmpty)
        assertTrue(S2Cap.full.rectBound.isFull)

        val kDegreeEps = 1e-13
        // Maximum allowable error for latitudes and longitudes measured in
        // degrees.  (assertEquals isn't sufficient.)

        // Cap that includes the south pole.
        var rect = S2Cap.fromCenterAngle(getLatLngPoint(-45, 57), S1Angle.degrees(50)).rectBound
        assertEquals(rect.latLo().degrees(), -90.0, kDegreeEps)
        assertEquals(rect.latHi().degrees(), 5.0, kDegreeEps)
        assertTrue(rect.lng.isFull)

        // Cap that is tangent to the north pole.
        rect = S2Cap.fromCenterAngle(S2Point(1, 0, 1).normalize(), S1Angle.radians(M_PI_4 + 1e-16)).rectBound
        assertEquals(rect.lat.lo, 0.0, kEps)
        assertEquals(rect.lat.hi, M_PI_2, kEps)
        assertTrue(rect.lng.isFull)

        rect = S2Cap.fromCenterAngle(S2Point(1, 0, 1).normalize(), S1Angle.degrees(45 + 5e-15)).rectBound
        assertEquals(rect.latLo().degrees(), 0.0, kDegreeEps)
        assertEquals(rect.latHi().degrees(), 90.0, kDegreeEps)
        assertTrue(rect.lng.isFull)

        // The eastern hemisphere.
        rect = S2Cap.fromCenterAngle(S2Point(0, 1, 0), S1Angle.radians(M_PI_2 + 2e-16)).rectBound
        assertEquals(rect.latLo().degrees(), -90.0, kDegreeEps)
        assertEquals(rect.latHi().degrees(), 90.0, kDegreeEps)
        assertTrue(rect.lng.isFull)

        // A cap centered on the equator.
        rect = S2Cap.fromCenterAngle(getLatLngPoint(0, 50), S1Angle.degrees(20)).rectBound
        assertEquals(rect.latLo().degrees(), -20.0, kDegreeEps)
        assertEquals(rect.latHi().degrees(), 20.0, kDegreeEps)
        assertEquals(rect.lngLo().degrees(), 30.0, kDegreeEps)
        assertEquals(rect.lngHi().degrees(), 70.0, kDegreeEps)

        // A cap centered on the north pole.
        rect = S2Cap.fromCenterAngle(getLatLngPoint(90, 123), S1Angle.degrees(10)).rectBound
        assertEquals(rect.latLo().degrees(), 80.0, kDegreeEps)
        assertEquals(rect.latHi().degrees(), 90.0, kDegreeEps)
        assertTrue(rect.lng.isFull)
    }

    fun testS2CellMethods() {
        // For each cube face, we construct some cells on
        // that face and some caps whose positions are relative to that face,
        // and then check for the expected intersection/containment results.

        // The distance from the center of a face to one of its vertices.
        val kFaceRadius = atan(sqrt(2.0))

        for (face in 0..5) {
            // The cell consisting of the entire face.
            val rootCell = S2Cell.fromFace(face)

            // A leaf cell at the midpoint of the v=1 edge.
            val edgeCell = S2Cell(S2Coords.faceUVtoXYZ(face, 0.0, 1 - kEps))

            // A leaf cell at the u=1, v=1 corner.
            val cornerCell = S2Cell(S2Coords.faceUVtoXYZ(face, 1 - kEps, 1 - kEps))

            // Quick check for full and empty caps.
            assertTrue(S2Cap.full.contains(rootCell))
            assertFalse(S2Cap.empty.mayIntersect(rootCell))

            // Check intersections with the bounding caps of the leaf cells that are
            // adjacent to 'corner_cell' along the Hilbert curve.  Because this corner
            // is at (u=1,v=1), the curve stays locally within the same cube face.
            var first = cornerCell.id()
            repeat(3) { first = first.previous() }
            var last = cornerCell.id()
            repeat(3) { last = last.next()}
            var id = first
            while (id.lessOrEquals(last)) {
                val cell = S2Cell(id)
                assertEquals(id == cornerCell.id(), cell.capBound.contains(cornerCell))
                assertEquals(id.parent().contains(cornerCell.id()), cell.capBound.mayIntersect(cornerCell))
                id = id.next()
            }

            val antiFace = (face + 3) % 6  // Opposite face.
            for (capFace in 0..5) {
                // A cap that barely contains all of 'capFace'.
                val center = S2Coords.getNorm(capFace)
                val covering = S2Cap.fromCenterAngle(center, S1Angle.radians(kFaceRadius + kEps))
                assertEquals(capFace == face, covering.contains(rootCell))
                assertEquals(capFace != antiFace, covering.mayIntersect(rootCell))
                assertEquals(center.dotProd(edgeCell.getCenter()) > 0.1, covering.contains(edgeCell))
                assertEquals(covering.mayIntersect(edgeCell), covering.contains(edgeCell))
                assertEquals(capFace == face, covering.contains(cornerCell))
                assertEquals(center.dotProd(cornerCell.getCenter()) > 0, covering.mayIntersect(cornerCell))

                // A cap that barely intersects the edges of 'capFace'.
                val bulging = S2Cap.fromCenterAngle(center, S1Angle.radians(M_PI_4 + kEps))
                assertFalse(bulging.contains(rootCell))
                assertEquals(capFace != antiFace, bulging.mayIntersect(rootCell))
                assertEquals("capFace = $capFace == $face != $bulging.contains($edgeCell) == ${bulging.contains(edgeCell)}", capFace == face, bulging.contains(edgeCell))
                assertEquals(center.dotProd(edgeCell.getCenter()) > 0.1, bulging.mayIntersect(edgeCell))
                assertFalse(bulging.contains(cornerCell))
                assertFalse(bulging.mayIntersect(cornerCell))

                /*
                capFace = 0 == 0 != [Center = (1.0, 0.0, 0.0) Radius = dilivia.s2.S1ChordAngle@ce336b0].contains([0, 30, 0, 0/211111111111111111111111111111]) == true
capFace = 1 == 0 != [Center = (0.0, 1.0, 0.0) Radius = dilivia.s2.S1ChordAngle@ce336b0].contains([0, 30, 0, 0/211111111111111111111111111111]) == false
capFace = 2 == 0 != [Center = (0.0, 0.0, 1.0) Radius = dilivia.s2.S1ChordAngle@ce336b0].contains([0, 30, 0, 0/211111111111111111111111111111]) == false
capFace = 3 == 0 != [Center = (-1.0, 0.0, 0.0) Radius = dilivia.s2.S1ChordAngle@ce336b0].contains([0, 30, 0, 0/211111111111111111111111111111]) == false
capFace = 4 == 0 != [Center = (0.0, -1.0, 0.0) Radius = dilivia.s2.S1ChordAngle@ce336b0].contains([0, 30, 0, 0/211111111111111111111111111111]) == false
capFace = 5 == 0 != [Center = (0.0, 0.0, -1.0) Radius = dilivia.s2.S1ChordAngle@ce336b0].contains([0, 30, 0, 0/211111111111111111111111111111]) == false
capFace = 0 == 1 != [Center = (1.0, 0.0, 0.0) Radius = dilivia.s2.S1ChordAngle@ce336b0].contains([1, 30, 2, 1/233333333333333333333333333333]) == false
capFace = 1 == 1 != [Center = (0.0, 1.0, 0.0) Radius = dilivia.s2.S1ChordAngle@ce336b0].contains([1, 30, 2, 1/233333333333333333333333333333]) == true
capFace = 2 == 1 != [Center = (0.0, 0.0, 1.0) Radius = dilivia.s2.S1ChordAngle@ce336b0].contains([1, 30, 2, 1/233333333333333333333333333333]) == false
capFace = 3 == 1 != [Center = (-1.0, 0.0, 0.0) Radius = dilivia.s2.S1ChordAngle@ce336b0].contains([1, 30, 2, 1/233333333333333333333333333333]) == false
capFace = 4 == 1 != [Center = (0.0, -1.0, 0.0) Radius = dilivia.s2.S1ChordAngle@ce336b0].contains([1, 30, 2, 1/233333333333333333333333333333]) == false
capFace = 5 == 1 != [Center = (0.0, 0.0, -1.0) Radius = dilivia.s2.S1ChordAngle@ce336b0].contains([1, 30, 2, 1/233333333333333333333333333333]) == false
capFace = 0 == 2 != [Center = (1.0, 0.0, 0.0) Radius = dilivia.s2.S1ChordAngle@ce336b0].contains([2, 30, 0, 2/211111111111111111111111111111]) == false
capFace = 1 == 2 != [Center = (0.0, 1.0, 0.0) Radius = dilivia.s2.S1ChordAngle@ce336b0].contains([2, 30, 0, 2/211111111111111111111111111111]) == false
capFace = 2 == 2 != [Center = (0.0, 0.0, 1.0) Radius = dilivia.s2.S1ChordAngle@ce336b0].contains([2, 30, 0, 2/211111111111111111111111111111]) == true
capFace = 3 == 2 != [Center = (-1.0, 0.0, 0.0) Radius = dilivia.s2.S1ChordAngle@ce336b0].contains([2, 30, 0, 2/211111111111111111111111111111]) == false
capFace = 4 == 2 != [Center = (0.0, -1.0, 0.0) Radius = dilivia.s2.S1ChordAngle@ce336b0].contains([2, 30, 0, 2/211111111111111111111111111111]) == false
capFace = 5 == 2 != [Center = (0.0, 0.0, -1.0) Radius = dilivia.s2.S1ChordAngle@ce336b0].contains([2, 30, 0, 2/211111111111111111111111111111]) == false
capFace = 0 == 3 != [Center = (1.0, 0.0, 0.0) Radius = dilivia.s2.S1ChordAngle@ce336b0].contains([3, 30, 2, 3/233333333333333333333333333333]) == false
capFace = 1 == 3 != [Center = (0.0, 1.0, 0.0) Radius = dilivia.s2.S1ChordAngle@ce336b0].contains([3, 30, 2, 3/233333333333333333333333333333]) == false
capFace = 2 == 3 != [Center = (0.0, 0.0, 1.0) Radius = dilivia.s2.S1ChordAngle@ce336b0].contains([3, 30, 2, 3/233333333333333333333333333333]) == false
capFace = 3 == 3 != [Center = (-1.0, 0.0, 0.0) Radius = dilivia.s2.S1ChordAngle@ce336b0].contains([3, 30, 2, 3/233333333333333333333333333333]) == true
capFace = 4 == 3 != [Center = (0.0, -1.0, 0.0) Radius = dilivia.s2.S1ChordAngle@ce336b0].contains([3, 30, 2, 3/233333333333333333333333333333]) == false
capFace = 5 == 3 != [Center = (0.0, 0.0, -1.0) Radius = dilivia.s2.S1ChordAngle@ce336b0].contains([3, 30, 2, 3/233333333333333333333333333333]) == false
capFace = 0 == 4 != [Center = (1.0, 0.0, 0.0) Radius = dilivia.s2.S1ChordAngle@ce336b0].contains([4, 30, 0, 4/211111111111111111111111111111]) == false
capFace = 1 == 4 != [Center = (0.0, 1.0, 0.0) Radius = dilivia.s2.S1ChordAngle@ce336b0].contains([4, 30, 0, 4/211111111111111111111111111111]) == false
capFace = 2 == 4 != [Center = (0.0, 0.0, 1.0) Radius = dilivia.s2.S1ChordAngle@ce336b0].contains([4, 30, 0, 4/211111111111111111111111111111]) == false
capFace = 3 == 4 != [Center = (-1.0, 0.0, 0.0) Radius = dilivia.s2.S1ChordAngle@ce336b0].contains([4, 30, 0, 4/211111111111111111111111111111]) == false
capFace = 4 == 4 != [Center = (0.0, -1.0, 0.0) Radius = dilivia.s2.S1ChordAngle@ce336b0].contains([4, 30, 0, 4/211111111111111111111111111111]) == true
capFace = 5 == 4 != [Center = (0.0, 0.0, -1.0) Radius = dilivia.s2.S1ChordAngle@ce336b0].contains([4, 30, 0, 4/211111111111111111111111111111]) == false
capFace = 0 == 5 != [Center = (1.0, 0.0, 0.0) Radius = dilivia.s2.S1ChordAngle@ce336b0].contains([5, 30, 2, 5/233333333333333333333333333333]) == false
capFace = 1 == 5 != [Center = (0.0, 1.0, 0.0) Radius = dilivia.s2.S1ChordAngle@ce336b0].contains([5, 30, 2, 5/233333333333333333333333333333]) == false
capFace = 2 == 5 != [Center = (0.0, 0.0, 1.0) Radius = dilivia.s2.S1ChordAngle@ce336b0].contains([5, 30, 2, 5/233333333333333333333333333333]) == false
capFace = 3 == 5 != [Center = (-1.0, 0.0, 0.0) Radius = dilivia.s2.S1ChordAngle@ce336b0].contains([5, 30, 2, 5/233333333333333333333333333333]) == false
capFace = 4 == 5 != [Center = (0.0, -1.0, 0.0) Radius = dilivia.s2.S1ChordAngle@ce336b0].contains([5, 30, 2, 5/233333333333333333333333333333]) == false
capFace = 5 == 5 != [Center = (0.0, 0.0, -1.0) Radius = dilivia.s2.S1ChordAngle@ce336b0].contains([5, 30, 2, 5/233333333333333333333333333333]) == true
                 */

                // A singleton cap.
                val singleton = S2Cap.fromCenterAngle (center, S1Angle.zero)
                assertEquals(capFace == face, singleton.mayIntersect(rootCell))
                assertFalse(singleton.mayIntersect(edgeCell))
                assertFalse(singleton.mayIntersect(cornerCell))
            }
        }
    }

    fun testGetCellUnionBoundLevel1Radius() {
        // Check that a cap whose radius is approximately the width of a level 1
        // S2Cell can be covered by only 3 faces.
        val cap = S2Cap.fromCenterAngle(S2Point(1, 1, 1).normalize(), S1Angle.radians(S2CellMetrics.kMinWidth.getValue(1)))
        val covering = mutableListOf<S2CellId>()
        cap.getCellUnionBound(covering)
        assertEquals(3, covering.size)
    }

    fun testExpanded() {
        assertTrue(S2Cap.empty.expanded(S1Angle.radians(2)).isEmpty)
        assertTrue(S2Cap.full.expanded(S1Angle.radians(2)).isFull)
        val cap50 = S2Cap.fromCenterAngle(S2Point(1, 0, 0), S1Angle.degrees(50))
        val cap51 = S2Cap.fromCenterAngle(S2Point(1, 0, 0), S1Angle.degrees(51))
        assertTrue(cap50.expanded(S1Angle.radians(0)).approxEquals(cap50))
        assertTrue(cap50.expanded(S1Angle.degrees(1)).approxEquals(cap51))
        assertFalse(cap50.expanded(S1Angle.degrees(129.99)).isFull)
        assertTrue(cap50.expanded(S1Angle.degrees(130.01)).isFull)
    }

    fun testGetCentroid() {
        // Empty and full caps.
        assertEquals(S2Point(), S2Cap.empty.centroid)
        assertTrue(S2Cap.full.centroid.norm() <= 1e-15)

        // Random caps.
        for (i in 0 until 100) {
            val center = S2Random.randomPoint()
            val height = Random.nextDouble(0.0, 2.0)
            val cap = S2Cap.fromCenterHeight(center, height)
            val centroid = cap.centroid
            val expected = center * (1.0 - height / 2.0) * cap.area
            assertTrue((expected - centroid).norm() <= 1e-15)
        }
    }

    fun testUnion() {
        // Two caps which have the same center but one has a larger radius.
        val a = S2Cap.fromCenterAngle(getLatLngPoint(50.0, 10.0), S1Angle.degrees(0.2))
        val b = S2Cap.fromCenterAngle(getLatLngPoint(50.0, 10.0), S1Angle.degrees(0.3))
        assertTrue(b.contains(a))
        assertEquals(b, a.union(b))

        // Two caps where one is the full cap.
        assertTrue(a.union(S2Cap.full).isFull)

        // Two caps where one is the empty cap.
        assertEquals(a, a.union(S2Cap.empty))

        // Two caps which have different centers, one entirely encompasses the other.
        val c = S2Cap.fromCenterAngle(getLatLngPoint(51.0, 11.0), S1Angle.degrees(1.5))
        assertTrue(c.contains(a))
        assertEquals(a.union(c).center, c.center)
        assertEquals(a.union(c).radius(), c.radius())

        // Two entirely disjoint caps.
        val d = S2Cap.fromCenterAngle(getLatLngPoint(51.0, 11.0), S1Angle.degrees(0.1))
        assertFalse(d.contains(a))
        assertFalse(d.intersects(a))
        assertTrue(a.union(d).approxEquals(d.union(a)))
        assertEquals(50.4588, S2LatLng.fromPoint(a.union(d).center).lat().degrees(), 0.001)
        assertEquals(10.4525, S2LatLng.fromPoint(a.union(d).center).lng().degrees(), 0.001)
        assertEquals(0.7425, a.union(d).radius().degrees(), 0.001)

        // Two partially overlapping caps.
        val e = S2Cap.fromCenterAngle(getLatLngPoint(50.3, 10.3), S1Angle.degrees(0.2))
        assertFalse(e.contains(a))
        assertTrue(e.intersects(a))
        assertTrue(a.union(e).approxEquals(e.union(a)))
        assertEquals(50.1500, S2LatLng.fromPoint(a.union(e).center).lat().degrees(), 0.001)
        assertEquals(10.1495, S2LatLng.fromPoint(a.union(e).center).lng().degrees(), 0.001)
        assertEquals(0.3781, a.union(e).radius().degrees(), 0.001)

        // Two very large caps, whose radius sums to in excess of 180 degrees, and
        // whose centers are not antipodal.
        val f = S2Cap.fromCenterAngle(S2Point(0, 0, 1).normalize(), S1Angle.degrees(150))
        val g = S2Cap.fromCenterAngle(S2Point(0, 1, 0).normalize(), S1Angle.degrees(150))
        assertTrue(f.union(g).isFull)

        // Two non-overlapping hemisphere caps with antipodal centers.
        val hemi = S2Cap.fromCenterHeight(S2Point(0, 0, 1).normalize(), 1.0)
        assertTrue(hemi.union(hemi.complement).isFull)
    }

    /*
    fun testEncodeDecode() {
        S2Cap cap = S2Cap . fromCenterHeight (S2Point(3, 2, 1).normalize(), 1)
        Encoder encoder
                cap.Encode(& encoder)
        Decoder decoder (encoder.base(), encoder.length())
        S2Cap decoded_cap
                assertTrue(decoded_cap.Decode(& decoder))
        assertEquals(cap, decoded_cap)
    }
    */

    companion object {
        // About 9 times the double-precision roundoff relative error.
        const val kEps = 1e-15
    }
}
