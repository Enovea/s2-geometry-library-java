package dilivia.s2

import junit.framework.TestCase

class S2RandomTest: TestCase() {

    fun testSamplePointInCap() {
        var cap = S2Cap.fromCenterAngle(S2Point(1.0, 0.0, 0.0), S1Angle.degrees(4.0))
        repeat(10000) {
            val point = S2Random.samplePoint(cap)
            assertTrue(cap.contains(point))
        }

        repeat(100) {
            val center = S2Random.randomPoint()
            val angle = S2Random.randomDouble(1.0, 10.0)
            cap = S2Cap.fromCenterAngle(center, S1Angle.degrees(angle))
            repeat(10000) {
                val point = S2Random.samplePoint(cap)
                assertTrue(cap.contains(point))
            }
        }
    }

    fun testSampleCapInRect() {
        repeat(100) {
            val p1 = S2LatLng.fromPoint(S2Random.randomPoint())
            val p2 = S2LatLng.fromPoint(S2Random.randomPoint())
            val rect = S2LatLngRect.fromPointPair(p1, p2)
            repeat(10000) {
                val point = S2Random.samplePoint(rect)
                assertTrue(rect.contains(point))
            }
        }
    }
}