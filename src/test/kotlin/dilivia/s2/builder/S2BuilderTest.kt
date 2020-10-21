package dilivia.s2.builder

import dilivia.s2.S2CellId
import dilivia.s2.S2GeometryTestCase
import dilivia.s2.S2TextParser
import dilivia.s2.builder.snap.S2CellIdSnapFunction
import dilivia.s2.region.S2Polygon
import dilivia.s2.region.S2PolygonTest

class S2BuilderTest : S2GeometryTestCase() {

    fun test() {
        val kNear0: String = "-1:0, 0:1, 1:0, 0:-1;"
        val kNear1: String = "-1:-1, -1:0, -1:1, 0:1, 1:1, 1:0, 1:-1, 0:-1;"
        val near10: S2Polygon = S2TextParser.makeVerbatimPolygon(kNear0 + kNear1)
        val builder = S2Builder(S2Builder.Options(snapFunction = S2CellIdSnapFunction(S2CellId.kMaxLevel), verbose = true))

        val polygon = S2Polygon()
        polygon.initFromBuilder(near10, builder)
    }
}
