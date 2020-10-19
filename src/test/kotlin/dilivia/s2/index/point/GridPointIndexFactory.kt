package dilivia.s2.index.point

import Matrix3x3
import dilivia.s2.S2Point
import dilivia.s2.S2Random
import dilivia.s2.region.S2Cap
import kotlin.math.ceil
import kotlin.math.sqrt
import kotlin.math.tan

// Generates vertices on a square grid that includes the entire given cap.
class GridPointIndexFactory : PointIndexFactory {

    override fun addPoints(indexCap: S2Cap, numPoints: Int, index: S2PointIndex<Int>) {
        val sqrt_num_points = ceil(sqrt(numPoints.toDouble())).toInt()
        val frame = S2Random.randomFrameAt(indexCap.center)
        val radius = indexCap.radius().radians
        val spacing = 2 * radius / sqrt_num_points;
        for (i in 0 until sqrt_num_points) {
            for (j in 0 until sqrt_num_points) {
                val point = S2Point(tan((i + 0.5) * spacing - radius), tan((j + 0.5) * spacing - radius), 1.0)
                index.add(S2Point.fromFrame(Matrix3x3.fromCols(frame), point.normalize()), i * sqrt_num_points + j)
            }
        }
    }
}
