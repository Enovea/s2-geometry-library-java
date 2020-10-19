package dilivia.s2.index.point

import Matrix3x3
import dilivia.s2.Fractal
import dilivia.s2.S2Random
import dilivia.s2.region.S2Cap

// Generates the vertices of a fractal whose convex hull approximately
// matches the given cap.
class FractalPointIndexFactory : PointIndexFactory {

    override fun addPoints(indexCap: S2Cap, numPoints: Int, index: S2PointIndex<Int>) {
        val fractal = Fractal()
        fractal.setLevelForApproxMaxEdges(numPoints)
        fractal.setDimension(1.5)
        val loop = fractal.makeLoop(Matrix3x3.fromCols(S2Random.randomFrameAt(indexCap.center)), indexCap.radius())
        for (i in 0 until loop.numVertices()) {
            index.add(loop.vertex(i), i)
        }
    }
}
