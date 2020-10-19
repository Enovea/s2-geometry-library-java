package dilivia.s2.index.shape

import dilivia.s2.Fractal
import dilivia.s2.S2Random
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
