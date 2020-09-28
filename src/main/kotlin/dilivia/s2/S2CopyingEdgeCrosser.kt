package dilivia.s2

// TODO(fmeurisse): Is it necessary in kotlin ?

// S2CopyingEdgeCrosser is exactly like S2EdgeCrosser, except that it makes its
// own copy of all arguments so that they do not need to persist between
// calls.  This is less efficient, but makes it possible to use points that
// are generated on demand and cannot conveniently be stored by the client.
class S2CopyingEdgeCrosser {

    private var a: S2Point
    private var b: S2Point
    private var c: S2Point
    // TODO(ericv): It would be more efficient to implement S2CopyingEdgeCrosser
    // directly rather than as a wrapper around S2EdgeCrosser.
    private val crosser: S2EdgeCrosser

    // These methods are all exactly like S2EdgeCrosser, except that the
    // arguments can be temporaries.
    constructor(a: S2Point, b: S2Point) {
        this.a = a
        this.b = b
        this.c = S2Point()
        this.crosser = S2EdgeCrosser(a, b)
    }
    constructor(a: S2Point, b: S2Point, c: S2Point) {
        this.a = a
        this. b = b
        this.c = c
        this.crosser = S2EdgeCrosser(a, b, c)
    }
    fun a(): S2Point { return a }
    fun b(): S2Point { return b }
    fun c(): S2Point { return c }
    fun init(a: S2Point, b: S2Point) {
        this.a = a
        this.b = b
        this.c = S2Point()
        this.crosser.init(a, b)
    }
    fun crossingSign(c: S2Point, d: S2Point): Int {
        if (c != this.c || crosser.c == null) restartAt(c)
        return crossingSign(d)
    }

    fun edgeOrVertexCrossing(c: S2Point, d: S2Point): Boolean {
        if (c != this.c || crosser.c == null) restartAt(c)
        return edgeOrVertexCrossing(d)
    }
    fun restartAt(c: S2Point) {
        this.c = c
        crosser.restartAt(this.c)
    }
    fun crossingSign(d: S2Point): Int {
        val result = crosser.crossingSign(d)
        this.c = d
        crosser.c = this.c
        return result
    }
    fun edgeOrVertexCrossing(d: S2Point): Boolean {
        val result = crosser.edgeOrVertexCrossing(d)
        this.c = d
        crosser.c = this.c
        return result
    }

}