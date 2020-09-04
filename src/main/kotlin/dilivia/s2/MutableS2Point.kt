package dilivia.s2

class MutableS2Point(override val coords: MutableList<Double>): S2Point(coords) {

    constructor(x: Double, y: Double, z: Double): this(mutableListOf(x, y, z))

    operator fun set(idx: Int, value: Double) {
        coords[idx] = value
    }

    fun x(x: Double) {
        set(0, x)
    }

    fun y(y: Double) {
        set(1, x)
    }

    fun z(z: Double) {
        set(2, x)
    }

    override fun newInstance(coords: List<Double>): S2Point = S2Point(if (coords is MutableList) coords else coords.toMutableList())


}