package dilivia.s2


// The XYZ and face,si,ti coordinates of an S2Point and, if this point is equal
// to the center of an S2Cell, the level of this cell (-1 otherwise).
data class S2XYZFaceSiTi(val xyz: S2Point, val face: Int, val si: UInt, val ti: UInt, val cellLevel: Int)