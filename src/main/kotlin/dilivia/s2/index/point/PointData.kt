package dilivia.s2.index.point

import dilivia.s2.S2Point

// PointData is essentially std::pair with named fields.  It stores an
// S2Point and its associated data, taking advantage of the "empty base
// optimization" to ensure that no extra space is used when Data is empty.
data class PointData<T : Comparable<T>>(val point: S2Point, val data: T) : Comparable<PointData<T>> {

    override fun compareTo(other: PointData<T>): Int {
        val pointComparison = point.compareTo(other.point)
        return if (pointComparison == 0) data.compareTo(other.data) else pointComparison
    }

}
