package dilivia.s2

fun <T> List<T>.isSorted(): Boolean where T:Comparable<T>{
    if (this.size <= 1) return true

    val iter = this.iterator()
    var current: T
    var previous = iter.next()
    while (iter.hasNext()) {
        current = iter.next();
        if (previous > current) return false
        previous = current;
    }
    return true
}

fun <T: Comparable<T>> List<T>.lowerBound(beginIdx: Int, endIdx: Int, value: T): Int {
    var i = endIdx
    val listIterator = this.listIterator(beginIdx).withIndex()
    while (listIterator.hasNext()) {
        val element = listIterator.next()
        val index = element.index + beginIdx
        if (index >= endIdx) {
            break
        }
        if (element.value >= value) {
            i = index
            break
        }
    }
    return i
}

fun <T: Comparable<T>> List<T>.upperBound(beginIdx: Int, endIdx: Int, value: T): Int {
    var i = endIdx
    val listIterator = this.listIterator(beginIdx).withIndex()
    while (listIterator.hasNext()) {
        val element = listIterator.next()
        val index = element.index + beginIdx
        if (index >= endIdx) {
            break
        }
        if (element.value > value) {
            i = index
            break
        }
    }
    return i
}