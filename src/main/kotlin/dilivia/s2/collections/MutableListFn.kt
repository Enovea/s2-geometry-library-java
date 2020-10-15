package dilivia.s2.collections

import java.util.ArrayList
import kotlin.math.min


fun <T> MutableList<T>.remove(fromIndex: Int, toIndex: Int) {
    require(fromIndex >= 0 && fromIndex <= this.lastIndex)
    require(toIndex >= fromIndex && toIndex <= this.size)
    val nbElementToRemove = toIndex - fromIndex
    repeat(nbElementToRemove) { this.removeAt(fromIndex) }
}

fun <T> MutableList<T>.assign(size: Int, value: T) {
    if (this is ArrayList) {
        this.ensureCapacity(size)
    }
    (0 until min(this.size, size)).forEach { i -> this[i] = value }
    while(this.size < size) this.add(value)
    while (this.size > size) this.removeLast()
}


fun <T : Comparable<T>> MutableList<T>.sortAndRemoveDuplicates() {
    this.sort()
    var idx = this.lastIndex - 1
    while(idx > 1) {
        if (this[idx] == this[idx - 1]) this.removeAt(idx)
        --idx
    }
}

fun <T> MutableList<T>.reverse(startIdx: Int, endIdx: Int) {
    this.subList(startIdx, endIdx).reverse()
}
