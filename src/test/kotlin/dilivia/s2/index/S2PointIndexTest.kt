package dilivia.s2.index

import com.google.common.collect.SortedMultiset
import com.google.common.collect.TreeMultiset
import dilivia.s2.S2CellId
import dilivia.s2.S2GeometryTestCase
import dilivia.s2.S2Point
import dilivia.s2.S2Random
import dilivia.s2.region.S2CellUnion
import mu.KotlinLogging


class S2PointIndexTest : S2GeometryTestCase() {

  private val logger = KotlinLogging.logger {  }

  private val index: S2PointIndex<Int> = S2PointIndex()
  private val contents: SortedMultiset<PointData<Int>> = TreeMultiset.create()


  fun add(point: S2Point, data: Int) {
    index.add(point, data)
    contents.add(PointData(point, data))
  }

  fun remove(point: S2Point, data: Int): Boolean {
    // If there are multiple copies, remove only one.
    contents.remove(PointData(point, data))
    return index.remove(point, data)  // Invalidates "point".
  }

  fun verify() {
    verifyContents()
    verifyIteratorMethods()
  }

  fun verifyContents() {
    logger.info { "Contents: ${contents.size} elements" }
    logger.info { "Index: ${index.numPoints()} elements" }
    val remaining = TreeMultiset.create(contents)
    val iter = index.Iterator()
    var i = 0
    while (!iter.done()) {
      val currentPointData = iter.pointData()
      logger.info { "Iteration $i" }
      logger.info { "Current point : $currentPointData" }
      logger.info { "Point remaining occurences: ${remaining.count(currentPointData)}" }
      val removed = remaining.remove(currentPointData)
      logger.info { "Removed from remaining: $removed" }
      logger.info { "Point remaining occurences after removal: ${remaining.count(currentPointData)}" }
      logger.info { "Remaining: ${remaining.size}" }
      assertTrue(removed)
      iter.next()
      ++i
    }
    assertTrue(remaining.isEmpty())
  }

  fun verifyIteratorMethods() {
    val iter = index.iterator()
    assertFalse(iter.prev())
    iter.finish()
    assertTrue(iter.done())

    // Iterate through all the cells in the index.
    var prev_cellid = S2CellId.none()
    var min_cellid = S2CellId.begin(S2CellId.kMaxLevel)
    iter.begin()
    while (!iter.done()) {
      val cellid = iter.id()
      assertEquals(cellid, S2CellId.fromPoint(iter.point()))
      assertTrue(cellid >= prev_cellid)

      var iter2 = index.iterator()
      if (cellid == prev_cellid) {
        iter2.seek(cellid)
      }

      // Generate a cellunion that covers the range of empty leaf cells between
      // the last cell and this one.  Then make sure that seeking to any of
      // those cells takes us to the immediately following cell.
      if (cellid > prev_cellid) {
        for (skipped in S2CellUnion.fromBeginEnd(min_cellid, cellid)) {
          iter2.seek(skipped)
          assertEquals(cellid, iter2.id())
        }

        // Test Prev(), Next(), and Seek().
        if (prev_cellid.isValid()) {
          iter2.seek(iter.id())
          assertTrue(iter2.prev())
          assertEquals(prev_cellid, iter2.id());
          iter2.next();
          assertEquals(cellid, iter2.id());
          iter2.seek(prev_cellid);
          assertEquals(prev_cellid, iter2.id());
        }
      }
      prev_cellid = cellid;
      min_cellid = cellid.next()
      iter.next()
    }
  }

fun testNoPoints() {
  verify()
}

fun testDuplicatePoints() {
  repeat (10) {
    add(S2Point(1, 0, 0), 123);  // All points have same Data argument.
  }
  assertEquals(10, index.numPoints())
  assertEquals(10, contents.size)
  verify()
  // Now remove half of the points.
  repeat (5) {
    assertTrue(remove(S2Point(1, 0, 0), 123))
  }
  verify()
  assertEquals(5, index.numPoints())
}

fun testRandomPoints() {
  repeat(100) {
    add(S2Random.randomPoint(), S2Random.randomInt(100))
  }
  verify();
  // Now remove some of the points.
  repeat(10) {
    val iter = index.iterator()
    do {
      iter.seek(S2Random.randomCellId(S2CellId.kMaxLevel))
    } while (iter.done())
    remove(iter.point(), iter.data())
    verify();
  }
}

  /*
fun testEmptyData() {
  // Verify that when Data is an empty class, no space is used.
 // assertEquals(sizeof(S2Point), sizeof(S2PointIndex<>::PointData));

  // Verify that points can be added and removed with an empty Data class.
  S2PointIndex<> index;
  index.Add(S2Point(1, 0, 0));
  index.Remove(S2Point(1, 0, 0));
  assertEquals(0, index.num_points());
}*/

}
