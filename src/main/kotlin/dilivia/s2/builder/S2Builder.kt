package dilivia.s2.builder

import com.google.common.geometry.S2
import com.google.common.geometry.S2.DBL_EPSILON
import dilivia.s2.Assertions.assertGE
import dilivia.s2.Assertions.assertLE
import dilivia.s2.S1Angle
import dilivia.s2.S1ChordAngle
import dilivia.s2.S2CellId
import dilivia.s2.S2EdgeCrossings
import dilivia.s2.S2EdgeDistances
import dilivia.s2.S2Error
import dilivia.s2.S2Point
import dilivia.s2.collections.assign
import dilivia.s2.builder.layers.Layer
import dilivia.s2.index.MutableS2ShapeIndex
import dilivia.s2.index.S2PointIndex
import dilivia.s2.region.S2Loop
import dilivia.s2.region.S2Polygon
import dilivia.s2.region.S2Polyline
import dilivia.s2.shape.S2Shape
import dilivia.s2.sin
import kotlin.math.acos


//////////////////////  Input Types  /////////////////////////
// All types associated with the S2Builder inputs are prefixed with "Input".

// Identifies an input vertex.
typealias InputVertexId = Int

// Defines an input edge.
typealias InputEdge = Pair<InputVertexId, InputVertexId>

// Identifies an input edge.
typealias InputEdgeId = Int

// Identifies the set of input edge ids that were snapped to a given edge.
typealias InputEdgeIdSetId = Int

// Sort key for prioritizing input vertices.  (Note that keys are *not*
// compared using std::less; see SortInputVertices for details.)
typealias InputVertexKey = Pair<S2CellId, InputVertexId>

//////////////////////  Output Types  /////////////////////////
// These types define the output vertices and edges.

// Identifies a snapped vertex ("snap site").  If there is only one layer,
// than SiteId is the same as Graph::VertexId, but if there are many layers
// then each Graph may contain only a subset of the sites.  Also see
// GraphOptions::allow_vertex_filtering().
typealias SiteId = Int

// Defines an output edge.
typealias Edge = Pair<SiteId, SiteId>
fun min(edge1: Edge, edge2: Edge) = when {
    edge1.first < edge2.first -> edge1
    edge1.first > edge2.first -> edge2
    edge1.second > edge2.second -> edge2
    else -> edge1
}

// Identifies an output edge.
typealias EdgeId = Int

// Identifies an output edge in a particular layer.
typealias LayerEdgeId = Pair<Int, EdgeId>

// Every edge can have a set of non-negative integer labels attached to it.
// When used with an appropriate layer type, you can then retrieve the
// labels associated with each output edge.  This can be useful when merging
// or combining data from several sources.  (Note that in many cases it is
// easier to use separate output layers rather than labels.)
//
// Labels are 32-bit non-negative integers.  To support other label types,
// you can use ValueLexicon to store the set of unique labels seen so far:
//
//   ValueLexicon<MyLabel> my_label_lexicon;
//   builder.set_label(my_label_lexicon.Add(label));
//
// The current set of labels is represented as a stack.  This makes it easy
// to add and remove labels hierarchically (e.g., polygon 5, loop 2).  Use
// set_label() and clear_labels() if you need at most one label per edge.
//
typealias Label = Int

// Each input edge has "label set id" (an int32) representing the set of
// labels attached to that edge.  This vector is populated only if at least
// one label is used.
typealias LabelSetId = Int;

// For output layers that represent polygons, there is an ambiguity inherent
// in spherical geometry that does not exist in planar geometry.  Namely, if
// a polygon has no edges, does it represent the empty polygon (containing
// no points) or the full polygon (containing all points)?  This ambiguity
// also occurs for polygons that consist only of degeneracies, e.g. a
// degenerate loop with only two edges could be either a degenerate shell in
// the empty polygon or a degenerate hole in the full polygon.
//
// To resolve this ambiguity, an IsFullPolygonPredicate may be specified for
// each output layer (see AddIsFullPolygonPredicate below).  If the output
// after snapping consists only of degenerate edges and/or sibling pairs
// (including the case where there are no edges at all), then the layer
// implementation calls the given predicate to determine whether the polygon
// is empty or full except for those degeneracies.  The predicate is given
// an S2Builder::Graph containing the output edges, but note that in general
// the predicate must also have knowledge of the input geometry in order to
// determine the correct result.
//
// This predicate is only needed by layers that are assembled into polygons.
// It is not used by other layer types.
interface IsFullPolygonPredicate {
    fun test(g: Graph, error: S2Error): Boolean
}

// Indicates whether the input edges are undirected.  Typically this is
// specified for each output layer (e.g., s2builderutil::S2PolygonLayer).
//
// Directed edges are preferred, since otherwise the output is ambiguous.
// For example, output polygons may be the *inverse* of the intended result
// (e.g., a polygon intended to represent the world's oceans may instead
// represent the world's land masses).  Directed edges are also somewhat
// more efficient.
//
// However even with undirected edges, most S2Builder layer types try to
// preserve the input edge direction whenever possible.  Generally, edges
// are reversed only when it would yield a simpler output.  For example,
// S2PolygonLayer assumes that polygons created from undirected edges should
// cover at most half of the sphere.  Similarly, S2PolylineVectorLayer
// assembles edges into as few polylines as possible, even if this means
// reversing some of the "undirected" input edges.
//
// For shapes with interiors, directed edges should be oriented so that the
// interior is to the left of all edges.  This means that for a polygon with
// holes, the outer loops ("shells") should be directed counter-clockwise
// while the inner loops ("holes") should be directed clockwise.  Note that
// S2Builder::AddPolygon() follows this convention automatically.
enum class EdgeType { DIRECTED, UNDIRECTED }

// An S2Shape used to represent the entire collection of S2Builder input edges.
// Vertices are specified as indices into a vertex vector to save space.
// Requires that "edges" is constant for the lifetime of this object.
class VertexIdEdgeVectorShape(private val edges: List<Pair<Int, Int>>, private val vertices: List<S2Point>) : S2Shape() {

    fun vertex0(e: Int): S2Point = vertex(edges[e].first)
    fun vertex1(e: Int): S2Point = vertex(edges[e].second)

    // S2Shape interface:

    override val dimension: Int = 1
    override val numEdges: Int = edges.size
    override fun edge(edgeId: Int): dilivia.s2.shape.Edge = dilivia.s2.shape.Edge(vertices[edges[edgeId].first], vertices[edges[edgeId].second])

    override fun getReferencePoint(): ReferencePoint = ReferencePoint(contained = false)

    override val numChains: Int = edges.size
    override fun chain(chain_id: Int): Chain = Chain(chain_id, 1)
    override fun chainEdge(chainId: Int, offset: Int): dilivia.s2.shape.Edge = edge(chainId)

    override fun chainPosition(edgeId: Int): ChainPosition = ChainPosition(edgeId, 0)

    private fun vertex(i: Int): S2Point = vertices[i]

}


//
// This class is a replacement for S2PolygonBuilder.  Once all clients have
// been updated to use this class, S2PolygonBuilder will be removed.


// S2Builder is a tool for assembling polygonal geometry from edges.  Here are
// some of the things it is designed for:
//
// 1. Building polygons, polylines, and polygon meshes from unsorted
//    collections of edges.
//
// 2. Snapping geometry to discrete representations (such as S2CellId centers
//    or E7 lat/lng coordinates) while preserving the input topology and with
//    guaranteed error bounds.
//
// 3. Simplifying geometry (e.g. for indexing, display, or storage).
//
// 4. Importing geometry from other formats, including repairing geometry
//    that has errors.
//
// 5. As a tool for implementing more complex operations such as polygon
//    intersections and unions.
//
// The implementation is based on the framework of "snap rounding".  Unlike
// most snap rounding implementations, S2Builder defines edges as geodesics on
// the sphere (straight lines) and uses the topology of the sphere (i.e.,
// there are no "seams" at the poles or 180th meridian).  The algorithm is
// designed to be 100% robust for arbitrary input geometry.  It offers the
// following properties:
//
//   - Guaranteed bounds on how far input vertices and edges can move during
//     the snapping process (i.e., at most the given "snap_radius").
//
//   - Guaranteed minimum separation between edges and vertices other than
//     their endpoints (similar to the goals of Iterated Snap Rounding).  In
//     other words, edges that do not intersect in the output are guaranteed
//     to have a minimum separation between them.
//
//   - Idempotency (similar to the goals of Stable Snap Rounding), i.e. if the
//     input already meets the output criteria then it will not be modified.
//
//   - Preservation of the input topology (up to the creation of
//     degeneracies).  This means that there exists a continuous deformation
//     from the input to the output such that no vertex crosses an edge.  In
//     other words, self-intersections won't be created, loops won't change
//     orientation, etc.
//
//   - The ability to snap to arbitrary discrete point sets (such as S2CellId
//     centers, E7 lat/lng points on the sphere, or simply a subset of the
//     input vertices), rather than being limited to an integer grid.
//
// Here are some of its other features:
//
//  - It can handle both directed and undirected edges.  Undirected edges can
//    be useful for importing data from other formats, e.g. where loops have
//    unspecified orientations.
//
//  - It can eliminate self-intersections by finding all edge pairs that cross
//    and adding a new vertex at each intersection point.
//
//  - It can simplify polygons to within a specified tolerance.  For example,
//    if two vertices are close enough they will be merged, and if an edge
//    passes nearby a vertex then it will be rerouted through that vertex.
//    Optionally, it can also detect nearly straight chains of short edges and
//    replace them with a single long edge, while maintaining the same
//    accuracy, separation, and topology guarantees ("simplify_edge_chains").
//
//  - It supports many different output types through the concept of "layers"
//    (polylines, polygons, polygon meshes, etc).  You can build multiple
//    layers at once in order to ensure that snapping does not create
//    intersections between different objects (for example, you can simplify a
//    set of contour lines without the risk of having them cross each other).
//
//  - It supports edge labels, which allow you to attach arbitrary information
//    to edges and have it preserved during the snapping process.  (This can
//    also be achieved using layers, at a coarser level of granularity.)
//
// Caveats:
//
//  - Because S2Builder only works with edges, it cannot distinguish between
//    the empty and full polygons.  If your application can generate both the
//    empty and full polygons, you must implement logic outside of this class.
//
// Example showing how to snap a polygon to E7 coordinates:
//
//  using s2builderutil::IntLatLngSnapFunction;
//  S2Builder builder(S2Builder::Options(IntLatLngSnapFunction(7)));
//  S2Polygon output;
//  builder.StartLayer(absl::make_unique<s2builderutil::S2PolygonLayer>(&output));
//  builder.AddPolygon(input);
//  S2Error error;
//  if (!builder.Build(&error)) {
//    S2_LOG(ERROR) << error;
//    ...
//  }
class S2Builder(val options: Options = Options()) {

        //////////// Parameters

        // The maximum distance (inclusive) that a vertex can move when snapped,
        // equal to S1ChordAngle(options_.snap_function().snap_radius()).
        private var siteSnapRadiusCa: S1ChordAngle

        // The maximum distance (inclusive) that an edge can move when snapping to a
        // snap site.  It can be slightly larger than the site snap radius when
        // edges are being split at crossings.
        private var edgeSnapRadiusCa: S1ChordAngle

        private var maxEdgeDeviation: S1Angle
        private var edgeSiteQueryRadiusCa: S1ChordAngle
        private var minEdgeLengthToSplitCa: S1ChordAngle

        private var minSiteSeparation: S1Angle
        private var minSiteSeparationCa: S1ChordAngle
        private var minEdgeSiteSeparationCa: S1ChordAngle
        private var minEdgeSiteSeparationCaLimit: S1ChordAngle

        private var maxAdjacentSiteSeparationCa: S1ChordAngle

        // The squared sine of the edge snap radius.  This is equivalent to the snap
        // radius (squared) for distances measured through the interior of the
        // sphere to the plane containing an edge.  This value is used only when
        // interpolating new points along edges (see GetSeparationSite).
        private var edgeSnapRadiusSin2: Double

        // A copy of the argument to Build().
        private var error: S2Error = S2Error(code = S2Error.OK, text = "")

        // True if snapping was requested.  This is true if either snap_radius() is
        // positive, or split_crossing_edges() is true (which implicitly requests
        // snapping to ensure that both crossing edges are snapped to the
        // intersection point).
        private var snappingRequested: Boolean

        // Initially false, and set to true when it is discovered that at least one
        // input vertex or edge does not meet the output guarantees (e.g., that
        // vertices are separated by at least snap_function.min_vertex_separation).
        private var snappingNeeded: Boolean

        //////////// Input Data /////////////

        // A flag indicating whether label_set_ has been modified since the last
        // time label_set_id_ was computed.
        private var labelSetModified: Boolean

        private var inputVertices: MutableList<S2Point> = mutableListOf()
        private var inputEdges: MutableList<InputEdge> = mutableListOf()

        private var layers: MutableList<Layer> = mutableListOf()
        private var layerOptions: MutableList<GraphOptions> = mutableListOf()
        private var layerBegins: MutableList<InputEdgeId> = mutableListOf()
        private var layerIsFullPolygonPredicates: MutableList<IsFullPolygonPredicate> = mutableListOf()

        private var labelSetIds: MutableList<LabelSetId> = mutableListOf()
        private var label_set_lexicon: IdSetLexicon = IdSetLexicon()

        // The current set of labels (represented as a stack).
        private var labelSet: MutableList<Label> = mutableListOf()

        // The LabelSetId corresponding to the current label set, computed on demand
        // (by adding it to label_set_lexicon()).
        private var labelSetId: LabelSetId

        ////////////// Data for Snapping and Simplifying //////////////

        // The number of sites specified using ForceVertex().  These sites are
        // always at the beginning of the sites_ vector.
        private var num_forced_sites: SiteId = -1

        // The set of snapped vertex locations ("sites").
        private var sites: MutableList<S2Point> = mutableListOf()

        // A map from each input edge to the set of sites "nearby" that edge,
        // defined as the set of sites that are candidates for snapping and/or
        // avoidance.  Note that compact_array will inline up to two sites, which
        // usually takes care of the vast majority of edges.  Sites are kept sorted
        // by increasing distance from the origin of the input edge.
        //
        // Once snapping is finished, this field is discarded unless edge chain
        // simplification was requested, in which case instead the sites are
        // filtered by removing the ones that each edge was snapped to, leaving only
        // the "sites to avoid" (needed for simplification).
        private var edgeSites: MutableList<MutableList<SiteId>> = mutableListOf()

        // Initializes an S2Builder with the given options.
        init {
            val snapFunction = options.snapFunction
            val snapRadius = snapFunction.snapRadius
            assertLE(snapRadius, SnapFunction.kMaxSnapRadius)

            // Convert the snap radius to an S1ChordAngle.  This is the "true snap
            // radius" used when evaluating exact predicates (s2predicates.h).
            siteSnapRadiusCa = S1ChordAngle(snapRadius)

            // When split_crossing_edges() is true, we need to use a larger snap radius
            // for edges than for vertices to ensure that both edges are snapped to the
            // edge intersection location.  This is because the computed intersection
            // point is not exact; it may be up to kIntersectionError away from its true
            // position.  The computed intersection point might then be snapped to some
            // other vertex up to snap_radius away.  So to ensure that both edges are
            // snapped to a common vertex, we need to increase the snap radius for edges
            // to at least the sum of these two values (calculated conservatively).
            var edgeSnapRadius = snapRadius
            if (!options.simplifyEdgeChains) {
                edgeSnapRadiusCa = siteSnapRadiusCa
            } else {
                edgeSnapRadius += S2EdgeCrossings.kIntersectionError
                edgeSnapRadiusCa = roundUp(edgeSnapRadius)
            }
            snappingRequested = (edgeSnapRadius > S1Angle.zero)

            // Compute the maximum distance that a vertex can be separated from an
            // edge while still affecting how that edge is snapped.
            maxEdgeDeviation = snapFunction.maxEdgeDeviation()
            edgeSiteQueryRadiusCa = S1ChordAngle(maxEdgeDeviation + snapFunction.minEdgeVertexSeparation())

            // Compute the maximum edge length such that even if both endpoints move by
            // the maximum distance allowed (i.e., snap_radius), the center of the edge
            // will still move by less than max_edge_deviation().  This saves us a lot
            // of work since then we don't need to check the actual deviation.
            minEdgeLengthToSplitCa = S1ChordAngle.radians(2 * acos(sin(snapRadius) / sin(maxEdgeDeviation)))

            // If the condition below is violated, then AddExtraSites() needs to be
            // modified to check that snapped edges pass on the same side of each "site
            // to avoid" as the input edge.  Currently it doesn't need to do this
            // because the condition below guarantees that if the snapped edge passes on
            // the wrong side of the site then it is also too close, which will cause a
            // separation site to be added.
            //
            // Currently max_edge_deviation() is at most 1.1 * snap_radius(), whereas
            // min_edge_vertex_separation() is at least 0.219 * snap_radius() (based on
            // S2CellIdSnapFunction, which is currently the worst case).
            assertLE(snapFunction.maxEdgeDeviation(),snapFunction.snapRadius + snapFunction.minEdgeVertexSeparation())

            // To implement idempotency, we check whether the input geometry could
            // possibly be the output of a previous S2Builder invocation.  This involves
            // testing whether any site/site or edge/site pairs are too close together.
            // This is done using exact predicates, which require converting the minimum
            // separation values to an S1ChordAngle.
            minSiteSeparation = snapFunction.minVertexSeparation()
            minSiteSeparationCa = S1ChordAngle(minSiteSeparation)
            minEdgeSiteSeparationCa = S1ChordAngle(snapFunction.minEdgeVertexSeparation())

            // This is an upper bound on the distance computed by S2ClosestPointQuery
            // where the true distance might be less than min_edge_site_separation_ca_.
            minEdgeSiteSeparationCaLimit = addPointToEdgeError(minEdgeSiteSeparationCa)

            // Compute the maximum possible distance between two sites whose Voronoi
            // regions touch.  (The maximum radius of each Voronoi region is
            // edge_snap_radius_.)  Then increase this bound to account for errors.
            maxAdjacentSiteSeparationCa = addPointToPointError(roundUp(edgeSnapRadius * 2.0))

            // Finally, we also precompute sin^2(edge_snap_radius), which is simply the
            // squared distance between a vertex and an edge measured perpendicular to
            // the plane containing the edge, and increase this value by the maximum
            // error in the calculation to compare this distance against the bound.
            val d = sin(edgeSnapRadius)
            edgeSnapRadiusSin2 = d * d;
            edgeSnapRadiusSin2 += ((9.5 * d + 2.5 + 2 * S2.M_SQRT3) * d + 9 * DBL_EPSILON) * DBL_EPSILON

            // Initialize the current label set.
            labelSetId = IdSetLexicon.emptySetId()
            labelSetModified = false

            // If snapping was requested, we try to determine whether the input geometry
            // already meets the output requirements.  This is necessary for
            // idempotency, and can also save work.  If we discover any reason that the
            // input geometry needs to be modified, snapping_needed_ is set to true.
            snappingNeeded = false;
        }

    data class Options(

            // Sets the desired snap function.  The snap function is copied
            // internally, so you can safely pass a temporary object.
            //
            // Note that if your input data includes vertices that were created using
            // S2::GetIntersection(), then you should use a "snap_radius" of
            // at least S2::kIntersectionSnapRadius, e.g. by calling
            //
            //  options.set_snap_function(s2builderutil::IdentitySnapFunction(
            //      S2::kIntersectionSnapRadius));
            //
            // DEFAULT: s2builderutil::IdentitySnapFunction(S1Angle::Zero())
            // [This does no snapping and preserves all input vertices exactly.]
            val snapFunction: SnapFunction = IdentitySnapFunction(S1Angle.zero),

            // If true, then detect all pairs of crossing edges and eliminate them by
            // adding a new vertex at their intersection point.
            //
            // When this option is true, the effective snap_radius() for edges is
            // increased by S2::kIntersectionError to take into account the
            // additional error when computing intersection points.  In other words,
            // edges may move by up to snap_radius() + S2::kIntersectionError.
            //
            // Undirected edges should always be used when the output is a polygon,
            // since splitting a directed loop at a self-intersection converts it into
            // two loops that don't define a consistent interior according to the
            // "interior is on the left" rule.  (On the other hand, it is fine to use
            // directed edges when defining a polygon *mesh* because in that case the
            // input consists of sibling edge pairs.)
            //
            // Self-intersections can also arise when importing data from a 2D
            // projection.  You can minimize this problem by subdividing the input
            // edges so that the S2 edges (which are geodesics) stay close to the
            // original projected edges (which are curves on the sphere).  This can
            // be done using s2builderutil::EdgeSplitter(), for example.
            //
            // DEFAULT: false
            val splitCrossingEdges: Boolean = false,

            // If true, then simplify the output geometry by replacing nearly straight
            // chains of short edges with a single long edge.
            //
            // The combined effect of snapping and simplifying will not change the
            // input by more than the guaranteed tolerances (see the list documented
            // with the SnapFunction class).  For example, simplified edges are
            // guaranteed to pass within snap_radius() of the *original* positions of
            // all vertices that were removed from that edge.  This is a much tighter
            // guarantee than can be achieved by snapping and simplifying separately.
            //
            // However, note that this option does not guarantee idempotency.  In
            // other words, simplifying geometry that has already been simplified once
            // may simplify it further.  (This is unavoidable, since tolerances are
            // measured with respect to the original geometry, which is no longer
            // available when the geometry is simplified a second time.)
            //
            // When the output consists of multiple layers, simplification is
            // guaranteed to be consistent: for example, edge chains are simplified in
            // the same way across layers, and simplification preserves topological
            // relationships between layers (e.g., no crossing edges will be created).
            // Note that edge chains in different layers do not need to be identical
            // (or even have the same number of vertices, etc) in order to be
            // simplified together.  All that is required is that they are close
            // enough together so that the same simplified edge can meet all of their
            // individual snapping guarantees.
            //
            // Note that edge chains are approximated as parametric curves rather than
            // point sets.  This means that if an edge chain backtracks on itself (for
            // example, ABCDEFEDCDEFGH) then such backtracking will be preserved to
            // within snap_radius() (for example, if the preceding point were all in a
            // straight line then the edge chain would be simplified to ACFCFH, noting
            // that C and F have degree > 2 and therefore can't be simplified away).
            //
            // Simplified edges are assigned all labels associated with the edges of
            // the simplified chain.
            //
            // For this option to have any effect, a SnapFunction with a non-zero
            // snap_radius() must be specified.  Also note that vertices specified
            // using ForceVertex are never simplified away.
            //
            // DEFAULT: false
            val simplifyEdgeChains: Boolean = false,

            // If true, then snapping occurs only when the input geometry does not
            // already meet the S2Builder output guarantees (see the SnapFunction
            // class description for details).  This means that if all input vertices
            // are at snapped locations, all vertex pairs are separated by at least
            // min_vertex_separation(), and all edge-vertex pairs are separated by at
            // least min_edge_vertex_separation(), then no snapping is done.
            //
            // If false, then all vertex pairs and edge-vertex pairs closer than
            // "snap_radius" will be considered for snapping.  This can be useful, for
            // example, if you know that your geometry contains errors and you want to
            // make sure that features closer together than "snap_radius" are merged.
            //
            // This option is automatically turned off by simplify_edge_chains(),
            // since simplifying edge chains is never guaranteed to be idempotent.
            //
            // DEFAULT: true
            val idempotent: Boolean = true
    )


  // Starts a new output layer.  This method must be called before adding any
  // edges to the S2Builder.  You may call this method multiple times to build
  // multiple geometric objects that are snapped to the same set of sites.
  //
  // For example, if you have a set of contour lines, then you could put each
  // contour line in a separate layer.  This keeps the contour lines separate
  // from each other, while also ensuring that no crossing edges are created
  // when they are snapped and/or simplified.  (This is not true if the
  // contour lines are snapped or simplified independently.)
  //
  // Similarly, if you have a set of polygons that share common boundaries
  // (e.g., countries), you can snap and/or simplify them at the same time by
  // putting them in different layers, while ensuring that their boundaries
  // remain consistent (i.e., no crossing edges or T-vertices are introduced).
  //
  // Ownership of the layer is transferred to the S2Builder.  Example usage:
  //
  // S2Polyline line1, line2;
  // builder.StartLayer(make_unique<s2builderutil::S2PolylineLayer>(&line1)));
  // ... Add edges using builder.AddEdge(), etc ...
  // builder.StartLayer(make_unique<s2builderutil::S2PolylineLayer>(&line2)));
  // ... Add edges using builder.AddEdge(), etc ...
  // S2Error error;
  // S2_CHECK(builder.Build(&error)) << error;  // Builds "line1" & "line2"
  fun startLayer(layer: Layer) {
      layerOptions.add(layer.graphOptions())
      layerBegins.add(inputEdges.size)
      layerIsFullPolygonPredicates.add(IsFullPolygon(false))
      layers.add(layer)
  }

  // Adds a degenerate edge (representing a point) to the current layer.
  fun addPoint(v: S2Point) {
      addEdge(v, v)
  }

  // Adds the given edge to the current layer.
  fun addEdge(v0: S2Point, v1: S2Point) {
      check(layers.isNotEmpty()) { "Call StartLayer before adding any edges" }

      if (v0 == v1 && (layerOptions.last().degenerate_edges == DegenerateEdges.DISCARD)) {
          return
      }
      val j0 = addVertex(v0)
      val j1 = addVertex(v1)
      inputEdges.add(InputEdge(j0, j1))

      // If there are any labels, then attach them to this input edge.
      if (labelSetModified) {
          if (labelSetIds.isEmpty()) {
              // Populate the missing entries with empty label sets.
              labelSetIds.assign(inputEdges.size - 1, labelSetId)
          }
          labelSetId = label_set_lexicon.add(labelSet)
          labelSetIds.add(labelSetId)
          labelSetModified = false
      } else if (labelSetIds.isNotEmpty()) {
          labelSetIds.add(labelSetId)
      }
  }

  // Adds the edges in the given polyline.  (Note that if the polyline
  // consists of 0 or 1 vertices, this method does nothing.)
  fun addPolyline(polyline: S2Polyline) {
      val n = polyline.numVertices()
      for (i in 1 until n) {
          addEdge(polyline.vertex(i - 1), polyline.vertex(i))
      }
  }

  // Adds the edges in the given loop.  If the sign() of the loop is negative
  // (i.e. this loop represents a hole within a polygon), the edge directions
  // are automatically reversed to ensure that the polygon interior is always
  // to the left of every edge.
  fun addLoop(loop: S2Loop) {
      // Ignore loops that do not have a boundary.
      if (loop.isEmptyOrFull()) return

      // For loops that represent holes, we add the edge from vertex n-1 to vertex
      // n-2 first.  This is because these edges will be assembled into a
      // clockwise loop, which will eventually be normalized in S2Polygon by
      // calling S2Loop::Invert().  S2Loop::Invert() reverses the order of the
      // vertices, so to end up with the original vertex order (0, 1, ..., n-1) we
      // need to build a clockwise loop with vertex order (n-1, n-2, ..., 0).
      // This is done by adding the edge (n-1, n-2) first, and then ensuring that
      // Build() assembles loops starting from edges in the order they were added.
      val n = loop.numVertices()
      for (i in 0 until n) {
          addEdge(loop.orientedVertex(i), loop.orientedVertex(i + 1))
      }
  }

  // Adds the loops in the given polygon.  Loops representing holes have their
  // edge directions automatically reversed as described for AddLoop().  Note
  // that this method does not distinguish between the empty and full polygons,
  // i.e. adding a full polygon has the same effect as adding an empty one.
  fun addPolygon(polygon: S2Polygon) {
      for (i in 0 until polygon.numLoops()) {
          addLoop(polygon.loop(i))
      }
  }

  // Adds the edges of the given shape to the current layer.
  fun addShape(shape: S2Shape) {
      var e = 0
      val n = shape.numEdges
      while (e < n) {
          val edge = shape.edge(e)
          addEdge(edge.v0, edge.v1)
          ++e
      }
  }

  // For layers that are assembled into polygons, this method specifies a
  // predicate that is called when the output consists entirely of degenerate
  // edges and/or sibling pairs.  The predicate is given an S2Builder::Graph
  // containing the output edges (if any) and is responsible for deciding
  // whether this graph represents the empty polygon (possibly with degenerate
  // shells) or the full polygon (possibly with degenerate holes).  Note that
  // this cannot be determined from the output edges alone; it also requires
  // knowledge of the input geometry.  (Also see IsFullPolygonPredicate above.)
  //
  // This method should be called at most once per layer; additional calls
  // simply overwrite the previous value for the current layer.
  //
  // The default predicate simply returns false (i.e., degenerate polygons are
  // assumed to be empty).  Arguably it would better to return an error in
  // this case, but the fact is that relatively few clients need to be able to
  // construct full polygons, and it is unreasonable to expect all such
  // clients to supply an appropriate predicate.
  //
  // The reason for having a predicate rather than a boolean value is that the
  // predicate is responsible for determining whether the output polygon is
  // empty or full.  In general the input geometry is not degenerate, but
  // rather collapses into a degenerate configuration due to snapping and/or
  // simplification.
  //
  // TODO(ericv): Provide standard predicates to handle common cases,
  // e.g. valid input geometry that becomes degenerate due to snapping.
  fun addIsFullPolygonPredicate(predicate: IsFullPolygonPredicate) {
      layerIsFullPolygonPredicates[layerIsFullPolygonPredicates.lastIndex] = predicate
  }
      // A predicate that returns an error indicating that no polygon predicate
      // has been specified.}
    class IsFullPolygonUnspecified: IsFullPolygonPredicate {
        override fun test(g: Graph, error: S2Error): Boolean {
            error.code = S2Error.BUILDER_IS_FULL_PREDICATE_NOT_SPECIFIED
            return false
        }
    }

    // Returns a predicate that returns a constant value (true or false);
    class IsFullPolygon(val isFull: Boolean): IsFullPolygonPredicate {
        override fun test(g: Graph, error: S2Error): Boolean {
            return isFull
        }
    }

  // Forces a vertex to be located at the given position.  This can be used to
  // prevent certain input vertices from moving.  However if you are trying to
  // preserve part of the input boundary, be aware that this option does not
  // prevent edges from being split by new vertices.
  //
  // Forced vertices are never snapped; if this is desired then you need to
  // call options().snap_function().SnapPoint() explicitly.  Forced vertices
  // are also never simplified away (if simplify_edge_chains() is used).
  //
  // Caveat: Since this method can place vertices arbitrarily close together,
  // S2Builder makes no minimum separation guaranteees with forced vertices.
  fun forceVertex(vertex: S2Point) {
      sites.add(vertex)
  }

  // Clear the stack of labels.
  fun clearLabels() {
      labelSet.clear();
      labelSetModified = true
  }

  // Add a label to the stack.
  // REQUIRES: label >= 0.
  fun pushLabel(label: Label) {
      assertGE(label, 0)
      labelSet.add(label)
      labelSetModified = true
  }

  // Remove a label from the stack.
  fun popLabel() {
      labelSet.removeLast()
      labelSetModified = true
  }

  // Convenience function that clears the stack and adds a single label.
  // REQUIRES: label >= 0.
  fun setLabel(label: Label) {
      assertGE(label, 0)
      labelSet.clear()
      labelSet.add(label)
      labelSetModified = true
  }

  // Performs the requested edge splitting, snapping, simplification, etc, and
  // then assembles the resulting edges into the requested output layers.
  //
  // Returns true if all edges were assembled; otherwise sets "error"
  // appropriately.  Depending on the error, some or all output layers may
  // have been created.  Automatically resets the S2Builder state so that it
  // can be reused.
  //
  // REQUIRES: error != nullptr.
  fun build(): S2Error {
      // S2_CHECK rather than S2_DCHECK because this is friendlier than crashing on the
      // "error->ok()" call below.  It would be easy to allow (error == nullptr)
      // by declaring a local "tmp_error", but it seems better to make clients
      // think about error handling.
      error.code = S2Error.OK
      error.text = ""

      // Mark the end of the last layer.
      layerBegins.add(inputEdges.size)

      // See the algorithm overview at the top of this file.
      if (snappingRequested && !options.idempotent) {
          snappingNeeded = true
      }
      chooseSites()
      buildLayers()
      reset();
      return error.copy()
  }

  // Clears all input data and resets the builder state.  Any options
  // specified are preserved.
  fun reset() {
      inputVertices.clear()
      inputEdges.clear()
      layers.clear()
      layerOptions.clear()
      layerBegins.clear()
      layerIsFullPolygonPredicates.clear()
      labelSetIds.clear()
      label_set_lexicon.clear()
      labelSet.clear()
      labelSetModified = false
      sites.clear()
      edgeSites.clear()
      snappingNeeded = false
  }

  private fun addVertex(v: S2Point): InputVertexId {
      // Remove duplicate vertices that follow the pattern AB, BC, CD.  If we want
      // to do anything more sophisticated, either use a ValueLexicon, or sort the
      // vertices once they have all been added, remove duplicates, and update the
      // edges.
      if (inputVertices.isEmpty() || v != inputVertices.last()) {
          inputVertices.add(v);
      }
      return inputVertices.size - 1
  }

  private fun chooseSites(): Unit {
      if (inputVertices.isEmpty()) return

      val inputEdgeIndex = MutableS2ShapeIndex()
      inputEdgeIndex.add(VertexIdEdgeVectorShape(inputEdges, inputVertices))
      if (options.splitCrossingEdges) {
          addEdgeCrossings(inputEdgeIndex)
      }
      if (snappingRequested) {
          val siteIndex = S2PointIndex<SiteId>()
          addForcedSites(siteIndex)
          chooseInitialSites(siteIndex)
          collectSiteEdges(siteIndex)
      }
      if (snappingNeeded) {
          addExtraSites(inputEdgeIndex)
      } else {
          copyInputEdges()
      }
  }

  private fun copyInputEdges(): Unit {
      // Sort the input vertices, discard duplicates, and update the input edges
      // to refer to the pruned vertex list.  (We sort in the same order used by
      // ChooseInitialSites() to avoid inconsistencies in tests.)
      val sorted = sortInputVertices()
      val vmap = mutableListOf<InputVertexId>()
      sites.clear()
      if (sites is ArrayList) {
          (sites as ArrayList<S2Point>).ensureCapacity(inputVertices.size)
      }
      var i = 0
      while (i < sorted.size) {
          val site = inputVertices[sorted[i].second]
          vmap[sorted[i].second] = sites.size
          while (++i < sorted.size && inputVertices[sorted[i].second] == site) {
          vmap[sorted[i].second] = sites.size
      }
          sites.add(site)
      }
      inputVertices = sites

      for (e in inputEdges.indices) {
          val edge = inputEdges[e]
          inputEdges[e] = Pair(vmap[edge.first], vmap[edge.second])
      }
  }

  private fun sortInputVertices(): List<InputVertexKey> = TODO()
  private fun addEdgeCrossings(input_edge_index: MutableS2ShapeIndex): Unit = TODO()
  private fun addForcedSites(site_index: S2PointIndex<SiteId>): Unit = TODO()
  private fun isForced(v: SiteId): Boolean =  v < num_forced_sites
  private fun chooseInitialSites(site_index: S2PointIndex<SiteId>): Unit = TODO()
  private fun snapSite(point: S2Point): S2Point = TODO()
  private fun collectSiteEdges(site_index: S2PointIndex<SiteId>): Unit = TODO()
  private fun sortSitesByDistance(x: S2Point, sites: MutableList<SiteId>): Unit = TODO()
  private fun addExtraSites(input_edge_index: MutableS2ShapeIndex): Unit = TODO()
  private fun maybeAddExtraSites(edge_id: InputEdgeId, max_edge_id: InputEdgeId, chain: List<SiteId>, 
                                 input_edge_index: MutableS2ShapeIndex, snap_queue: MutableList<InputEdgeId>): Unit = TODO()
  private fun addExtraSite(new_site: S2Point, max_edge_id: InputEdgeId, input_edge_index: MutableS2ShapeIndex,
                           snap_queue: MutableList<InputEdgeId>): Unit = TODO()
  private fun getSeparationSite(site_to_avoid: S2Point, v0: S2Point, v1: S2Point, input_edge_id: InputEdgeId): S2Point = TODO()
  private fun getCoverageEndpoint(p: S2Point, x: S2Point, y: S2Point, n: S2Point): S2Point = TODO()
  private fun snapEdge(e: InputEdgeId, chain: MutableList<SiteId>): Unit = TODO()

  private fun buildLayers(): Unit = TODO()
  private fun buildLayerEdges(
      layer_edges: MutableList<MutableList<Edge>>,
      layer_input_edge_ids: MutableList<MutableList<InputEdgeIdSetId>>,
      input_edge_id_set_lexicon: IdSetLexicon): Unit = TODO()
  private fun addSnappedEdges(begin: InputEdgeId, end: InputEdgeId, options: GraphOptions, edges: MutableList<Edge>,
                             input_edge_ids:  MutableList<InputEdgeIdSetId>, input_edge_id_set_lexicon: IdSetLexicon,
                              site_vertices: List<InputVertexId>): Unit = TODO()
  private fun maybeAddInputVertex(v: InputVertexId, id: SiteId, site_vertices: MutableList<MutableList<InputVertexId>>): Unit = TODO()
  private fun addSnappedEdge(src: SiteId, dst: SiteId, id: InputEdgeIdSetId, edge_type: EdgeType,
                             edges: MutableList<Edge>, input_edge_ids: MutableList<InputEdgeIdSetId>): Unit = TODO()
  private fun simplifyEdgeChains(site_vertices: MutableList<MutableList<InputVertexId>>, layer_edges: MutableList<MutableList<Edge>>,
                                 layer_input_edge_ids: MutableList<MutableList<InputEdgeIdSetId>>,
                                 input_edge_id_set_lexicon: IdSetLexicon): Unit = TODO()
  private fun mergeLayerEdges(layer_edges: MutableList<MutableList<Edge>>, layer_input_edge_ids: MutableList<MutableList<InputEdgeIdSetId>>,
                              edges: MutableList<Edge>, input_edge_ids: MutableList<InputEdgeIdSetId>,
                              edge_layers: MutableList<Int>): Unit = TODO()


    companion object {

        fun stableLessThan(a: Edge, b: Edge, ai: LayerEdgeId, bi: LayerEdgeId): Boolean = TODO()

        private fun roundUp(a: S1Angle): S1ChordAngle {
            val ca = S1ChordAngle(a)
            return ca.plusError(ca.getS1AngleConstructorMaxError())
        }

        private fun addPointToPointError(ca: S1ChordAngle): S1ChordAngle {
            return ca.plusError(ca.getS2PointConstructorMaxError())
        }

        private fun addPointToEdgeError(ca: S1ChordAngle): S1ChordAngle {
            return ca.plusError(S2EdgeDistances.getUpdateMinDistanceMaxError(ca))
        }
    }

}


