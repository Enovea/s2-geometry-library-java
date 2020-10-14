package dilivia.s2.builder.layers

import dilivia.s2.S2Error
import dilivia.s2.builder.Graph
import dilivia.s2.builder.GraphOptions


// This class is not needed by ordinary S2Builder clients.  It is only
// necessary if you wish to implement a new S2Builder::Layer subtype.
abstract class Layer {

  // Defines options for building the edge graph that is passed to Build().
  abstract fun graphOptions(): GraphOptions

  // Assembles a graph of snapped edges into the geometry type implemented by
  // this layer.  If an error is encountered, sets "error" appropriately.
  //
  // Note that when there are multiple layers, the Graph objects passed to all
  // layers are guaranteed to be valid until the last Build() method returns.
  // This makes it easier to write algorithms that gather the output graphs
  // from several layers and process them all at once (such as
  // s2builderutil::ClosedSetNormalizer).
  abstract fun build(g: Graph, error: S2Error)

}

