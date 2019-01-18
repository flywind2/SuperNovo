package org.pankratzlab.supernovo;

import com.google.common.collect.Range;
import htsjdk.variant.variantcontext.VariantContext;

public class Segment {

  private final String contig;
  private final Range<Integer> range;
  private final int start;
  private final int stop;

  /**
   * @param contig contig this {@link Segment} resides on
   * @param start 1-based inclusive start position
   * @param stop 1-based closed stop position
   */
  public Segment(String contig, int start, int stop) {
    this.contig = contig;
    range = Range.closed(start, stop);
    this.start = start;
    this.stop = stop;
  }

  public Segment(String contig, int start) {
    this(contig, start, start);
  }

  public Segment(VariantContext vc) {
    this(vc.getContig(), vc.getStart(), vc.getEnd());
  }
}
