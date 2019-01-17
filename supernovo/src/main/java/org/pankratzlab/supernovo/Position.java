package org.pankratzlab.supernovo;

import htsjdk.variant.variantcontext.VariantContext;

public class Position {

  private final String contig;
  private final int position;

  /**
   * @param contig
   * @param position
   */
  public Position(String contig, int position) {
    super();
    this.contig = contig;
    this.position = position;
  }

  public Position(VariantContext vc) {
    super();
    this.contig = vc.getContig();
    this.position = vc.getStart();
  }

  /**
   * @return the contig
   */
  public String getContig() {
    return contig;
  }

  /**
   * @return the position
   */
  public int getPosition() {
    return position;
  }

}
