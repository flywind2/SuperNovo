package org.pankratzlab.supernovo;

import java.util.Optional;
import htsjdk.variant.variantcontext.VariantContext;

public class Position implements Comparable<Position> {

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

  /** @return the contig */
  public String getContig() {
    return contig;
  }

  /** @return the position */
  public int getPosition() {
    return position;
  }

  private Optional<Integer> locateChrNumber() {
    int i;
    for (i = 0; i < contig.length(); i++) {
      if (Character.isDigit(contig.charAt(i))) break;
    }
    StringBuilder digitBuilder = new StringBuilder();
    for (; i < contig.length(); i++) {
      char charAt = contig.charAt(i);
      if (Character.isDigit(charAt)) digitBuilder.append(charAt);
    }
    try {
      return Optional.of(Integer.parseInt(digitBuilder.toString()));
    } catch (NumberFormatException e) {
      return Optional.empty();
    }
  }

  @Override
  public int compareTo(Position o) {
    if (contig.equals(contig)) return Integer.compare(position, o.position);
    int cmpChrNum =
        locateChrNumber()
            .orElse(Integer.MAX_VALUE)
            .compareTo(o.locateChrNumber().orElse(Integer.MAX_VALUE));
    if (cmpChrNum != 0) return cmpChrNum;
    return contig.compareTo(o.contig);
  }

  /* (non-Javadoc)
   * @see java.lang.Object#hashCode()
   */
  @Override
  public int hashCode() {
    final int prime = 31;
    int result = 1;
    result = prime * result + ((contig == null) ? 0 : contig.hashCode());
    result = prime * result + position;
    return result;
  }

  /* (non-Javadoc)
   * @see java.lang.Object#equals(java.lang.Object)
   */
  @Override
  public boolean equals(Object obj) {
    if (this == obj) return true;
    if (obj == null) return false;
    if (!(obj instanceof Position)) return false;
    Position other = (Position) obj;
    if (contig == null) {
      if (other.contig != null) return false;
    } else if (!contig.equals(other.contig)) return false;
    if (position != other.position) return false;
    return true;
  }
}
