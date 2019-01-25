package org.pankratzlab.supernovo;

import java.util.Optional;
import javax.annotation.Nullable;
import com.google.common.collect.ImmutableList;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;

public class ReferencePosition extends Position {

  final ImmutableList<Byte> refAllele;
  final Optional<ImmutableList<Byte>> altAllele;

  /**
   * @param contig
   * @param position
   * @param refAllele
   * @param altAllele
   */
  public ReferencePosition(
      String contig,
      int position,
      ImmutableList<Byte> refAllele,
      @Nullable ImmutableList<Byte> altAllele) {
    super(contig, position);
    this.refAllele = refAllele;
    this.altAllele = Optional.ofNullable(altAllele);
  }

  /**
   * @param contig
   * @param position
   * @param refAllele
   * @param altAllele
   */
  public ReferencePosition(
      String contig, int position, byte refAllele, @Nullable ImmutableList<Byte> altAllele) {
    this(contig, position, ImmutableList.of(refAllele), altAllele);
  }

  /**
   * @param contig
   * @param position
   * @param refAllele
   * @param altAllele
   */
  public ReferencePosition(
      String contig, int position, ImmutableList<Byte> refAllele, byte altAllele) {
    this(contig, position, refAllele, ImmutableList.of(altAllele));
  }

  /**
   * @param contig
   * @param position
   * @param refAllele
   * @param altAllele
   */
  public ReferencePosition(String contig, int position, byte refAllele, byte altAllele) {
    this(contig, position, refAllele, ImmutableList.of(altAllele));
  }

  /**
   * @param contig
   * @param position
   * @param refAllele
   */
  public ReferencePosition(String contig, int position, byte refAllele) {
    this(contig, position, refAllele, null);
  }

  public ReferencePosition(VariantContext vc, Allele ref, Allele alt) {
    this(vc.getContig(), vc.getStart(), alleleToBaseList(ref), alleleToBaseList(alt));
  }

  private static ImmutableList<Byte> alleleToBaseList(Allele allele) {
    ImmutableList.Builder<Byte> baseListBuilder =
        ImmutableList.builderWithExpectedSize(allele.length());
    for (byte base : allele.getBases()) {
      baseListBuilder.add(base);
    }
    return baseListBuilder.build();
  }

  /* (non-Javadoc)
   * @see java.lang.Object#hashCode()
   */
  @Override
  public int hashCode() {
    final int prime = 31;
    int result = 1;
    result = prime * result + ((altAllele == null) ? 0 : altAllele.hashCode());
    result = prime * result + ((contig == null) ? 0 : contig.hashCode());
    result = prime * result + position;
    result = prime * result + ((refAllele == null) ? 0 : refAllele.hashCode());
    return result;
  }

  /* (non-Javadoc)
   * @see java.lang.Object#equals(java.lang.Object)
   */
  @Override
  public boolean equals(Object obj) {
    if (this == obj) return true;
    if (obj == null) return false;
    if (!(obj instanceof ReferencePosition)) return false;
    ReferencePosition other = (ReferencePosition) obj;
    if (altAllele == null) {
      if (other.altAllele != null) return false;
    } else if (!altAllele.equals(other.altAllele)) return false;
    if (contig == null) {
      if (other.contig != null) return false;
    } else if (!contig.equals(other.contig)) return false;
    if (position != other.position) return false;
    if (refAllele == null) {
      if (other.refAllele != null) return false;
    } else if (!refAllele.equals(other.refAllele)) return false;
    return true;
  }
}
