package org.pankratzlab.supernovo.metrics;

import java.util.Iterator;
import java.util.Optional;
import java.util.Set;
import java.util.function.Function;
import org.pankratzlab.supernovo.Pileup;
import com.google.common.collect.ImmutableSet;

public class Depth {

  public enum Allele {
    A1(Depth::getA1),
    A2(Depth::getA2);

    private final Function<Depth, Optional<Byte>> getterFunc;

    /** @param getterFunc */
    private Allele(Function<Depth, Optional<Byte>> getterFunc) {
      this.getterFunc = getterFunc;
    }

    private Optional<Byte> getAllele(Depth depth) {
      return getterFunc.apply(depth);
    }
  }

  private final Pileup pileup;
  private final Optional<Byte> a1;
  private final Optional<Byte> a2;
  private final Set<Byte> biAlleles;

  /** @param pileup */
  public Depth(Pileup pileup) {
    super();
    this.pileup = pileup;
    Iterator<Byte> alleleIter = pileup.getWeightedBaseCounts().keySet().iterator();
    a1 = Optional.ofNullable(alleleIter.hasNext() ? alleleIter.next() : null);
    a2 = Optional.ofNullable(alleleIter.hasNext() ? alleleIter.next() : null);
    ImmutableSet.Builder<Byte> allelesBuilder = ImmutableSet.builderWithExpectedSize(2);
    a1.ifPresent(allelesBuilder::add);
    a2.ifPresent(allelesBuilder::add);
    biAlleles = allelesBuilder.build();
  }

  public double weightedBiallelicDepth() {
    return biAlleles.stream().mapToDouble(pileup.getWeightedBaseCounts()::get).sum();
  }

  public double weightedTotalDepth() {
    return pileup.getWeightedBaseCounts().values().stream().mapToDouble(Double::doubleValue).sum();
  }

  public int rawBiallelicDepth() {
    return biAlleles.stream().mapToInt(pileup.getBaseCounts()::count).sum();
  }

  public int rawTotalDepth() {
    return pileup.getBaseCounts().size();
  }

  /** @return the a1 */
  public Optional<Byte> getA1() {
    return a1;
  }

  /** @return the a2 */
  public Optional<Byte> getA2() {
    return a2;
  }

  /** @return the biAlleles */
  public Set<Byte> getBiAlleles() {
    return biAlleles;
  }

  public double allelicWeightedDepth(byte allele) {
    return pileup.getWeightedBaseCounts().get(allele);
  }

  public double allelicWeightedDepth(Allele allele) {
    return allele
        .getAllele(this)
        .map(this::allelicWeightedDepth)
        .orElse(Double.valueOf(0))
        .doubleValue();
  }

  public int allelicRawDepth(byte allele) {
    return pileup.getBaseCounts().count(allele);
  }

  public int allelicRawDepth(Allele allele) {
    return allele.getAllele(this).map(this::allelicRawDepth).orElse(Integer.valueOf(0)).intValue();
  }
}
