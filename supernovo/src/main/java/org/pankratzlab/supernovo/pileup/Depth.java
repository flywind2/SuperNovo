package org.pankratzlab.supernovo.pileup;

import java.util.Iterator;
import java.util.Optional;
import java.util.Set;
import java.util.function.Function;
import org.pankratzlab.supernovo.PileAllele;
import com.google.common.collect.ImmutableSet;

public class Depth {

  public enum Allele {
    A1(Depth::getA1),
    A2(Depth::getA2);

    private final Function<Depth, Optional<PileAllele>> getterFunc;

    /** @param getterFunc */
    private Allele(Function<Depth, Optional<PileAllele>> getterFunc) {
      this.getterFunc = getterFunc;
    }

    private Optional<PileAllele> getAllele(Depth depth) {
      return getterFunc.apply(depth);
    }
  }

  private final Pileup pileup;
  private final Optional<PileAllele> a1;
  private final Optional<PileAllele> a2;
  private final Set<PileAllele> biAlleles;

  /** @param pileup */
  public Depth(Pileup pileup) {
    super();
    this.pileup = pileup;
    Iterator<PileAllele> alleleIter = pileup.getWeightedBaseCounts().keySet().iterator();
    a1 = Optional.ofNullable(alleleIter.hasNext() ? alleleIter.next() : null);
    a2 = Optional.ofNullable(alleleIter.hasNext() ? alleleIter.next() : null);
    ImmutableSet.Builder<PileAllele> allelesBuilder = ImmutableSet.builderWithExpectedSize(2);
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

  public double rawMinorAlleleFraction() {
    return biAlleles.stream().mapToDouble(pileup.getBaseFractions()::get).min().orElse(0.0);
  }

  public double weightedMinorAlleleFraction() {
    return biAlleles.stream().mapToDouble(pileup.getWeightedBaseFractions()::get).min().orElse(0.0);
  }

  /** @return the a1 */
  public Optional<PileAllele> getA1() {
    return a1;
  }

  /** @return the a2 */
  public Optional<PileAllele> getA2() {
    return a2;
  }

  /** @return the biAlleles */
  public Set<PileAllele> getBiAlleles() {
    return biAlleles;
  }

  public double allelicWeightedDepth(PileAllele allele) {
    return pileup.getWeightedBaseCounts().getOrDefault(allele, 0.0);
  }

  public double allelicWeightedDepth(Allele allele) {
    return allele
        .getAllele(this)
        .map(this::allelicWeightedDepth)
        .orElse(Double.valueOf(0))
        .doubleValue();
  }

  public int allelicRawDepth(PileAllele allele) {
    return pileup.getBaseCounts().count(allele);
  }

  public int allelicRawDepth(Allele allele) {
    return allele.getAllele(this).map(this::allelicRawDepth).orElse(Integer.valueOf(0)).intValue();
  }

  public ImmutableSet<Integer> allelicRecords(Allele allele) {
    return allele.getAllele(this).map(pileup.getRecordsByBase()::get).orElse(ImmutableSet.of());
  }
}
