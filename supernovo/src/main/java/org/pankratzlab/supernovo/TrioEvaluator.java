package org.pankratzlab.supernovo;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Arrays;
import java.util.List;
import java.util.Optional;
import java.util.Set;
import org.pankratzlab.supernovo.output.DeNovoResult;
import org.pankratzlab.supernovo.output.OutputFields;
import org.pankratzlab.supernovo.pileup.Depth;
import org.pankratzlab.supernovo.pileup.Pileup;
import org.pankratzlab.supernovo.pileup.SAMPositionOverlap;
import com.google.common.base.Predicates;
import com.google.common.cache.CacheBuilder;
import com.google.common.cache.CacheLoader;
import com.google.common.cache.LoadingCache;
import com.google.common.collect.ImmutableList;
import com.google.common.collect.ImmutableSet;
import com.google.common.collect.MoreCollectors;
import com.google.common.collect.Multiset;
import com.google.common.collect.Sets;
import htsjdk.samtools.SamReader;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;

public class TrioEvaluator {

  private static final int READ_LENGTH = 150;
  private static final int MIN_DEPTH = 10;
  private static final int MIN_ALLELIC_DEPTH = 4;
  private static final double MAX_MISCALL_RATIO = 0.05;
  private static final double MAX_MISCALL_WEIGHT = 1.0;
  private static final CacheBuilder<Object, Object> PILEUP_CACHE_BUILDER =
      CacheBuilder.newBuilder().maximumSize(READ_LENGTH * 2L);

  private final String childID;
  private final String parent1ID;
  private final String parent2ID;

  private final LoadingCache<ReferencePosition, Pileup> childPileups;
  private final LoadingCache<ReferencePosition, Pileup> p1Pileups;
  private final LoadingCache<ReferencePosition, Pileup> p2Pileups;

  /**
   * @param child {@link SamReader} of child to evluate for de novo variants
   * @param parent1 {@link SamReader} of one parent for child
   * @param parent2 {@link SamReader} of second parent for child
   */
  public TrioEvaluator(
      SamReader child,
      String childID,
      SamReader parent1,
      String parent1ID,
      SamReader parent2,
      String parent2ID) {
    super();
    this.childID = childID;
    this.parent1ID = parent1ID;
    this.parent2ID = parent2ID;

    this.childPileups =
        PILEUP_CACHE_BUILDER.build(
            CacheLoader.from(
                pos -> new Pileup(new SAMPositionOverlap(child, pos).getRecords(), pos)));
    this.p1Pileups =
        PILEUP_CACHE_BUILDER.build(
            CacheLoader.from(
                pos -> new Pileup(new SAMPositionOverlap(parent1, pos).getRecords(), pos)));
    this.p2Pileups =
        PILEUP_CACHE_BUILDER.build(
            CacheLoader.from(
                pos -> new Pileup(new SAMPositionOverlap(parent2, pos).getRecords(), pos)));
  }

  public void reportDeNovos(VCFFileReader queriedVariants, File output) throws IOException {
    try (PrintWriter writer = new PrintWriter(new BufferedWriter(new FileWriter(output)))) {
      writer.println(OutputFields.generateHeader(DeNovoResult.class));
      queriedVariants
          .iterator()
          .stream()
          .filter(this::keepVariant)
          .map(this::generatePosition)
          .map(this::evaluate)
          .filter(Optional::isPresent)
          .map(Optional::get)
          .map(DeNovoResult::generateLine)
          .forEachOrdered(writer::println);
    }
  }

  private boolean keepVariant(VariantContext vc) {
    Genotype geno = vc.getGenotype(childID);
    return geno.isHet()
        && !geno.isHetNonRef()
        && geno.getAlleles().stream().mapToInt(Allele::length).anyMatch(i -> i == 1);
  }

  private ReferencePosition generatePosition(VariantContext vc) {
    Allele ref = vc.getReference();
    Genotype geno = vc.getGenotype(childID);
    Allele alt =
        geno.getAlleles()
            .stream()
            .filter(Predicates.not(vc.getReference()::equals))
            .collect(MoreCollectors.onlyElement());
    return ReferencePosition.fromVariantContext(vc, ref, alt);
  }

  private Optional<DeNovoResult> evaluate(ReferencePosition pos) {
    Pileup childPile = childPileups.getUnchecked(pos);
    if (looksVariant(childPile.getDepth())) {
      return Optional.of(
          new DeNovoResult(
              pos,
              new HaplotypeEvaluator(
                      pos, childPile, p1Pileups.getUnchecked(pos), p2Pileups.getUnchecked(pos))
                  .haplotypeConcordance(),
              generateSample(childID, pos, childPile, childPile),
              generateSample(parent1ID, pos, p1Pileups.getUnchecked(pos), childPile),
              generateSample(parent2ID, pos, p2Pileups.getUnchecked(pos), childPile)));
    }
    return Optional.empty();
  }

  public static boolean looksBiallelic(Pileup pileup) {
    return looksVariant(pileup.getDepth()) && !moreThanTwoViableAlleles(pileup);
  }

  public static boolean looksVariant(Depth depth) {
    return depth.getBiAlleles().size() == 2
        && depth.weightedBiallelicDepth() >= MIN_DEPTH
        && Arrays.stream(Depth.Allele.values())
            .mapToDouble(depth::allelicWeightedDepth)
            .allMatch(d -> d >= MIN_ALLELIC_DEPTH);
  }

  public static boolean moreThanTwoViableAlleles(Pileup pileup) {
    return possibleAlleles(pileup).size() > 2;
  }

  private static Set<PileAllele> possibleAlleles(Pileup pileup) {
    return pileup
        .getBaseCounts()
        .entrySet()
        .stream()
        .filter(e -> e.getCount() > MAX_MISCALL_WEIGHT)
        .map(Multiset.Entry::getElement)
        .filter(b -> pileup.getBaseFractions().get(b) > MAX_MISCALL_RATIO)
        .collect(ImmutableSet.toImmutableSet());
  }

  private static DeNovoResult.Sample generateSample(
      String id, ReferencePosition pos, Pileup pileup, Pileup childPile) {
    return new DeNovoResult.Sample(
        id, pileup, pos, childPile.getDepth().getA1(), childPile.getDepth().getA2());
  }

  public static boolean looksDenovo(Pileup childPileup, Pileup p1Pileup, Pileup p2Pileup) {
    List<Pileup> parentPileups = ImmutableList.of(p1Pileup, p2Pileup);
    Set<PileAllele> parentalAlleles =
        parentPileups
            .stream()
            .map(TrioEvaluator::possibleAlleles)
            .flatMap(Set::stream)
            .collect(ImmutableSet.toImmutableSet());
    return !Sets.difference(childPileup.getDepth().getBiAlleles(), parentalAlleles).isEmpty();
  }
}
