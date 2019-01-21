package org.pankratzlab.supernovo;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.Optional;
import java.util.Set;
import org.pankratzlab.supernovo.metrics.Depth;
import com.google.common.cache.CacheBuilder;
import com.google.common.cache.CacheLoader;
import com.google.common.cache.LoadingCache;
import com.google.common.collect.ImmutableList;
import com.google.common.collect.ImmutableSet;
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
  private static final CacheBuilder<Object, Object> PILEUP_CACHE_BUILDER =
      CacheBuilder.newBuilder().maximumSize(READ_LENGTH * 2L);

  private final String childID;
  private final String parent1ID;
  private final String parent2ID;

  private final LoadingCache<Position, Pileup> childPileups;
  private final LoadingCache<Position, Pileup> p1Pileups;
  private final LoadingCache<Position, Pileup> p2Pileups;

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
            CacheLoader.from(pos -> SAMRecordPileup.pileOverlapping(child, pos)));
    this.p1Pileups =
        PILEUP_CACHE_BUILDER.build(
            CacheLoader.from(pos -> SAMRecordPileup.pileOverlapping(parent1, pos)));
    this.p2Pileups =
        PILEUP_CACHE_BUILDER.build(
            CacheLoader.from(pos -> SAMRecordPileup.pileOverlapping(parent2, pos)));
  }

  public void reportDeNovos(VCFFileReader queriedVariants, File output) throws IOException {
    try (PrintWriter writer = new PrintWriter(new BufferedWriter(new FileWriter(output)))) {
      writer.println(DeNovoResult.Output.generateHeaderLine());
      queriedVariants
          .iterator()
          .stream()
          .filter(this::keepVariant)
          .map(Position::new)
          .map(this::evaluate)
          .filter(Optional::isPresent)
          .map(Optional::get)
          .map(DeNovoResult.Output::generateOutputLine)
          .forEachOrdered(writer::println);
    }
  }

  private boolean keepVariant(VariantContext vc) {
    Genotype geno = vc.getGenotype(childID);
    return geno.isHet()
        && geno.getAlleles().stream().mapToInt(Allele::length).allMatch(i -> i == 1);
  }

  private Optional<DeNovoResult> evaluate(Position pos) {
    Pileup childPile = childPileups.getUnchecked(pos);
    if (looksBiallelic(childPile) && looksDenovo(pos, childPile.getDepth())) {

      return Optional.of(
          new DeNovoResult(
              pos,
              generateSample(childID, childPile),
              generateSample(parent1ID, p1Pileups.getUnchecked(pos)),
              generateSample(parent2ID, p2Pileups.getUnchecked(pos))));
    }
    return Optional.empty();
  }

  private static boolean looksBiallelic(Pileup pileup) {
    Depth depth = pileup.getDepth();
    return depth.getBiAlleles().size() == 2
        && depth.weightedBiallelicDepth() >= MIN_DEPTH
        && pileup
            .getWeightedBaseCounts()
            .entrySet()
            .stream()
            .filter(e -> !depth.getBiAlleles().contains(e.getKey()))
            .mapToDouble(Map.Entry::getValue)
            .map(d -> d / depth.weightedTotalDepth())
            .allMatch(f -> f <= MAX_MISCALL_RATIO)
        && Arrays.stream(Depth.Allele.values())
            .mapToDouble(depth::allelicWeightedDepth)
            .allMatch(d -> d >= MIN_ALLELIC_DEPTH);
  }

  private static DeNovoResult.Sample generateSample(String id, Pileup pileup) {
    return new DeNovoResult.Sample(id, pileup.getDepth());
  }

  private boolean looksDenovo(Position position, Depth childDepth) {
    List<Pileup> parentPileups =
        ImmutableList.of(p1Pileups, p2Pileups)
            .stream()
            .map(c -> c.getUnchecked(position))
            .collect(ImmutableList.toImmutableList());
    Set<Byte> parentalAlleles =
        parentPileups
            .stream()
            .map(Pileup::getDepth)
            .map(Depth::getBiAlleles)
            .flatMap(Set::stream)
            .collect(ImmutableSet.toImmutableSet());
    return parentPileups.stream().allMatch(TrioEvaluator::looksBiallelic)
        && !Sets.difference(childDepth.getBiAlleles(), parentalAlleles).isEmpty();
  }
}
