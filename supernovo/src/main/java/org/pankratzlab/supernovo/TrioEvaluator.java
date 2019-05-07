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
import java.util.stream.IntStream;
import java.util.stream.Stream;
import org.pankratzlab.supernovo.output.DeNovoResult;
import org.pankratzlab.supernovo.output.OutputFields;
import org.pankratzlab.supernovo.pileup.Depth;
import org.pankratzlab.supernovo.pileup.Pileup;
import org.pankratzlab.supernovo.pileup.PileupCache;
import org.pankratzlab.supernovo.pileup.SAMPositionQueryOverlap;
import com.google.common.base.Predicates;
import com.google.common.cache.CacheBuilder;
import com.google.common.cache.CacheLoader;
import com.google.common.cache.LoadingCache;
import com.google.common.collect.ImmutableList;
import com.google.common.collect.ImmutableSet;
import com.google.common.collect.MoreCollectors;
import com.google.common.collect.Multiset;
import com.google.common.collect.RangeSet;
import com.google.common.collect.Sets;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFFileReader;

public class TrioEvaluator implements AutoCloseable {

  private static final int READ_LENGTH = 150;
  private static final int MIN_DEPTH = 10;
  private static final int MIN_ALLELIC_DEPTH = 4;
  private static final double MIN_ALLELIC_FRAC = 0.1;
  private static final double MAX_MISCALL_RATIO = 0.05;
  private static final double MAX_MISCALL_WEIGHT = 1.0;
  private static final CacheBuilder<Object, Object> PILEUP_CACHE_BUILDER =
      CacheBuilder.newBuilder().maximumSize(READ_LENGTH * 2L);

  private final String childID;
  private final String parent1ID;
  private final String parent2ID;

  private final SAMSequenceDictionary dict;

  private final PileupCache childPileups;
  private final LoadingCache<GenomePosition, Pileup> p1Pileups;
  private final LoadingCache<GenomePosition, Pileup> p2Pileups;

  /**
   * @param child {@link SamReader} of child to evluate for de novo variants
   * @param childReader2 TODO
   * @param parent1 {@link SamReader} of one parent for child
   * @param parent2 {@link SamReader} of second parent for child
   * @param intervals Genomic intervals to interrogate
   */
  public TrioEvaluator(
      SamReader child,
      SamReader childReader2,
      String childID,
      SamReader parent1,
      String parent1ID,
      SamReader parent2,
      String parent2ID,
      RangeSet<GenomePosition> intervals) {
    super();
    this.childID = childID;
    this.parent1ID = parent1ID;
    this.parent2ID = parent2ID;

    dict = child.getFileHeader().getSequenceDictionary();
    if (Stream.of(parent1, parent2)
        .map(SamReader::getFileHeader)
        .map(SAMFileHeader::getSequenceDictionary)
        .anyMatch(p -> !p.isSameDictionary(dict))) {
      throw new IllegalArgumentException(
          "Parent sequence dictionaries don't match child sequence dictionary");
    }

    this.childPileups = new PileupCache(child, childReader2, intervals);
    this.p1Pileups =
        PILEUP_CACHE_BUILDER.build(
            CacheLoader.from(
                pos -> new Pileup(new SAMPositionQueryOverlap(parent1, pos).getRecords(), pos)));
    this.p2Pileups =
        PILEUP_CACHE_BUILDER.build(
            CacheLoader.from(
                pos -> new Pileup(new SAMPositionQueryOverlap(parent2, pos).getRecords(), pos)));
  }

  public void reportAllDeNovos(IndexedFastaSequenceFile genome, File output) throws IOException {
    try (PrintWriter writer = new PrintWriter(new BufferedWriter(new FileWriter(output)))) {
      writer.println(OutputFields.generateHeader(DeNovoResult.class));
      childPileups.forEachRemaining(
          pileup ->
              evaluate(pileup)
                  .filter(d -> d.deNovo)
                  .map(DeNovoResult::generateLine)
                  .ifPresent(writer::println));
    }
  }

  private static Stream<ReferencePosition> generateRefPositions(ReferenceSequence refSeq) {
    Set<Byte> validBases = ImmutableSet.of((byte) 'A', (byte) 'T', (byte) 'C', (byte) 'G');
    byte[] bases = refSeq.getBases();
    String contig = refSeq.getName();
    return IntStream.range(0, bases.length)
        .filter(i -> validBases.contains(bases[i]))
        .mapToObj(i -> new ReferencePosition(contig, i + 1, SNPAllele.of(bases[i])));
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
    GenomePosition cleanup = new GenomePosition(pos.getContig(), pos.getPosition() - READ_LENGTH);
    childPileups
        .asMap()
        .keySet()
        .stream()
        .filter(p -> p.compareTo(cleanup) < 0)
        .forEach(childPileups::invalidate);
    if (looksVariant(childPile.getDepth())) {
      return Optional.of(
          new DeNovoResult(
              pos,
              new HaplotypeEvaluator(
                      childPile,
                      childPileups::getUnchecked,
                      p1Pileups::getUnchecked,
                      p2Pileups::getUnchecked)
                  .haplotypeConcordance(),
              generateSample(childID, pos, childPile, childPile),
              generateSample(parent1ID, pos, p1Pileups.getUnchecked(pos), childPile),
              generateSample(parent2ID, pos, p2Pileups.getUnchecked(pos), childPile)));
    }
    return Optional.empty();
  }

  private Optional<DeNovoResult> evaluate(Pileup childPile) {
    GenomePosition pos = childPile.getPosition();
    if (looksVariant(childPile.getDepth())) {
      return Optional.of(
          new DeNovoResult(
              pos,
              new HaplotypeEvaluator(
                      childPile,
                      childPileups::getUnchecked,
                      p1Pileups::getUnchecked,
                      p2Pileups::getUnchecked)
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
        && depth.weightedMinorAlleleFraction() >= MIN_ALLELIC_FRAC
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
      String id, GenomePosition pos, Pileup pileup, Pileup childPile) {
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

  @Override
  public void close() {
    childPileups.close();
  }
}
