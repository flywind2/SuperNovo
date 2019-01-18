package org.pankratzlab.supernovo;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Iterator;
import java.util.Optional;
import java.util.Set;
import com.google.common.cache.CacheBuilder;
import com.google.common.cache.CacheLoader;
import com.google.common.cache.LoadingCache;
import com.google.common.collect.ImmutableList;
import com.google.common.collect.ImmutableMultiset;
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
  private static final CacheBuilder<Object, Object> PILEUP_CACHE_BUILDER = CacheBuilder.newBuilder()
                                                                                       .maximumSize(READ_LENGTH
                                                                                                    * 2L);

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
  public TrioEvaluator(SamReader child, String childID, SamReader parent1, String parent1ID,
                       SamReader parent2, String parent2ID) {
    super();
    this.childID = childID;
    this.parent1ID = parent1ID;
    this.parent2ID = parent2ID;

    this.childPileups = PILEUP_CACHE_BUILDER.build(CacheLoader.from(pos -> new SAMRecordPileup(child,
                                                                                               pos)));
    this.p1Pileups = PILEUP_CACHE_BUILDER.build(CacheLoader.from(pos -> new SAMRecordPileup(parent1,
                                                                                            pos)));
    this.p2Pileups = PILEUP_CACHE_BUILDER.build(CacheLoader.from(pos -> new SAMRecordPileup(parent2,
                                                                                            pos)));
  }

  public void reportDeNovos(VCFFileReader queriedVariants, File output) throws IOException {
    try (PrintWriter writer = new PrintWriter(new BufferedWriter(new FileWriter(output)))) {
      writer.println(DeNovoResult.Output.generateHeaderLine());
      queriedVariants.iterator().stream().filter(this::keepVariant).map(Position::new)
                     .map(this::evaluate).filter(Optional::isPresent).map(Optional::get)
                     .map(DeNovoResult.Output::generateOutputLine).forEachOrdered(writer::println);
    }
  }

  private boolean keepVariant(VariantContext vc) {
    Genotype geno = vc.getGenotype(childID);
    return geno.isHet()
           && geno.getAlleles().stream().mapToInt(Allele::length).allMatch(i -> i == 1);
  }

  private Optional<DeNovoResult> evaluate(Position pos) {
    Pileup childPile = childPileups.getUnchecked(pos);
    Multiset<Byte> baseCounts = passingBaseCounts(childPile.getBaseCounts());
    if (baseCountBiallelicVariant(baseCounts) && looksDenovo(pos, baseCounts)) {
      Iterator<Byte> alleleIter = baseCounts.elementSet().iterator();
      Byte a1 = alleleIter.next();
      Byte a2 = alleleIter.next();
      return Optional.of(new DeNovoResult(pos, a1.byteValue(),
                                          a2.byteValue(),
                                          generateSample(childID, childPile, a1, a2),
                                          generateSample(parent1ID, p1Pileups.getUnchecked(pos), a1,
                                                         a2), generateSample(parent2ID, p2Pileups.getUnchecked(pos), a1,
                                                         a2)));
    }
    return Optional.empty();
  }

  private static boolean baseCountBiallelicVariant(Multiset<Byte> baseCounts) {
    return baseCounts.size() >= MIN_DEPTH && baseCounts.elementSet().size() == 2
           && baseCounts.entrySet().stream().mapToInt(Multiset.Entry::getCount)
                        .allMatch(i -> i >= MIN_ALLELIC_DEPTH);
  }

  private static DeNovoResult.Sample generateSample(String id, Pileup pileup, Byte a1, Byte a2) {
    Multiset<Byte> rawBaseCounts = pileup.getBaseCounts();
    return new DeNovoResult.Sample(id, rawBaseCounts.size(), rawBaseCounts.count(a1),
                                   rawBaseCounts.count(a2));
  }

  private boolean looksDenovo(Position position, Multiset<Byte> childBaseCounts) {
    return ImmutableList.of(p1Pileups, p2Pileups).stream().map(c -> c.getUnchecked(position))
                        .map(Pileup::getBaseCounts).map(TrioEvaluator::passingBaseCounts)
                        .map(Multiset::elementSet)
                        .map(p -> Sets.difference(childBaseCounts.elementSet(), p))
                        .noneMatch(Set::isEmpty);
  }

  private static Multiset<Byte> passingBaseCounts(Multiset<Byte> allBaseCounts) {
    double totalDepth = allBaseCounts.size();
    return allBaseCounts.entrySet().stream()
                        .filter(e -> e.getCount() / totalDepth > MAX_MISCALL_RATIO)
                        .collect(ImmutableMultiset.toImmutableMultiset(Multiset.Entry::getElement,
                                                                       Multiset.Entry::getCount));
  }

}
