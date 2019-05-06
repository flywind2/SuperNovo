package org.pankratzlab.supernovo.pileup;

import java.util.Comparator;
import java.util.List;
import java.util.Map;
import java.util.Optional;
import org.pankratzlab.supernovo.GenomePosition;
import org.pankratzlab.supernovo.PileAllele;
import org.pankratzlab.supernovo.ReferencePosition;
import org.pankratzlab.supernovo.SNPAllele;
import com.google.common.collect.ImmutableList;
import com.google.common.collect.ImmutableMap;
import com.google.common.collect.ImmutableMultiset;
import com.google.common.collect.ImmutableSetMultimap;
import com.google.common.collect.Maps;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.SAMRecord;

public class Pileup {

  private static final double MIN_PERCENT_BASES_MATCH = 0.5;

  private final ImmutableSetMultimap<PileAllele, Integer> basePiles;
  private final ImmutableMap<PileAllele, Double> weightedBaseCounts;
  private final ImmutableList<SAMRecord> queriedRecords;
  private final ImmutableMultiset<PileAllele> clippedReadCounts;
  private final ImmutableMultiset<PileAllele> apparentMismapReadCounts;
  private final ImmutableMultiset<PileAllele> unmappedMateCounts;

  private Optional<Depth> depth = Optional.empty();

  public Pileup(ImmutableList<SAMRecord> queriedRecords, GenomePosition position) {
    super();
    List<PileAllele> queriedAlleles = generateQueriedAlleles(position);
    ImmutableSetMultimap.Builder<PileAllele, Integer> basePilesBuilder =
        ImmutableSetMultimap.builder();
    Map<PileAllele, Double> weightedDepth = Maps.newHashMap();
    ImmutableMultiset.Builder<PileAllele> clippedReadCountsBuilder = ImmutableMultiset.builder();
    ImmutableMultiset.Builder<PileAllele> apparentMismapReadCountsBuilder =
        ImmutableMultiset.builder();
    ImmutableMultiset.Builder<PileAllele> unmappedMateCountsBuilder = ImmutableMultiset.builder();
    for (int i = 0; i < queriedRecords.size(); i++) {
      SAMRecord samRecord = queriedRecords.get(i);
      int readPos = samRecord.getReadPositionAtReferencePosition(position.getPosition()) - 1;
      if (readPos != -1) {
        PileAllele allele =
            queriedAlleles
                .stream()
                .filter(a -> a.supported(samRecord, readPos))
                .findFirst()
                .orElseGet(() -> getAppropriateAllele(samRecord, readPos));
        basePilesBuilder.put(allele, i);
        boolean countWeight = true;
        if (samRecord.getCigar().isClipped()) {
          clippedReadCountsBuilder.add(allele);
          countWeight = false;
        } else if (calcPercentReadMatchesRef(samRecord) >= MIN_PERCENT_BASES_MATCH) {
          apparentMismapReadCountsBuilder.add(allele);
          countWeight = false;
        }
        if (samRecord.getMateUnmappedFlag()) {
          unmappedMateCountsBuilder.add(allele);
          countWeight = false;
        }
        if (countWeight) {
          weightedDepth.put(
              allele,
              weightedDepth.getOrDefault(allele, 0.0) + allele.weightedDepth(samRecord, readPos));
        }
      }
    }
    basePiles = basePilesBuilder.build();
    weightedBaseCounts =
        ImmutableMap.<PileAllele, Double>builderWithExpectedSize(weightedDepth.size())
            .putAll(weightedDepth)
            .orderEntriesByValue(Comparator.reverseOrder())
            .build();
    clippedReadCounts = clippedReadCountsBuilder.build();
    apparentMismapReadCounts = apparentMismapReadCountsBuilder.build();
    unmappedMateCounts = unmappedMateCountsBuilder.build();
    this.queriedRecords = queriedRecords;
  }

  private double calcPercentReadMatchesRef(SAMRecord samRecord) {
    return samRecord
            .getCigar()
            .getCigarElements()
            .stream()
            .filter(c -> c.getOperator().equals(CigarOperator.EQ))
            .mapToInt(CigarElement::getLength)
            .sum()
        / (double) samRecord.getReadLength();
  }

  private static PileAllele getAppropriateAllele(SAMRecord samRecord, int readPos) {
    byte base = samRecord.getReadBases()[readPos];
    return SNPAllele.of(base);
  }

  private static ImmutableList<PileAllele> generateQueriedAlleles(GenomePosition pos) {
    if (!(pos instanceof ReferencePosition)) return ImmutableList.of();
    ReferencePosition refPos = (ReferencePosition) pos;
    ImmutableList.Builder<PileAllele> queriedAllelesBuilder =
        ImmutableList.builderWithExpectedSize(2);
    queriedAllelesBuilder.add(refPos.getRefAllele());
    refPos.getAltAllele().ifPresent(queriedAllelesBuilder::add);
    return queriedAllelesBuilder.build();
  }

  /** @return Multiset of {@link PileAllele} counts */
  public ImmutableMultiset<PileAllele> getBaseCounts() {
    return basePiles.keys();
  }

  /**
   * @return Map from {@link PileAllele} to fraction of total depth for that {@link PileAllele},
   *     iteration order is in descending order of base fraction
   */
  public ImmutableMap<PileAllele, Double> getBaseFractions() {

    return ImmutableMap.copyOf(
        Maps.<PileAllele, Double>asMap(
            getBaseCounts().elementSet(),
            b -> getBaseCounts().count(b) / (double) getBaseCounts().size()));
  }

  /**
   * @return Map from {@link PileAllele} to weighted depth for that PileAllele, iteration order is
   *     in descending order of weighted base counts
   */
  public ImmutableMap<PileAllele, Double> getWeightedBaseCounts() {
    return weightedBaseCounts;
  }

  /** @return Multimap from {@link PileAllele} to index for the piled read */
  public ImmutableSetMultimap<PileAllele, Integer> getRecordsByBase() {
    return basePiles;
  }

  /** @return List of {@link SAMRecord}s piled up */
  public ImmutableList<SAMRecord> getRecords() {
    return queriedRecords;
  }

  /** @return the clippedReadCounts */
  public ImmutableMultiset<PileAllele> getClippedReadCounts() {
    return clippedReadCounts;
  }

  /** @return the apparentMismapReadCounts */
  public ImmutableMultiset<PileAllele> getApparentMismapReadCounts() {
    return apparentMismapReadCounts;
  }

  /** @return the unmappedMateCounts */
  public ImmutableMultiset<PileAllele> getUnmappedMateCounts() {
    return unmappedMateCounts;
  }

  /**
   * @return Map from {@link PileAllele} to weighted fraction of total weighted depth for that
   *     {@link PileAllele}, iteration order is in descending order of weighted base fraction
   */
  public ImmutableMap<PileAllele, Double> getWeightedBaseFractions() {
    final double weightedDepth =
        weightedBaseCounts.values().stream().mapToDouble(Double::valueOf).sum();
    return ImmutableMap.copyOf(Maps.transformValues(weightedBaseCounts, c -> c / weightedDepth));
  }

  private Depth setDepth() {
    depth = Optional.of(new Depth(this));
    return depth.get();
  }

  public Depth getDepth() {
    return depth.orElseGet(this::setDepth);
  }

  @Override
  public String toString() {
    return "Pileup [basePiles="
        + basePiles
        + ", weightedBaseCounts="
        + weightedBaseCounts
        + ", clippedReadCounts="
        + clippedReadCounts
        + ", apparentMismapReadCounts="
        + apparentMismapReadCounts
        + ", unmappedMateCounts="
        + unmappedMateCounts
        + ", position="
        + position
        + ", depth="
        + depth
        + "]";
  }
}
