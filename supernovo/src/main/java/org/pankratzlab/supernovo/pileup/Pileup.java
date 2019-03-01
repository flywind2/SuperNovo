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
import htsjdk.samtools.SAMRecord;

public class Pileup {

  private final ImmutableSetMultimap<PileAllele, Integer> basePiles;
  private final ImmutableMap<PileAllele, Double> weightedBaseCounts;
  private final ImmutableList<SAMRecord> queriedRecords;
  private final ImmutableMultiset<PileAllele> clippedReadCounts;
  private final ImmutableMultiset<PileAllele> unmappedMateCounts;

  private Optional<Depth> depth = Optional.empty();

  public Pileup(ImmutableList<SAMRecord> queriedRecords, GenomePosition position) {
    super();
    List<PileAllele> queriedAlleles = generateQueriedAlleles(position);
    ImmutableSetMultimap.Builder<PileAllele, Integer> basePilesBuilder =
        ImmutableSetMultimap.builder();
    Map<PileAllele, Double> weightedDepth = Maps.newHashMap();
    ImmutableMultiset.Builder<PileAllele> clippedReadCountsBuilder = ImmutableMultiset.builder();
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
        if (!allele.clipped(samRecord, readPos))
          weightedDepth.put(
              allele,
              weightedDepth.getOrDefault(allele, 0.0) + allele.weightedDepth(samRecord, readPos));
        if (samRecord.getCigar().isClipped()) clippedReadCountsBuilder.add(allele);
        if (samRecord.getMateUnmappedFlag()) unmappedMateCountsBuilder.add(allele);
      }
    }
    basePiles = basePilesBuilder.build();
    weightedBaseCounts =
        ImmutableMap.<PileAllele, Double>builderWithExpectedSize(weightedDepth.size())
            .putAll(weightedDepth)
            .orderEntriesByValue(Comparator.reverseOrder())
            .build();
    clippedReadCounts = clippedReadCountsBuilder.build();
    unmappedMateCounts = unmappedMateCountsBuilder.build();
    this.queriedRecords = queriedRecords;
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
}
