package org.pankratzlab.supernovo.pileup;

import java.util.Comparator;
import java.util.Map;
import java.util.Optional;
import org.pankratzlab.supernovo.GenomePosition;
import org.pankratzlab.supernovo.utilities.Phred;
import com.google.common.collect.ImmutableList;
import com.google.common.collect.ImmutableMap;
import com.google.common.collect.ImmutableMultiset;
import com.google.common.collect.ImmutableSetMultimap;
import com.google.common.collect.Maps;
import htsjdk.samtools.SAMRecord;

public class Pileup {

  private final ImmutableSetMultimap<Byte, Integer> basePiles;
  private final ImmutableMap<Byte, Double> weightedBaseCounts;
  private final ImmutableList<SAMRecord> queriedRecords;
  private final ImmutableMultiset<Byte> clippedReadCounts;
  private final ImmutableMultiset<Byte> unmappedMateCounts;

  private Optional<Depth> depth = Optional.empty();

  public Pileup(ImmutableList<SAMRecord> queriedRecords, GenomePosition position) {
    super();
    ImmutableSetMultimap.Builder<Byte, Integer> basePilesBuilder = ImmutableSetMultimap.builder();
    Map<Byte, Double> weightedDepth = Maps.newHashMap();
    ImmutableMultiset.Builder<Byte> clippedReadCountsBuilder = ImmutableMultiset.builder();
    ImmutableMultiset.Builder<Byte> unmappedMateCountsBuilder = ImmutableMultiset.builder();
    for (int i = 0; i < queriedRecords.size(); i++) {
      SAMRecord samRecord = queriedRecords.get(i);
      int readPos = samRecord.getReadPositionAtReferencePosition(position.getPosition()) - 1;
      if (readPos != -1) {
        Byte base = samRecord.getReadBases()[readPos];
        basePilesBuilder.put(base, i);
        double accuracy =
            Phred.getAccuracy(samRecord.getBaseQualities()[readPos])
                * Phred.getAccuracy(samRecord.getMappingQuality());
        weightedDepth.put(base, weightedDepth.getOrDefault(base, 0.0) + accuracy);
        if (samRecord.getCigar().isClipped()) clippedReadCountsBuilder.add(base);
        if (samRecord.getMateUnmappedFlag()) unmappedMateCountsBuilder.add(base);
      }
    }
    basePiles = basePilesBuilder.build();
    weightedBaseCounts =
        ImmutableMap.<Byte, Double>builderWithExpectedSize(weightedDepth.size())
            .putAll(weightedDepth)
            .orderEntriesByValue(Comparator.reverseOrder())
            .build();
    clippedReadCounts = clippedReadCountsBuilder.build();
    unmappedMateCounts = unmappedMateCountsBuilder.build();
    this.queriedRecords = queriedRecords;
  }

  /** @return Multiset of byte values of bases */
  public ImmutableMultiset<Byte> getBaseCounts() {
    return basePiles.keys();
  }

  /**
   * @return Map from byte value of base to weighted depth for that base, iteration order is in
   *     descending order of weighted base counts
   */
  public ImmutableMap<Byte, Double> getWeightedBaseCounts() {
    return weightedBaseCounts;
  }

  /** @return Multimap from byte value of bases to index for the piled read */
  public ImmutableSetMultimap<Byte, Integer> getRecordsByBase() {
    return basePiles;
  }

  /** @return List of {@link SAMRecord}s piled up */
  public ImmutableList<SAMRecord> getRecords() {
    return queriedRecords;
  }

  /** @return the clippedReadCounts */
  public ImmutableMultiset<Byte> getClippedReadCounts() {
    return clippedReadCounts;
  }

  /** @return the unmappedMateCounts */
  public ImmutableMultiset<Byte> getUnmappedMateCounts() {
    return unmappedMateCounts;
  }

  /**
   * @return Map from byte value of base to weighted fraction of total weighted depth for that base,
   *     iteration order is in descending order of weighted base fraction
   */
  public ImmutableMap<Byte, Double> getWeightedBaseFractions() {
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
