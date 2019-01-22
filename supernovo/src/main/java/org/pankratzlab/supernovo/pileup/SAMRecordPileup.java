package org.pankratzlab.supernovo.pileup;

import java.util.Collection;
import java.util.Comparator;
import java.util.Map;
import java.util.stream.Collector;
import org.pankratzlab.supernovo.Position;
import org.pankratzlab.supernovo.utilities.Phred;
import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.ImmutableList;
import com.google.common.collect.ImmutableMap;
import com.google.common.collect.ImmutableMultiset;
import com.google.common.collect.ImmutableSetMultimap;
import com.google.common.collect.ListMultimap;
import com.google.common.collect.Maps;
import htsjdk.samtools.SAMRecord;

public class SAMRecordPileup extends AbstractPileup {

  private final ImmutableSetMultimap<Byte, Integer> basePiles;
  private final ImmutableMap<Byte, Double> weightedBaseCounts;
  private final ImmutableList<SAMRecord> queriedRecords;

  public SAMRecordPileup(ImmutableList<SAMRecord> queriedRecords, Position position) {
    super();
    ImmutableSetMultimap.Builder<Byte, Integer> basePilesBuilder = ImmutableSetMultimap.builder();
    ListMultimap<Byte, Integer> basePhreds = ArrayListMultimap.create();
    for (int i = 0; i < queriedRecords.size(); i++) {
      SAMRecord samRecord = queriedRecords.get(i);
      int readPos = samRecord.getReadPositionAtReferencePosition(position.getPosition()) - 1;
      if (readPos != -1) {
        Byte base = samRecord.getReadBases()[readPos];
        basePilesBuilder.put(base, i);
        basePhreds.put(base, Integer.valueOf(samRecord.getBaseQualities()[readPos]));
      }
    }
    basePiles = basePilesBuilder.build();
    weightedBaseCounts =
        basePhreds
            .asMap()
            .entrySet()
            .stream()
            .collect(
                Collector.of(
                    () ->
                        ImmutableMap.<Byte, Double>builder()
                            .orderEntriesByValue(Comparator.reverseOrder()),
                    (b, e) -> b.put(phredScoresToWeightedDepth(e)),
                    (b1, b2) -> b1.putAll(b2.build()),
                    ImmutableMap.Builder::build));
    this.queriedRecords = queriedRecords;
  }

  private static <K> Map.Entry<K, Double> phredScoresToWeightedDepth(
      Map.Entry<K, Collection<Integer>> phredScores) {
    return Maps.immutableEntry(
        phredScores.getKey(),
        phredScores.getValue().stream().mapToDouble(Phred::getAccuracy).sum());
  }

  @Override
  public ImmutableMultiset<Byte> getBaseCounts() {
    return basePiles.keys();
  }

  @Override
  public ImmutableMap<Byte, Double> getWeightedBaseCounts() {
    return weightedBaseCounts;
  }

  @Override
  public ImmutableSetMultimap<Byte, Integer> getRecordsByBase() {
    return basePiles;
  }

  @Override
  public ImmutableList<SAMRecord> getRecords() {
    return queriedRecords;
  }
}
