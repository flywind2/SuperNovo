package org.pankratzlab.supernovo;

import java.util.Collection;
import java.util.Comparator;
import java.util.Map;
import java.util.stream.Collector;
import org.pankratzlab.supernovo.utilities.Phred;
import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.ImmutableMap;
import com.google.common.collect.ImmutableMultiset;
import com.google.common.collect.ImmutableSetMultimap;
import com.google.common.collect.ListMultimap;
import com.google.common.collect.Maps;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;

public class SAMRecordPileup extends AbstractPileup {

  private final ImmutableSetMultimap<Byte, PiledRecord> basePiles;
  private final ImmutableMap<Byte, Double> weightedBaseCounts;

  public SAMRecordPileup(SamReader samReader, Position position) {
    super();
    try (SAMRecordIterator iterator =
        samReader.queryOverlapping(
            position.getContig(), position.getPosition(), position.getPosition())) {
      ImmutableSetMultimap.Builder<Byte, PiledRecord> basePilesBuilder =
          ImmutableSetMultimap.builder();
      ListMultimap<Byte, Integer> basePhreds = ArrayListMultimap.create();
      while (iterator.hasNext()) {
        SAMRecord samRecord = iterator.next();
        int readPos = samRecord.getReadPositionAtReferencePosition(position.getPosition()) - 1;
        if (readPos != -1) {
          Byte base = samRecord.getReadBases()[readPos];
          basePilesBuilder.put(base, new PiledRecord(samRecord));
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
    }
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
  public ImmutableSetMultimap<Byte, PiledRecord> getRecordsByBase() {
    return basePiles;
  }
}
