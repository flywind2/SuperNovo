package org.pankratzlab.supernovo.pileup;

import com.google.common.collect.ImmutableList;
import com.google.common.collect.ImmutableMap;
import com.google.common.collect.ImmutableMultiset;
import com.google.common.collect.ImmutableSetMultimap;
import htsjdk.samtools.SAMRecord;

public interface Pileup {

  /**
   * @return Map from byte value of base to weighted depth for that base, iteration order is in
   *     descending order of weighted base counts
   */
  public ImmutableMap<Byte, Double> getWeightedBaseCounts();

  /** @return Multiset of byte values of bases */
  public ImmutableMultiset<Byte> getBaseCounts();

  /** @return Multimap from byte value of bases to index for the piled read */
  public ImmutableSetMultimap<Byte, Integer> getRecordsByBase();

  /** @return List of {@link SAMRecord}s piled up */
  public ImmutableList<SAMRecord> getRecords();

  /** @return the {@link Depth} for this {@link Pileup} */
  public Depth getDepth();
}
