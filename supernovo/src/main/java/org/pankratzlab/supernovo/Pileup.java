package org.pankratzlab.supernovo;

import org.pankratzlab.supernovo.SAMRecordPileup.PiledRecord;
import com.google.common.collect.ImmutableMap;
import com.google.common.collect.ImmutableMultimap;
import com.google.common.collect.ImmutableMultiset;

public interface Pileup {

  public ImmutableMap<Byte, Double> getWeightedBaseCounts();

  public ImmutableMultiset<Byte> getBaseCounts();

  public ImmutableMultimap<Byte, PiledRecord> getRecordsByBase();
}
