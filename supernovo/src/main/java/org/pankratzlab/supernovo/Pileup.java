package org.pankratzlab.supernovo;

import org.pankratzlab.supernovo.SAMRecordPileup.PiledRecord;
import com.google.common.collect.ImmutableMultimap;
import com.google.common.collect.ImmutableMultiset;

public interface Pileup {

  public ImmutableMultiset<Byte> getBaseCounts();

  public ImmutableMultimap<Byte, PiledRecord> getRecordsByBase();

}
