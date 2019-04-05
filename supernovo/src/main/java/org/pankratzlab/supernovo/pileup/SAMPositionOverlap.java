package org.pankratzlab.supernovo.pileup;

import com.google.common.collect.ImmutableList;
import htsjdk.samtools.SAMRecord;

public interface SAMPositionOverlap {

  /** @return the records */
  ImmutableList<SAMRecord> getRecords();
}