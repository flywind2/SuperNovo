package org.pankratzlab.supernovo;

import htsjdk.samtools.SAMRecord;

public interface PileAllele {

  /**
   * @param record {@link SAMRecord} to test
   * @return true if this read supports this {@link PileAllele}
   */
  boolean supported(SAMRecord record, int readPos);

  /** @return String representation of the allele */
  @Override
  String toString();
}
