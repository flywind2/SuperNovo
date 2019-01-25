package org.pankratzlab.supernovo;

import htsjdk.samtools.SAMRecord;

public interface PileAllele {

  /**
   * @param record to test
   * @param readPos to query
   * @return true if this read supports this {@link PileAllele}
   */
  boolean supported(SAMRecord record, int readPos);

  /**
   * @param record to test
   * @param readPos to query
   * @return weighted depth for allele
   */
  double weightedDepth(SAMRecord record, int readPos);

  /** @return String representation of the allele */
  @Override
  String toString();
}
