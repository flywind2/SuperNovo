package org.pankratzlab.supernovo.pileup;

import org.pankratzlab.supernovo.ReferencePosition;
import com.google.common.collect.ImmutableList;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;

public class SAMPositionQueryOverlap implements SAMPositionOverlap {

  private final ImmutableList<SAMRecord> records;

  public SAMPositionQueryOverlap(SamReader samReader, ReferencePosition position) {
    try (SAMRecordIterator iterator =
        samReader.queryOverlapping(
            position.getContig(), position.getPosition(), position.getPosition())) {
      records = iterator.stream().collect(ImmutableList.toImmutableList());
    }
  }

  /* (non-Javadoc)
   * @see org.pankratzlab.supernovo.pileup.SAMPositionOverlap#getRecords()
   */
  @Override
  public ImmutableList<SAMRecord> getRecords() {
    return records;
  }
}
