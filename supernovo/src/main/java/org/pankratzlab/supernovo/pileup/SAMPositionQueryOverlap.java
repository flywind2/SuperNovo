package org.pankratzlab.supernovo.pileup;

import java.util.stream.Stream;
import org.pankratzlab.supernovo.GenomePosition;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;

public class SAMPositionQueryOverlap implements SAMPositionOverlap {

  private final Stream<SAMRecord> records;

  public SAMPositionQueryOverlap(SamReader samReader, GenomePosition position) {
    try (SAMRecordIterator iterator =
        samReader.queryOverlapping(
            position.getContig(), position.getPosition(), position.getPosition())) {
      records = iterator.stream();
    }
  }

  /* (non-Javadoc)
   * @see org.pankratzlab.supernovo.pileup.SAMPositionOverlap#getRecords()
   */
  @Override
  public Stream<SAMRecord> getRecords() {
    return records;
  }
}
