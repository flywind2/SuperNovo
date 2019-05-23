package org.pankratzlab.supernovo.pileup;

import java.io.Closeable;
import java.util.stream.Stream;
import org.pankratzlab.supernovo.GenomePosition;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;

public class SAMPositionQueryOverlap implements SAMPositionOverlap, Closeable {

  private final Stream<SAMRecord> records;
  private final SAMRecordIterator iterator;

  public SAMPositionQueryOverlap(SamReader samReader, GenomePosition position) {
    iterator =
        samReader.queryOverlapping(
            position.getContig(), position.getPosition(), position.getPosition());
    records = iterator.stream();
  }

  /* (non-Javadoc)
   * @see org.pankratzlab.supernovo.pileup.SAMPositionOverlap#getRecords()
   */
  @Override
  public Stream<SAMRecord> getRecords() {
    return records;
  }

  @Override
  public void close() {
    iterator.close();
  }
}
