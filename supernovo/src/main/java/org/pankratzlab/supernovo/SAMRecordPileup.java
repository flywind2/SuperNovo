package org.pankratzlab.supernovo;

import com.google.common.collect.ImmutableMultimap;
import com.google.common.collect.ImmutableMultiset;
import com.google.common.collect.ImmutableSetMultimap;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;

public class SAMRecordPileup implements Pileup {

  public static class PiledRecord {

    private final int id;
    private final int mq;
    private final boolean clipped;

    private PiledRecord(SAMRecord samRecord) {
      id = samRecord.hashCode();
      mq = samRecord.getMappingQuality();
      clipped =
          samRecord.getUnclippedStart() != samRecord.getStart()
              || samRecord.getUnclippedEnd() != samRecord.getStart();
    }

    /*
     * (non-Javadoc)
     * @see java.lang.Object#hashCode()
     */
    @Override
    public int hashCode() {
      final int prime = 31;
      int result = 1;
      result = prime * result + (clipped ? 1231 : 1237);
      result = prime * result + id;
      result = prime * result + mq;
      return result;
    }

    /*
     * (non-Javadoc)
     * @see java.lang.Object#equals(java.lang.Object)
     */
    @Override
    public boolean equals(Object obj) {
      if (this == obj) return true;
      if (obj == null) return false;
      if (!(obj instanceof PiledRecord)) return false;
      PiledRecord other = (PiledRecord) obj;
      if (clipped != other.clipped) return false;
      if (id != other.id) return false;
      if (mq != other.mq) return false;
      return true;
    }
  }

  private final ImmutableSetMultimap<Byte, PiledRecord> basePiles;

  public SAMRecordPileup(SamReader samReader, Position position) {
    try (SAMRecordIterator iterator =
        samReader.queryOverlapping(
            position.getContig(), position.getPosition(), position.getPosition())) {
      ImmutableSetMultimap.Builder<Byte, PiledRecord> basePilesBuilder =
          ImmutableSetMultimap.builder();
      while (iterator.hasNext()) {
        SAMRecord samRecord = iterator.next();
        int readPos = samRecord.getReadPositionAtReferencePosition(position.getPosition()) - 1;
        if (readPos != -1)
          basePilesBuilder.put(samRecord.getReadBases()[readPos], new PiledRecord(samRecord));
      }
      basePiles = basePilesBuilder.build();
    }
  }

  @Override
  public ImmutableMultiset<Byte> getBaseCounts() {
    return basePiles.keys();
  }

  @Override
  public ImmutableMultimap<Byte, PiledRecord> getRecordsByBase() {
    return basePiles;
  }
}
