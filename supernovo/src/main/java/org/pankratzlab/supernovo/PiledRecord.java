package org.pankratzlab.supernovo;

import htsjdk.samtools.SAMRecord;

public class PiledRecord {

  private final int id;
  private final int mq;
  private final boolean clipped;

  PiledRecord(SAMRecord samRecord) {
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