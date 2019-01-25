package org.pankratzlab.supernovo;

import com.google.common.cache.CacheBuilder;
import com.google.common.cache.CacheLoader;
import com.google.common.cache.LoadingCache;
import htsjdk.samtools.SAMRecord;

public class SNPAllele extends AbstractPileAllele {

  private static LoadingCache<Byte, SNPAllele> cache =
      CacheBuilder.newBuilder().softValues().build(CacheLoader.from(SNPAllele::new));

  private final byte base;

  /** @param base */
  private SNPAllele(byte base) {
    super(String.valueOf((char) base));
    this.base = base;
  }

  public static SNPAllele of(byte base) {
    return cache.getUnchecked(base);
  }

  @Override
  public boolean supported(SAMRecord record, int readPos) {
    return readPos != -1 && record.getReadBases()[readPos] == base;
  }

  /* (non-Javadoc)
   * @see java.lang.Object#hashCode()
   */
  @Override
  public int hashCode() {
    final int prime = 31;
    int result = super.hashCode();
    result = prime * result + base;
    return result;
  }

  /* (non-Javadoc)
   * @see java.lang.Object#equals(java.lang.Object)
   */
  @Override
  public boolean equals(Object obj) {
    if (this == obj) return true;
    if (!super.equals(obj)) return false;
    if (!(obj instanceof SNPAllele)) return false;
    SNPAllele other = (SNPAllele) obj;
    if (base != other.base) return false;
    return true;
  }
}
