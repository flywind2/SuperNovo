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
}
