package org.pankratzlab.supernovo;

import com.google.common.collect.ImmutableList;
import htsjdk.samtools.SAMRecord;

public class InsertionAllele extends AbstractPileAllele {

  private class NonInsertionAllele extends AbstractPileAllele {

    private NonInsertionAllele() {
      super(preInsertionBase.toString());
    }

    @Override
    public boolean supported(SAMRecord record, int readPos) {
      return InsertionAllele.this.supportType(record, readPos).equals(Support.NO_INSERTION);
    }
  }

  private enum Support {
    INSERTION,
    NO_INSERTION,
    OTHER;
  }

  private final SNPAllele preInsertionBase;
  private final ImmutableList<Byte> insertedBases;
  private final NonInsertionAllele nonInsertionAllele;

  /**
   * @param preInsertionBase
   * @param insertedBases
   */
  public InsertionAllele(SNPAllele preInsertionBase, ImmutableList<Byte> insertedBases) {
    super(alleleString(preInsertionBase, insertedBases));
    this.preInsertionBase = preInsertionBase;
    this.insertedBases = insertedBases;
    this.nonInsertionAllele = new NonInsertionAllele();
  }

  @Override
  public boolean supported(final SAMRecord record, final int readPos) {
    return supportType(record, readPos).equals(Support.INSERTION);
  }

  private Support supportType(final SAMRecord record, final int readPos) {
    if (preInsertionBase.supported(record, readPos)) {
      boolean support = false;
      byte[] readBases = record.getReadBases();
      int offset = readPos + 1;
      for (int i = 0; i < insertedBases.size() && i + offset < readBases.length; i++) {
        if (readBases[i + offset] == insertedBases.get(i)) support = true;
        else {
          if (support) return Support.OTHER;
          return Support.NO_INSERTION;
        }
      }
      if (support) return Support.INSERTION;
      return Support.OTHER;
    }
    return Support.OTHER;
  }

  private static String alleleString(
      SNPAllele preInsertionBase, ImmutableList<Byte> insertedBases) {
    String preInsertion = preInsertionBase.toString();
    StringBuilder alleleBuilder = new StringBuilder(preInsertion.length() + insertedBases.size());
    alleleBuilder.append(preInsertion);
    insertedBases.forEach(b -> alleleBuilder.append((char) b.byteValue()));
    return alleleBuilder.toString();
  }

  /** @return a {@link PileAllele} that represents the alternative to this insertion */
  public NonInsertionAllele getNonInsertionAllele() {
    return nonInsertionAllele;
  }
}
