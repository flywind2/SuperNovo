package org.pankratzlab.supernovo.pileup;

import java.util.Iterator;
import java.util.NavigableMap;
import java.util.NoSuchElementException;
import org.pankratzlab.supernovo.GenomePosition;
import com.google.common.base.Predicates;
import com.google.common.collect.Iterators;
import com.google.common.collect.Maps;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;

public class IteratingPileupGenerator implements Iterator<Pileup>, AutoCloseable {

  private final NavigableMap<GenomePosition, Pileup.Builder> pileupBuilders;
  private final Iterator<SAMRecord> mappedRecordsIter;
  private SAMRecord curRecord;
  private final Runnable onClose;

  public IteratingPileupGenerator(SamReader samReader) {
    @SuppressWarnings("resource") // Closed by onClose in close method
    final SAMRecordIterator samRecordIter =
        samReader.iterator().assertSorted(SAMFileHeader.SortOrder.coordinate);
    onClose = samRecordIter::close;
    mappedRecordsIter =
        Iterators.filter(samRecordIter, Predicates.not(SAMRecord::getReadUnmappedFlag));
    curRecord = mappedRecordsIter.next();
    pileupBuilders = Maps.newTreeMap();
  }

  public IteratingPileupGenerator(final SAMRecordIterator samRecordIter) {
    samRecordIter.assertSorted(SAMFileHeader.SortOrder.coordinate);
    onClose = samRecordIter::close;
    mappedRecordsIter =
        Iterators.filter(samRecordIter, Predicates.not(SAMRecord::getReadUnmappedFlag));
    curRecord = mappedRecordsIter.next();
    pileupBuilders = Maps.newTreeMap();
  }

  @Override
  public boolean hasNext() {
    return !pileupBuilders.isEmpty() || mappedRecordsIter.hasNext();
  }

  @Override
  public Pileup next() {
    if (!hasNext()) throw new NoSuchElementException();
    final GenomePosition targetGP;
    if (pileupBuilders.isEmpty()) {
      targetGP = new GenomePosition(curRecord.getContig(), curRecord.getAlignmentStart());
    } else {
      targetGP = pileupBuilders.firstKey();
    }
    final int targetPos = targetGP.getPosition();
    final String targetContig = targetGP.getContig();
    while (curRecord.getContig().equals(targetContig)
        && curRecord.getAlignmentStart() <= targetPos) {
      addRecord(curRecord);
      if (mappedRecordsIter.hasNext()) curRecord = mappedRecordsIter.next();
      else break;
    }
    Pileup targetPileup = pileupBuilders.remove(targetGP).build();
    pileupBuilders.headMap(targetGP).clear();
    return targetPileup;
  }

  private void addRecord(SAMRecord record) {
    final String contig = record.getContig();
    //    for (int i = 1; i <= record.getReadLength(); i++) {
    //      int refPos = record.getReferencePositionAtReadPosition(i);
    //      if (refPos != 0) {
    //        pileupBuilders
    //            .computeIfAbsent(new GenomePosition(contig, refPos), Pileup.Builder::new)
    //            .addRecord(record);
    //      }
    //    }
    for (int i = record.getAlignmentStart(); i <= record.getAlignmentEnd(); i++) {
      pileupBuilders
          .computeIfAbsent(new GenomePosition(contig, i), Pileup.Builder::new)
          .addRecord(record);
    }
  }

  @Override
  public void close() {
    onClose.run();
  }
}
