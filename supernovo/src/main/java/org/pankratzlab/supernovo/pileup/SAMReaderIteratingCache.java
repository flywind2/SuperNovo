package org.pankratzlab.supernovo.pileup;

import java.util.Optional;
import java.util.TreeMap;
import org.pankratzlab.supernovo.GenomePosition;
import com.google.common.collect.ImmutableList;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;

public class SAMReaderIteratingCache {

  private final SamReader reader;
  private SAMRecordIterator contigIter;
  private final TreeMap<Integer, ImmutableList.Builder<SAMRecord>> unpiledPositions;
  private String curContig;

  public class SAMPositionIterationOverlap implements SAMPositionOverlap {

    private final ImmutableList<SAMRecord> records;

    public SAMPositionIterationOverlap(GenomePosition position) {
      if (position.getContig() != curContig) {
        curContig = position.getContig();
        unpiledPositions.clear();
        contigIter = reader.queryContained(curContig, 0, Integer.MAX_VALUE);
      }
      int searchPos = position.getPosition();
      while (contigIter.hasNext()) {
        SAMRecord curRecord = contigIter.next();
        addRecord(curRecord);
        if (curRecord.getAlignmentStart() > searchPos) break;
      }
      records =
          Optional.ofNullable(unpiledPositions.remove(searchPos))
              .map(ImmutableList.Builder::build)
              .orElse(ImmutableList.of());
      unpiledPositions.headMap(searchPos).clear();
    }

    @Override
    public ImmutableList<SAMRecord> getRecords() {
      return records;
    }
  }

  /** @param reader */
  public SAMReaderIteratingCache(SamReader reader) {
    super();
    this.reader = reader;
    this.unpiledPositions = new TreeMap<>();
    this.curContig = null;
  }

  private void addRecord(SAMRecord record) {
    for (int i = record.getAlignmentStart(); i <= record.getAlignmentEnd(); i++) {
      unpiledPositions.computeIfAbsent(i, k -> ImmutableList.builder()).add(record);
    }
  }
}
