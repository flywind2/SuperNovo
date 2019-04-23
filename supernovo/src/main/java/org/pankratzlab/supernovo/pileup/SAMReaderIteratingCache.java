package org.pankratzlab.supernovo.pileup;

import java.util.Optional;
import java.util.TreeMap;
import org.pankratzlab.supernovo.GenomePosition;
import com.google.common.collect.ImmutableList;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamReader;

public class SAMReaderIteratingCache implements AutoCloseable {

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
        if (contigIter != null) contigIter.close();
        contigIter = reader.queryContained(curContig, 0, 0);
      }
      int searchPos = position.getPosition();
      unpiledPositions.headMap(searchPos).clear();
      while (contigIter.hasNext()) {
        SAMRecord curRecord = contigIter.next();
        if (curRecord.getAlignmentEnd() >= searchPos) addRecord(curRecord);
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

  @Override
  public void close() {
    if (contigIter != null) contigIter.close();
  }
}