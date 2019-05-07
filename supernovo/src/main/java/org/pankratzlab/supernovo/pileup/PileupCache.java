package org.pankratzlab.supernovo.pileup;

import java.util.Deque;
import java.util.Iterator;
import java.util.NoSuchElementException;
import java.util.stream.IntStream;
import java.util.stream.Stream;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.pankratzlab.supernovo.GenomePosition;
import com.google.common.cache.CacheBuilder;
import com.google.common.cache.CacheLoader;
import com.google.common.cache.ForwardingLoadingCache;
import com.google.common.cache.LoadingCache;
import com.google.common.collect.Lists;
import com.google.common.collect.Range;
import com.google.common.collect.RangeSet;
import com.google.common.collect.TreeRangeSet;
import htsjdk.samtools.QueryInterval;
import htsjdk.samtools.SamReader;

public class PileupCache extends ForwardingLoadingCache<GenomePosition, Pileup>
    implements AutoCloseable, Iterator<Pileup> {
  private static final Logger LOG = LogManager.getLogger(PileupCache.class);

  private final RangeSet<GenomePosition> intervals;

  private final IteratingPileupGenerator pileupGenerator;
  private final LoadingCache<GenomePosition, Pileup> cachedPileups;
  private final Deque<Pileup> readAheadPileups;
  private final SamReader samReaderAlt;
  private GenomePosition lastIterated;
  private long lastLogged = System.currentTimeMillis();

  public PileupCache(
      SamReader samReader, SamReader samReaderAlt, RangeSet<GenomePosition> intervals) {
    this.intervals = TreeRangeSet.create(intervals);
    pileupGenerator =
        new IteratingPileupGenerator(
            samReader.queryOverlapping(generateQueryIntervals(samReader, intervals)));
    this.samReaderAlt = samReaderAlt;
    cachedPileups = CacheBuilder.newBuilder().softValues().build(CacheLoader.from(this::load));
    readAheadPileups = Lists.newLinkedList();
  }

  private static QueryInterval[] generateQueryIntervals(
      SamReader samReader, RangeSet<GenomePosition> intervals) {
    final Stream<QueryInterval> qiStream;
    if (intervals.encloses(Range.all()))
      qiStream =
          IntStream.range(0, samReader.getFileHeader().getSequenceDictionary().size())
              .mapToObj(i -> new QueryInterval(i, 0, 0));
    else qiStream = intervals.asRanges().stream().map(r -> convertRangeToQI(samReader, r));
    return QueryInterval.optimizeIntervals(qiStream.sorted().toArray(QueryInterval[]::new));
  }

  private static QueryInterval convertRangeToQI(SamReader samReader, Range<GenomePosition> range) {
    GenomePosition start = range.lowerEndpoint();
    GenomePosition stop = range.upperEndpoint();
    String contig = start.getContig();
    if (!stop.getContig().equals(contig))
      throw new IllegalArgumentException("Illegal query range includes more than one contig");
    int seqIndex = samReader.getFileHeader().getSequenceIndex(contig);
    return new QueryInterval(seqIndex, start.getPosition(), stop.getPosition());
  }

  private Pileup load(GenomePosition loadPosition) {
    if (intervals.contains(loadPosition)) {
      if (loadPosition.compareTo(lastIterated) > 0 && hasNext()) {
        Pileup nextPileup;
        do {
          nextPileup = generateNext();
          readAheadPileups.addLast(nextPileup);
        } while (loadPosition.compareTo(nextPileup.getPosition()) > 0);
        return nextPileup;
      }
      return new Pileup(
          new SAMPositionQueryOverlap(samReaderAlt, loadPosition).getRecords(), loadPosition);
    }
    return new Pileup.Builder(loadPosition).build();
  }

  @Override
  public boolean hasNext() {
    return !readAheadPileups.isEmpty() || pileupGenerator.hasNext();
  }

  @Override
  public Pileup next() {
    if (!hasNext()) throw new NoSuchElementException();
    maybeLog();
    return updateLastIterated(iterateNext());
  }

  private void maybeLog() {
    if (System.currentTimeMillis() - lastLogged > 10000) {
      lastLogged = System.currentTimeMillis();
      LOG.info(
          "Last piled position: {}, current read-ahead position: {}",
          lastIterated,
          readAheadPileups.isEmpty() ? "None" : readAheadPileups.peekLast().getPosition());
    }
  }

  private Pileup updateLastIterated(final Pileup nextPileup) {
    GenomePosition newPosition = nextPileup.getPosition();
    // Remove any skipped over intervals from the targeted intervals
    if (lastIterated != null && newPosition.getContig().equals(lastIterated.getContig())) {
      String contig = newPosition.getContig();
      int start = lastIterated.getPosition() + 1;
      int stop = newPosition.getPosition() - 1;
      if (start <= stop)
        intervals.remove(
            Range.closed(new GenomePosition(contig, start), new GenomePosition(contig, stop)));
    }

    lastIterated = newPosition;
    return nextPileup;
  }

  private Pileup iterateNext() {
    if (readAheadPileups.isEmpty()) return generateNext();
    return readAheadPileups.removeFirst();
  }

  private Pileup generateNext() {
    Pileup nextPileup = pileupGenerator.next();
    put(nextPileup.getPosition(), nextPileup);
    return nextPileup;
  }

  @Override
  protected LoadingCache<GenomePosition, Pileup> delegate() {
    return cachedPileups;
  }

  @Override
  public void close() {
    pileupGenerator.close();
  }
}
