package org.pankratzlab.supernovo.pileup;

import java.util.Deque;
import java.util.Iterator;
import java.util.NoSuchElementException;
import org.pankratzlab.supernovo.GenomePosition;
import com.google.common.cache.CacheBuilder;
import com.google.common.cache.CacheLoader;
import com.google.common.cache.ForwardingLoadingCache;
import com.google.common.cache.LoadingCache;
import com.google.common.collect.Lists;
import htsjdk.samtools.SamReader;

public class PileupCache extends ForwardingLoadingCache<GenomePosition, Pileup>
    implements AutoCloseable, Iterator<Pileup> {

  private final IteratingPileupGenerator pileupGenerator;
  private final LoadingCache<GenomePosition, Pileup> cachedPileups;
  private final Deque<Pileup> readAheadPileups;
  private final SamReader samReaderAlt;
  private GenomePosition lastIterated;

  public PileupCache(SamReader samReader, SamReader samReaderAlt) {
    pileupGenerator = new IteratingPileupGenerator(samReader);
    this.samReaderAlt = samReaderAlt;
    cachedPileups = CacheBuilder.newBuilder().softValues().build(CacheLoader.from(this::load));
    readAheadPileups = Lists.newLinkedList();
  }

  private Pileup load(GenomePosition loadPosition) {
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

  @Override
  public boolean hasNext() {
    return !readAheadPileups.isEmpty() || pileupGenerator.hasNext();
  }

  @Override
  public Pileup next() {
    if (!hasNext()) throw new NoSuchElementException();
    Pileup nextPileup = iterateNext();
    GenomePosition newPosition = nextPileup.getPosition();
    // Fill any voids with empty pileups to prevent attempted loading when queried
    if (lastIterated != null && newPosition.getContig().equals(lastIterated.getContig())) {
      for (int pos = lastIterated.getPosition() + 1; pos < newPosition.getPosition(); pos++) {
        GenomePosition emptyPos = new GenomePosition(newPosition.getContig(), pos);
        put(emptyPos, new Pileup.Builder(emptyPos).build());
      }
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
