package org.pankratzlab.supernovo.metrics;

import java.util.function.Function;
import org.pankratzlab.supernovo.Pileup;
import org.pankratzlab.supernovo.Position;

public abstract class SamplePileupMetric extends SampleMetric {

  private final Function<Position, Pileup> pileupFunc;

  /**
   * @param type
   * @param pileupFunc
   */
  public SamplePileupMetric(Type type, Function<Position, Pileup> pileupFunc) {
    super(type);
    this.pileupFunc = pileupFunc;
  }
}
