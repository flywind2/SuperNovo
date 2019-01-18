package org.pankratzlab.supernovo;

import java.util.function.Function;

public abstract class SampleMetric implements Metric {

  private final Function<Position, Pileup> pileupFunc;

  /** @param pileupFunc */
  public SampleMetric(Function<Position, Pileup> pileupFunc) {
    super();
    this.pileupFunc = pileupFunc;
  }

  protected Pileup getPileup(Position pos) {
    return pileupFunc.apply(pos);
  }
}
