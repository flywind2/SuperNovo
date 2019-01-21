package org.pankratzlab.supernovo;

import java.util.Optional;
import org.pankratzlab.supernovo.metrics.Depth;

public abstract class AbstractPileup implements Pileup {

  private Optional<Depth> depth = Optional.empty();

  AbstractPileup() {}

  private Depth setDepth() {
    depth = Optional.of(new Depth(this));
    return depth.get();
  }

  @Override
  public Depth getDepth() {
    return depth.orElseGet(this::setDepth);
  }
}
