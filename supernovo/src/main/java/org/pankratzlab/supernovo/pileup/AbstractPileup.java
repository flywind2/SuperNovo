package org.pankratzlab.supernovo.pileup;

import java.util.Optional;

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
