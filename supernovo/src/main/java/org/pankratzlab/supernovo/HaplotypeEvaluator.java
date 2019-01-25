package org.pankratzlab.supernovo;

import java.util.List;
import java.util.Set;
import org.pankratzlab.supernovo.pileup.Depth.Allele;
import org.pankratzlab.supernovo.pileup.Pileup;
import com.google.common.collect.ImmutableList;
import com.google.common.collect.Sets;

public class HaplotypeEvaluator {

  public static class Result {
    private final int otherVariants;
    private final int otherTriallelics;
    private final int otherBiallelics;
    private final int otherDeNovos;
    private final List<Double> concordances;
    /**
     * @param otherVariants
     * @param otherTriallelics
     * @param otherBiallelics
     * @param otherDeNovos
     * @param concordances
     */
    public Result(
        int otherVariants,
        int otherTriallelics,
        int otherBiallelics,
        int otherDeNovos,
        List<Double> concordances) {
      super();
      this.otherVariants = otherVariants;
      this.otherTriallelics = otherTriallelics;
      this.otherBiallelics = otherBiallelics;
      this.otherDeNovos = otherDeNovos;
      this.concordances = concordances;
    }
    /** @return the otherVariants */
    public int getOtherVariants() {
      return otherVariants;
    }
    /** @return the otherTriallelics */
    public int getOtherTriallelics() {
      return otherTriallelics;
    }
    /** @return the otherBiallelics */
    public int getOtherBiallelics() {
      return otherBiallelics;
    }
    /** @return the otherDeNovos */
    public int getOtherDeNovos() {
      return otherDeNovos;
    }
    /** @return the concordances */
    public List<Double> getConcordances() {
      return concordances;
    }
  }

  private static final int HAPLOTYPE_SEARCH_DISTANCE = 150;
  private static final double MIN_HAPLOTYPE_CONCORDANCE = 0.75;

  private final ReferencePosition pos;
  private final Pileup child;
  private final Pileup p1;
  private final Pileup p2;
  /**
   * @param child
   * @param p1
   * @param p2
   */
  public HaplotypeEvaluator(ReferencePosition pos, Pileup child, Pileup p1, Pileup p2) {
    super();
    this.pos = pos;
    this.child = child;
    this.p1 = p1;
    this.p2 = p2;
  }

  public Result haplotypeConcordance() {
    int startSearch = Integer.max(0, pos.getPosition() - HAPLOTYPE_SEARCH_DISTANCE);
    int stopSearch = pos.getPosition() + HAPLOTYPE_SEARCH_DISTANCE;

    int otherDenovos = 0;
    int otherTriallelics = 0;
    int otherBiallelics = 0;
    int otherVariants = 0;
    ImmutableList.Builder<Double> concordances = ImmutableList.builder();

    for (int searchPos = startSearch; searchPos < stopSearch; searchPos++) {
      if (searchPos == pos.getPosition()) continue;
      Position searchPosition = new Position(pos.getContig(), searchPos);
      Pileup searchPileup = searchPileup(child, searchPosition);
      if (TrioEvaluator.looksVariant(searchPileup.getDepth())) {
        otherVariants++;
        if (TrioEvaluator.moreThanTwoViableAlleles(searchPileup)) {
          otherTriallelics++;
        } else {
          otherBiallelics++;
          concordances.add(concordance(child, searchPileup));
          if (TrioEvaluator.looksDenovo(
              searchPileup, searchPileup(p1, searchPosition), searchPileup(p2, searchPosition))) {
            otherDenovos++;
          }
        }
      }
    }
    return new Result(
        otherVariants, otherTriallelics, otherBiallelics, otherDenovos, concordances.build());
  }

  private static double concordance(Pileup base, Pileup search) {
    Set<Integer> h1 = base.getDepth().allelicRecords(Allele.A1);
    Set<Integer> h2 = base.getDepth().allelicRecords(Allele.A2);

    Set<Integer> search1 = search.getDepth().allelicRecords(Allele.A1);
    Set<Integer> search2 = search.getDepth().allelicRecords(Allele.A2);

    double totalOverlap = search.getDepth().rawTotalDepth();
    int maxOverlap =
        Integer.max(
            Sets.intersection(search1, h1).size() + Sets.intersection(search2, h2).size(),
            Sets.intersection(search1, h2).size() + Sets.intersection(search2, h1).size());
    return maxOverlap / totalOverlap;
  }

  private static Pileup searchPileup(Pileup base, Position searchPos) {
    return new Pileup(base.getRecords(), searchPos);
  }
}
