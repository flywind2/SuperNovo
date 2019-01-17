package org.pankratzlab.supernovo;

import java.io.File;
import java.io.IOException;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.variant.vcf.VCFFileReader;
import picocli.CommandLine;
import picocli.CommandLine.Option;

public class App implements Runnable {

  @Option(names = {"--vcf", "-v"}, paramLabel = "VCF",
          description = "VCF with variants to query for de novo mutations", required = true)
  private File vcf;

  @Option(names = {"--childBam", "--bam"}, paramLabel = "BAM", description = "BAM of child",
          required = true)
  private File childBam;

  @Option(names = {"--childID", "--cID"}, paramLabel = "ID", description = "Sample ID of child",
          required = true)
  private String childID;

  @Option(names = {"--parent1Bam", "--p1Bam"}, paramLabel = "BAM", description = "BAM of parent 1",
          required = true)
  private File p1Bam;

  @Option(names = {"--parent1ID", "--p1ID"}, paramLabel = "ID",
          description = "Sample ID of parent 1", required = true)
  private String p1ID;

  @Option(names = {"--parent2Bam", "--p2Bam"}, paramLabel = "BAM", description = "BAM of parent 2",
          required = true)
  private File p2Bam;

  @Option(names = {"--parent2ID", "--p2ID"}, paramLabel = "ID",
          description = "Sample ID of parent 2", required = true)
  private String p2ID;

  @Option(names = {"--output", "-o"}, paramLabel = "FILE",
          description = "Output file for parsed de novo variants", required = true)
  private File output;

  public static void main(String[] args) {
    CommandLine.run(new App(), args);
  }

  @Override
  public void run() {
    SamReaderFactory srFactory = SamReaderFactory.make();
    try (SamReader child = srFactory.open(childBam);
         SamReader p1 = srFactory.open(p1Bam);
         SamReader p2 = srFactory.open(p2Bam);
         VCFFileReader vcfReader = new VCFFileReader(vcf)) {
      new TrioEvaluator(child, childID, p1, p1ID, p2, p2ID).reportDeNovos(vcfReader, output);
    } catch (IOException e) {
      e.printStackTrace();
    }

  }

}
